import argparse
import warnings
import subprocess
from io import StringIO
import math
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
import json
import struct

# R = 1.98720425864083 * math.pow(10, -1)  # gas contant in dagcal⋅K−1⋅mol−1
R = 1.98720425864083 * math.pow(10, -3)  # gas contant in kcal⋅K−1⋅mol−1
T = 273.15  # 0 Celsius in K
MIN_RATE = 10 ** (
    -18
)  # 10 ** (-15) # depends on what float precission is used in treekin
DISSOCIATED_STATE_ENERGY = 0  # in dagcal⋅K−1⋅mol−1
ASSOCIATION_FACTOR = 1000000000000


class State:
    def __init__(
        self,
        index: int,
        base_pairs: (int, int),
        absorbing: bool,
        energy: float,
    ):
        """Initialize a State by index. Base pairs and name are derived based on interaction length."""
        self.index = index
        self.base_pairs = base_pairs
        self.absorbing = absorbing
        self.energy = energy

    def is_full_interaction(self, interaction_length):
        if self.base_pair[0] == 0 and self.base_pair[1] == interaction_length - 1:
            return True
        else:
            return False

    def name(self, one_based=False):
        prefix = "s"
        if self.absorbing:
            prefix = "a"
        if (
            one_based
            and type(self.base_pairs[0]) == int
            and type(self.base_pairs[1]) == int
        ):
            return f"{prefix}:{self.base_pairs[0]+1}:{self.base_pairs[1]+1}"
        else:
            return f"{prefix}:{self.base_pairs[0]}:{self.base_pairs[1]}"

    def __repr__(self):
        return f"State(index={self.index}, base_pairs={self.base_pairs}, absorbing={self.absorbing}, energy={self.energy})"


class MarkovProcess:
    def __init__(
        self,
        states_file,
        rates_file,
        absorbing_states=[],
        absorbin_full_interaction=False,
        dissociation_at=2,
        absorbing_dissociated_state=True,
        energy_penalty_absorbing_state=19,
        binary=True,
    ):
        (
            self._rates,
            self._states,
        ) = MarkovProcess.generate_treekin_rates_file(
            states_file,
            rates_file,
            energy_type="E",
            absorbing_states=absorbing_states,
            absorbin_full_interaction=absorbin_full_interaction,
            dissociation_at=dissociation_at,
            absorbing_dissociated_state=absorbing_dissociated_state,
            energy_penalty_absorbing_state=energy_penalty_absorbing_state,
            binary=binary,
            state_names_file=None,
            one_based_state_names=True,
        )
        # TODO set up class and write class functions that use static functions

    """
    def get_states(self):
        #TODO
        return states
    """

    @staticmethod
    def format_rates(rate):
        """Format single rate entry according to rates file format for treekin."""
        if rate == 0:
            rate_string = "{:>10}".format(rate)
        else:
            rate_string = "{:10.8f}".format(rate)
        return rate_string

    @staticmethod
    def get_rate(energy_i, energy_j, t=37):
        """Compute transition rate between two states from their free energies.
        energy_i: free energy of first state in kcal⋅mol−1
        energy_j: free energy of second state in kcal⋅mol−1
        t: temperature in Celsius
        """
        deltaE = max(energy_i, energy_j) - energy_i
        rate = math.exp(-deltaE / (R * (T + t)))
        return rate

    @staticmethod
    def check_rate(
        rate,
        states,
        min_rate=MIN_RATE,
        strict=False,
    ):
        """Print warning if rate to small."""
        if rate < min_rate:
            if strict:
                raise Exception(
                    f"Transition {states[0]} to {states[1]} has a to small rate.\n"
                    + f"rate: {rate}"
                )
            else:
                warnings.warn(
                    f"Transition {states[0]} to {states[1]} has small rate.\n"
                    + f"rate: {rate}"
                )
        return rate

    @staticmethod
    def get_interaction_length(states_file):
        """Return interaction length from directpaths states output."""
        df = pd.read_csv(states_file, sep="\t")
        return df["k"].max() + 1

    @staticmethod
    def two2oneD(k, l, interaction_length):
        """Get index of state that starts at bp k and ends at bp l."""
        i = (k + 1) * interaction_length  # all states (0:k,0:interaction_length)
        i -= (k + 1) * k / 2  # minus all states where l<k
        i -= interaction_length - l - 1  # minus states (k,l:interaction_length)
        return int(i) - 1  # -1 to get zero based index

    @staticmethod
    def generate_treekin_rates_file(
        states_file,
        rate_file,
        energy_type="E",
        absorbing_states=[],
        absorbin_full_interaction=False,
        dissociation_at=2,
        absorbing_dissociated_state=True,
        association_scaling_factor=ASSOCIATION_FACTOR,
        energy_penalty_absorbing_state=19,
        binary=True,
        state_names_file=None,
        one_based_state_names=True,
    ):
        """
        Create a rates input file for Treekin from RRIkinDP states output, enabling RNA-RNA interaction kinetics analysis.

        Parameters
        ----------
        states_file : str
            Path to the input states file generated by RRIkinDP. This file contains
            the states and their respective free energies.
        rate_file : str
            Output path for the computed rates file formatted for Treekin.
        energy_type : str, optional
            Specifies the energy column in the states file to use (e.g., 'E', 'ED1').
            Default is 'E'.
        absorbing_states : list of States, optional
            List of States to which absorbing states should be attached.
            (e.g., [State(0,(0,0),False,0.8), State(1,(0,1),False,-1)]).
        absorbin_full_interaction : bool, optional
            If True, attaches an absorbing state to the full interaction without requiring
            the specific index of the full interaction state. Default is False.
        dissociation_at : int or None, optional
            Specifies an interaction length threshold at which an interaction can directly
            dissociate. Use None to exclude a dissociated state. For example, set this to 2
            to skip a single base pair as a barrier to dissociation.
        absorbing_dissociated_state : bool, optional
            If True, attaches an absorbing state to the dissociated state. Default is False.
        association_scaling_factor : float/int, optional
            To somewhat account for the assumtion that association is slower than the folding
            process, association rates are mutiplied with a prefactor k = 1/association_scaling_factor.
        energy_penalty_absorbing_state : float, optional
            Energy difference in kcal/mol between the absorbing state and its connected state.
            This penalty should be large enough to ensure that the out-rate is negligible
            within the simulation time but also numerically stable. If `binary` is False,
            avoid values above 10 kcal/mol, as they may lead to rounding of rates to zero
            in the output file. Default is 19 kcal/mol.
        binary : bool, optional
            If True, the output rate file is written in binary format, allowing for higher
            precision. Default is True.
        state_names_file : str or None, optional
            Path to save a JSON file with state names. Defaults to `<rate_file>.json` if None.

        Returns
        -------
        tuple of (list of lists of float, list of str)
            A tuple containing:
            - matrix (list of lists): The rate matrix used for the Treekin analysis.
            - state_names (list of str): Names of the states, suitable for plotting and analysis.

        Notes
        -----
        The function generates a rate matrix based on provided states and configurations.
        States are read from the `states_file`, and transition rates are calculated based on
        their free energies. The matrix includes additional rows for any defined absorbing
        or dissociated states.

        - **Absorbing State Caution**: If `energy_penalty_absorbing_state` is too large in
        non-binary mode (over ~10 kcal/mol), resulting rates may be rounded to zero, making
        the output unusable. For long simulations, consider using high-precision Treekin
        settings (`--mlapack-method`).

        - the rate files first contains all regular states, than putative dissociated, absorbing
        dissociated, general absorbing states (including full interaction).

        Examples
        --------
        ```
        matrix, state_names = treekin_rates_from_RRIkinDP_states(
            states_file="input_states.tsv",
            rate_file="output_rates.txt",
            energy_type="E",
            absorbing_states=[12, 1, 24],
            absorbin_full_interaction=True,
            dissociation_at=2,
            absorbing_dissociated_state=True,
            energy_penalty_absorbing_state=15,
            binary=False
        )
        ```

        Todo
        ----
        - Modularize.
        - Current version avoids some dependencies. As eg pandas is anyway used in evaluation
        scripts, it might also be used here.
        """

        # throw exception if text base rate file is used association scaling factor is not compatible with non-binary rate files.
        if not binary and absorbing_dissociated_state:
            if association_scaling_factor > 1 / MIN_RATE:
                raise Exception(
                    "Association rates are to small to be saved in non binary rates file."
                    + " Please use a smaller association scaling factor or binary rate file."
                )

        # set minimum rate depending on wether a rate file format is binary
        if binary:
            min_rate = MIN_RATE
        else:
            min_rate = 0.00000001

        # Warning on absorbing state rates
        if (not binary) and (
            MarkovProcess.get_rate(0, energy_penalty_absorbing_state) < min_rate
        ):
            raise Exception(
                "The set energy penalty for absorbing states  leads to small rates that"
                + " can not be represented in non binary rate file. "
                + "The current energy penalty "
                + f"is {energy_penalty_absorbing_state} kcal/mol, corresponding "
                + f"to a rate of {MarkovProcess.get_rate(0,energy_penalty_absorbing_state)}.\n"
                + "Recomended energy penalties for absorbing states when using "
                + "non binary rate files are <= 10kcal/mol."
            )

        # set up data structures with states info
        states = []
        states_dict = {}
        interaction_length = 0
        with open(states_file, "r") as f:
            i = 0
            for line in f:
                if i == 0:
                    lables = line.strip().split("\t")
                    index_e = lables.index(energy_type)
                    index_i = lables.index("k")
                    index_j = lables.index("l")
                else:
                    data = line.strip().split("\t")
                    states.append(
                        State(
                            index=i - 1,
                            base_pairs=(
                                int(data[index_i]),
                                int(data[index_j]),
                            ),
                            energy=float(data[index_e]) / 100,
                            absorbing=False,
                        )
                    )
                    states_dict[
                        (
                            int(data[index_i]),
                            int(data[index_j]),
                        )
                    ] = {
                        "energy": float(data[index_e]) / 100,
                        "index": i - 1,
                    }
                    if int(data[index_j]) > interaction_length:
                        interaction_length = int(data[index_j])
                i = i + 1
        interaction_length += 1
        energies = [state.energy for state in states]

        # set up absorbing state
        ## convert to index based
        absorbing_states = [
            states_dict[state[0], state[1]]["index"] for state in absorbing_states
        ]

        ## add absorbing full interaction
        if absorbin_full_interaction:
            absorbing_states.append(states_dict[0, interaction_length - 1]["index"])
        absorbing_states = list(set(absorbing_states))
        absorbing_states.sort()

        # count states
        number_of_states = len(states) + len(absorbing_states)

        if dissociation_at is not None:
            number_of_states += 1
            # if absorbing_dissociated_state:
            # number_of_states += 1

        # build rate matrix
        matrix = []

        # add rows for each RRIkinDP state
        for k in range(len(states)):
            row = [0.0] * (len(states))

            current_i, current_j = states[k].base_pairs

            connected_states_ij = [
                (current_i, current_j - 1),
                (current_i, current_j + 1),
                (current_i - 1, current_j),
                (current_i + 1, current_j),
            ]
            connected_states = [
                states_dict[state]["index"]
                for state in connected_states_ij
                if state in states_dict.keys()
            ]
            for l in connected_states:
                row[l] = MarkovProcess.check_rate(
                    MarkovProcess.get_rate(energies[k], energies[l]),
                    (
                        states[k].name(),
                        states[l].name(),
                    ),
                    min_rate=min_rate,
                    strict=not (binary),
                )

            ## add column entries for dissociated state
            if dissociation_at is not None:
                col_entry = 0.0
                if current_j - current_i < dissociation_at:
                    col_entry = MarkovProcess.check_rate(
                        MarkovProcess.get_rate(
                            energies[k],
                            DISSOCIATED_STATE_ENERGY,
                        ),
                        (
                            states[k].name(),
                            "dissociated-state",
                        ),
                        min_rate=min_rate,
                        strict=not (binary),
                    )
                row.append(col_entry)

            ## add absorbing states column entries
            rates_to_absorbing = [0.0] * len(absorbing_states)
            if k in absorbing_states:

                index_absorbing = absorbing_states.index(k)
                rates_to_absorbing[index_absorbing] = MarkovProcess.check_rate(
                    MarkovProcess.get_rate(
                        energies[k],
                        energies[k] - energy_penalty_absorbing_state,
                    ),
                    (
                        states[k].name(),
                        "absorbing-state",
                    ),
                    min_rate=min_rate,
                    strict=not (binary),
                )
            row += rates_to_absorbing

            matrix.append(row)

        # add rate row for dissociated state
        if dissociation_at is not None:
            states.append(
                State(
                    index=len(states),
                    base_pairs=("d", "d"),
                    energy=DISSOCIATED_STATE_ENERGY,
                    absorbing=absorbing_dissociated_state,
                )
            )
            row = [0.0] * (len(states))
            for i in range(interaction_length):
                for j in range(i, i + dissociation_at):
                    if j >= interaction_length:
                        continue
                    k = states_dict[(i, j)]["index"]
                    association_scaling = 1
                    if absorbing_dissociated_state:
                        association_scaling = association_scaling_factor
                    row[k] = MarkovProcess.check_rate(
                        MarkovProcess.get_rate(
                            DISSOCIATED_STATE_ENERGY,
                            states[k].energy,
                        )
                        / association_scaling,
                        (
                            "dissociated-state",
                            states[k].name(),
                        ),
                        min_rate=min_rate,
                        strict=not (binary),
                    )

            row += [0.0] * len(absorbing_states)
            matrix.append(row)

        # add rate rows for absorbing states (except for absorbing dissociated state)
        for a in absorbing_states:
            states.append(
                State(
                    index=len(states),
                    base_pairs=states[a].base_pairs,
                    energy=states[a].energy - energy_penalty_absorbing_state,
                    absorbing=True,
                )
            )
            row = [0.0] * number_of_states
            row[a] = MarkovProcess.check_rate(
                MarkovProcess.get_rate(
                    states[a].energy - energy_penalty_absorbing_state,
                    states[a].energy,
                ),
                ("absorbing_state", a),
                min_rate=min_rate,
                strict=not (binary),
            )
            matrix.append(row)

        # write matrix to file
        if binary:
            # transpose matrix
            t_matrix = [
                [matrix[j][i] for j in range(len(matrix))]
                for i in range(len(matrix[0]))
            ]
            out = open(rate_file, "w+b")
            out.write(struct.pack("<i", number_of_states))
            for row in t_matrix:
                for e in row:
                    out.write(struct.pack("<d", e))
            out.close()
        else:
            out = open(rate_file, "w")
            for row in matrix:
                for e in row:
                    out.write(MarkovProcess.format_rates(e) + " ")
                out.write("\n")
            out.close()

        # generate and save state names
        state_names = [state.name(one_based=one_based_state_names) for state in states]
        if state_names_file is not None:
            with open(state_names_file, "w") as json_file:
                json.dump(state_names, json_file)

        return matrix, states

    @staticmethod
    def run_treekin(
        rate_file,
        initial_distribution,
        binary=True,
        treekin_executable="treekin",
        write_treekin_output_files=True,
        treekin_output_file=None,
        time_increment=1.02,
        sim_start_time=0.1,
        sim_end_time="1E8",
        temperature=37.0,
        verbose=False,
    ):
        """
        Run the treekin executable to compute RNA-RNA interaction dynamics.

        Parameters
        ----------
        rate_file : str
            Path to the rates input file in the specified format (binary or plain text).
        initial_distribution : list of lists, [[state1, state1_prob], [state2, state2_prob], ...]
            Initial state distribution.
        binary : bool, optional
            If True, the rates input file is in binary format, which offers higher precision (default is True).
        treekin_executable : str, optional
            Path to the treekin executable. Default assumes 'treekin' is available in the system PATH.
        write_treekin_output_files : bool, optional
            If True, writes the treekin output to a file specified by 'treekin_output_file' (default is True).
        treekin_output_file : str or None, optional
            Path to save the raw output from treekin. Required if 'write_treekin_output_files' is True.
        time_increment: float, optional
            Time scaling factor for logarithmic time scale (default 1.02).
        sim_start_time: float, optional
            Set simulation start time in internal units (default 0.1).
        sim_end_time: string or float, optional
            Set simulation stop time in internal units (default 1e+10).
        temperature: float, optional
            Set the simulation temperature in Celsius to temp (default 37.0).
        verbose : bool, optional
            If True, prints detailed output and errors from treekin execution (default is False).

        Returns
        -------
        None
            Executes treekin, returns treekin output and saves it to output file if specified.

        """

        treekin_args = [
            treekin_executable,
            "-m",
            "I",
            "--ratesfile",
            rate_file,
            "--p0",
            " ".join(
                [
                    f"{state.index+1}={probability}"
                    for state, probability in initial_distribution
                ]
            ),
            "--mlapack-method",
            "DD",
            "-T",
            f"{temperature}",
            "--tinc",
            f"{time_increment}",
            "--t0",
            f"{sim_start_time}",
            "--t8",
            f"{sim_end_time}",
            # "MPFR",  # "LD", "QD",  "DD", "DOUBLE", "GMP", "MPFR", "FLOAT128"
            # "--mlapack-precision",  # necessary if "GMP", "MPFR"
            # "128",
        ]
        if binary:
            treekin_args.append("--bin")

        treekin_process = subprocess.Popen(
            treekin_args,
            # stdin=cat_process.stdout,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        (
            stdout_data,
            stderr_data,
        ) = treekin_process.communicate()
        stderr_data = stderr_data.decode()
        treekin_output = StringIO(stdout_data.decode())

        if write_treekin_output_files:
            if treekin_output_file is not None:
                with open(treekin_output_file, "w") as out_handle:
                    out_handle.write(treekin_output.read())

        if verbose:
            print(stderr_data)

        return treekin_output

    @staticmethod
    # read state names from json file
    def state_names_from_json(json_file):
        """Return content of a json file."""
        with open(json_file) as f:
            state_names = json.load(f)
        return state_names

    @staticmethod
    # get name of full interaction state from states_names
    def get_full_interaction_state(states):
        """Get full interaction name via get max interaction base pair index."""
        max_bp = max(
            [
                int(state.base_pairs[1])
                for state in states
                if type(state.base_pairs[1]) in [int, float]
            ]
        )
        full_state = [
            state
            for state in states
            if state.base_pairs == (0, max_bp) and state.absorbing == False
        ][0]
        return full_state

    @staticmethod
    # get kinetic features from treekin output
    def get_treekin_features(
        treekin_out_file,
        states,
        target_states=[],
        features_json=None,
        eval_full_interaction=False,
        eval_dissociated_state=False,
        eval_sum_absorbing=False,
    ):
        """
        Extract kinetic features from the Treekin output, including the times to
        reach specific population thresholds and state populations at key time points.

        Parameters
        ----------
        treekin_out_file : str
            Path to the Treekin output file containing population data for each state
            over time.
        states : list of States
            List of State corresponding to the columns in the Treekin output file.
        target_states : list of str, optional
            List of states (identified by name) to evaluate. Defaults to an empty list,
            which evaluates no specific states.
        features_json : str or None, optional
            Path to save the resulting features dictionary in JSON format. If None, no
            JSON file is saved.
        eval_full_interaction : bool, optional
            If True, includes the state representing the full interaction for analysis.
            Default is False.
        eval_dissociated_state : bool, optional
            If True, includes the dissociated state for evaluation. Default is False.
        eval_sum_absorbing : bool, optional
            If True, includes the sum of all dissociated states for evaluation.
            Default is False.

        Returns
        -------
        dict
            A dictionary where each key corresponds to a feature
            name (e.g., `state_t99`, `state_p1E10`)for each evaluated state. Values
            represent time points or population fractions for each feature.

        Notes
        -----
        - The function identifies specific time points when the population for each
          state reaches 50% and 99% (`_t50` and `_t99`).
        - Population fractions at fixed time points are also
          extracted (`_p1E10`, `_p1E5`, etc.).
        - If `features_json` is provided, results are saved as a JSON file.

        Example
        -------
        ```
        features = get_treekin_features(
            "output_treekin.txt",
            state_names=["state_0", "state_1", "a:2:5", a:f:f"],
            target_states=["state_0", "d:d", "a:2:5"],
            features_json="features.json",
            eval_full_interaction=True,
            eval_dissociated_state=True
        )
        ```

        Todo
        ----
        - Make time points and population fractions into parameters.
        """

        state_names = [state.name() for state in states]

        df = pd.read_csv(
            treekin_out_file,
            # index_col=0,
            header=None,
            names=["time"] + state_names + ["empty"],
            sep=" ",
            comment="#",
        )

        if eval_dissociated_state:
            target_states.append("s:d:d")
            if "a:d:d" in state_names:
                target_states.append("a:d:d")

        if eval_full_interaction:
            target_states.append("s:f:f")
            full_state = MarkovProcess.get_full_interaction_state(states)
            df["s:f:f"] = df[full_state.name()]
            if full_state.name().replace("s", "a") in state_names:
                target_states.append("a:f:f")
                df["a:f:f"] = df[full_state.name().replace("s", "a")]

        if eval_sum_absorbing:
            absorbing_states = [state.name() for state in states if state.absorbing]
            df["a:a:a"] = df[absorbing_states].sum(axis=1)
            target_states.append("a:a:a")

        target_states = list(set(target_states))

        times = df["time"].to_list()

        data_dict = {}
        for state in target_states:

            # if state does not exist return nan
            if state not in df.columns:
                data_dict[f"{state}_t99"] = "nan"
                data_dict[f"{state}_t50"] = "nan"
                data_dict[f"{state}_logt99"] = "nan"
                data_dict[f"{state}_logt50"] = "nan"
                data_dict[f"{state}_p1E8"] = "nan"
                data_dict[f"{state}_p1E5"] = "nan"
                data_dict[f"{state}_p1E3"] = "nan"
                data_dict[f"{state}_p1E2"] = "nan"
                data_dict[f"{state}_p1E1"] = "nan"
                continue

            # get probability of state per simulation time step
            probs = df[state].to_list()

            # get time point at which >= 99 percent of population are in the state the first time
            if probs[-1] < 0.99:
                data_dict[f"{state}_t99"] = "nan"
                data_dict[f"{state}_logt99"] = "nan"
            else:
                step_99_absorbed = next(x for x, val in enumerate(probs) if val >= 0.99)
                data_dict[f"{state}_t99"] = times[step_99_absorbed]
                data_dict[f"{state}_logt99"] = math.log10(times[step_99_absorbed])

            # get time point at which >= 50 percent of population are in the state the first time
            if probs[-1] < 0.50:
                data_dict[f"{state}_t50"] = "nan"
                data_dict[f"{state}_logt50"] = "nan"
            else:
                step_50_absorbed = next(x for x, val in enumerate(probs) if val >= 0.50)
                data_dict[f"{state}_t50"] = times[step_50_absorbed]
                data_dict[f"{state}_logt50"] = math.log10(times[step_50_absorbed])

            # get the population fraction (state probability) at the given time points
            step_1E8 = next(x for x, val in enumerate(times) if val >= 100000000)
            step_1E5 = next(x for x, val in enumerate(times) if val > 100000)
            step_1E3 = next(x for x, val in enumerate(times) if val > 1000)
            step_1E2 = next(x for x, val in enumerate(times) if val > 100)
            step_1E1 = next(x for x, val in enumerate(times) if val > 10)
            data_dict[f"{state}_p1E8"] = probs[step_1E8]
            data_dict[f"{state}_p1E5"] = probs[step_1E5]
            data_dict[f"{state}_p1E3"] = probs[step_1E3]
            data_dict[f"{state}_p1E2"] = probs[step_1E2]
            data_dict[f"{state}_p1E1"] = probs[step_1E1]

        # save features
        if features_json is not None:
            with open(features_json, "w") as json_file:
                json.dump(data_dict, json_file)

        return data_dict

    @staticmethod
    def plot_E_mean(
        treekin_out_file,
        states,
        figure_path=None,
        figsize=(2.7, 1.4),
        x_lim=(None, None),
        y_lim=(None, None),
        title=None,
        enable_tex_fonts=True,
    ):

        times, energies = MarkovProcess.get_E_mean(treekin_out_file, states)

        # set up fonts
        font_family = "sans-serif"
        if enable_tex_fonts:
            font_family = "serif"
        font_params = {
            "font.family": font_family,  # use serif/main font for text elements
            "font.size": 8,
            "text.usetex": enable_tex_fonts,  # use inline math for ticks
            "pgf.rcfonts": False,  # don't setup fonts from rc parameters
            # "pgf.preamble": [
            # "\\usepackage{units}",
            # "\\usepackage{metalogo}",
            # "\\usepackage{unicode-math}",  # unicode math setup
            # r"\setmathfont{xits-math.otf}",
            # r"\setmainfont{DejaVu Serif}", # serif font via preamble
            #    ]
        }

        plt.rcParams.update(font_params)

        # set figure size
        f, ax = plt.subplots(figsize=figsize, layout="constrained")

        # set axis labels
        if enable_tex_fonts:
            ax.set_ylabel("$\hat{E}$ (kcal/mol)")
        else:
            ax.set_ylabel("E_mean (kcal/mol)")

        ax.set_xlabel("Time (a.u.)")

        # set time axis to log scale
        ax.set_xscale("log")

        # plot E_mean
        ax.plot(times, energies)

        # set plot range
        ax.set_xlim(x_lim[0], x_lim[1])
        ax.set_ylim(y_lim[0], y_lim[1])

        # include title
        if title is not None:
            ax.set_title(title)

        # save figure
        if figure_path is not None:
            f.savefig(figure_path, bbox_inches="tight")

    @staticmethod
    def get_E_mean(
        treekin_out_file,
        states,
    ):
        state_names = [state.name() for state in states]
        df = pd.read_csv(
            treekin_out_file,
            index_col=0,
            header=None,
            names=["time"] + state_names + ["empty"],
            sep=" ",
            comment="#",
        )
        df.drop(columns=["empty"], inplace=True)
        for state in states:
            if state.absorbing and state.name() != "a:d:d":
                non_absorbing_state = [
                    na_state
                    for na_state in states
                    if (
                        na_state.base_pairs == state.base_pairs
                        and na_state.absorbing == False
                    )
                ][0]
                df[state.name()] = df[state.name()] * (non_absorbing_state.energy)
            else:
                df[state.name()] = df[state.name()] * state.energy
        df = df.copy()  # for defragmentation
        df["E_mean"] = df.sum(axis=1)
        # f, ax = plt.subplots(figsize=(40,40), layout="constrained")
        # sn.heatmap(df[state_names+['E_mean']], ax = ax)
        # f.savefig("test.pdf", bbox_inches="tight")
        return (
            df.index.to_list(),
            df["E_mean"].to_list(),
        )

    @staticmethod
    def get_E_mean_features(
        treekin_out_file,
        states,
        eval_times=[
            1,
            10,
            100,
            1000,
            10000,
            100000,
            1000000,
            100000000,
        ],
    ):
        times, energies = MarkovProcess.get_E_mean(treekin_out_file, states)
        features = {}
        for time in eval_times:
            if time > times[-1]:
                features[f"E_mean({time:.0E})"] = "nan"
            else:
                time_index = next(
                    i for i, val in enumerate(times) if val >= float(time)
                )
                features[f"E_mean({time:.1E})"] = energies[time_index]
        return features

    @staticmethod
    def plot_treekin(
        treekin_output,
        treekin_plot,
        states=None,
        labels=True,
        label_cutoff_fraction=0.1,
        figsize=(7, 4),
        x_lim=(None, None),
        y_lim=(-0.05, 1.05),
        title=None,
        enable_tex_fonts=True,
        one_based_state_names=True,
    ):
        """
        Plot the Treekin output, showing state probabilities over time.

        Parameters
        ----------
        treekin_output : str
            Path to the Treekin output file to be plotted. This file should contain
            the population data for each state over time.
        treekin_plot : str
            Output path for the saved plot file. The file format is inferred from the
            file extension, e.g., '.png', '.pdf', etc.
        state : list of States or None, optional
            List of States corresponding to the treekin output columns used to get state
            labels. If None, default numbering will be used.
        labels : bool, optional
            If True, labels for individual states will be displayed on the plot. Labels
            will only appear for states with a maximum population exceeding the
            `label_cutoff_fraction`. Default is True.
        label_cutoff_fraction : float, optional
            Minimum population fraction (0-1) that a state must reach at least once
            to be labeled in the plot. States with populations below this threshold
            will not be labeled, even if `labels` is True. Default is 0.1.
        figsize : tuple of float, optional
            Size of the figure in inches, given as (width, height). Default is (7, 4).
        x_lim : tuple of float or None, optional
            Plot range on the x-axis. Set to (None, None) to allow automatic limits
            based on data. Default is (None, None).
        y_lim : tuple of float, optional
            Plot range on the y-axis, usually set between -0.05 and 1.05 for population
            values between 0 and 1 with a small margin. Default is (-0.05, 1.05).
        title: string, optional
            Title to be included into the plot.
        enable_tex_fonts: bool, optional
            Wether to load predefined pgf preamble. Default is False.
        one_based_state_names: bool, optional
            Use one based base pair indices instead of zero based base pair indeices
            to label states.
        Returns
        -------
        None
            Saves the plot to the specified `treekin_plot` path.

        Todo
        ----
        - Consider adding labels as text annotations instead of markers for improved
          readability.
        - Evaluate using external libraries like `matplotlib-label-lines` to automate
          label positioning.

        Notes
        -----
        - **Logarithmic Time Axis**: The x-axis is displayed on a logarithmic scale to
          capture population dynamics over time.
        - **Labeling States**: If `labels` is enabled, states with populations that
          peak above `label_cutoff_fraction` are labeled on the plot. Absorbing states
          are marked with a
        distinct edge color (gray), and state names are derived from `state_names` if provided.
        - **File Format**: The plot's file format is determined by the file extension of
        `treekin_plot`. Ensure the extension matches the desired format (e.g., `.png`, `.pdf`).
        - **Tex support**: Requires a a working LaTeX installation.

        Example
        -------
        ```
        plot_treekin(
            treekin_output="treekin_output.txt",
            treekin_plot="state_populations.pdf",
            labels=True,
            label_cutoff_fraction=0.1,
            figsize=(8, 5),
            x_lim=(1, 1e6),
            y_lim=(0, 1)
        )
        ```
        """

        # set up fonts
        font_family = "sans-serif"
        if enable_tex_fonts:
            font_family = "serif"
        font_params = {
            "font.family": font_family,  # use serif/main font for text elements
            "font.size": 8,
            "text.usetex": enable_tex_fonts,  # use inline math for ticks
            "pgf.rcfonts": False,  # don't setup fonts from rc parameters
            # "pgf.preamble": [
            # "\\usepackage{units}",
            # "\\usepackage{metalogo}",
            # "\\usepackage{unicode-math}",  # unicode math setup
            # r"\setmathfont{xits-math.otf}",
            # r"\setmainfont{DejaVu Serif}", # serif font via preamble
            #    ]
        }

        plt.rcParams.update(font_params)

        # read treekin output file
        if states is not None:
            state_names = [
                state.name(one_based=one_based_state_names) for state in states
            ]
            df = pd.read_csv(
                treekin_output,
                index_col=0,
                header=None,
                names=state_names + ["nan"],
                sep=" ",
                comment="#",
            )

        else:
            df = pd.read_csv(
                treekin_output,
                index_col=0,
                header=None,
                sep=" ",
                comment="#",
            )

        df = df.iloc[:, :-1]  # remove empty column (tailing spaces in input)

        # set figure size
        f, ax = plt.subplots(figsize=figsize, layout="constrained")

        # set axis labels
        ax.set_ylabel("Population")
        ax.set_xlabel("Time (a.u.)")

        # set y_lim
        min_time = df.index.to_list()[0]
        max_time = df.index.to_list()[-1]

        if x_lim[0] is None:
            x_lim = (min_time / 4, x_lim[1])
        if x_lim[1] is None:
            x_lim = (x_lim[0], max_time * 4)

        # set min and max label x coordinates for markers
        # Todo: set based on marker size, figsize and x_lim
        min_x_label = x_lim[0] * 4
        max_x_label = x_lim[1] / 4

        # call plot function for each state
        for col in df.columns:

            # plot state population
            p = sn.lineplot(
                x=df.index,
                y=df[col],
                ax=ax,
                legend=False,
            )
            color = p.get_lines()[-1].get_color()
            max_population = df[col].max()

            # continue if no lable should be plotted
            if not labels:
                continue

            # check if state passes population cutoff to get labled
            if max_population < label_cutoff_fraction:
                continue
            max_population_time = df[[col]].idxmax()

            # set label position at max_population
            # marker_y = max_population
            marker_x = max_population_time.iloc[0]

            # reset marker y  coordinates if outside or close to y-limits
            if marker_x < min_x_label:
                marker_x = next(t for t in df.index.to_list() if t > min_x_label)
            elif marker_x > max_x_label:
                marker_x = next(
                    t for t in reversed(df.index.to_list()) if t < max_x_label
                )
            marker_y = df.at[marker_x, col]

            # set marker text and marker edge color
            edge_col = "white"
            text = ":".join(col.split(":")[1:])
            background_marker_edge_width = 0

            ## mark absorbing states with grey marker edge and remove leading "a:" in state name
            if col.startswith("a"):
                edge_col = "grey"
                background_marker_edge_width = 0.7

            # set text size in state labels
            text_markersize = 13
            if len(text) < 4:
                text_markersize = 11
                if len(text) < 2:
                    text_markersize = 8

            # plot state label background
            plt.plot(
                marker_x,
                marker_y,
                marker="o",
                color="white",
                markeredgecolor=edge_col,
                markeredgewidth=background_marker_edge_width,
                alpha=0.65,
                markersize=15,
            )

            # plot state label text
            plt.plot(
                marker_x,
                marker_y,
                marker="$%s$" % text,
                color=color,
                linewidth=0.5,
                markersize=text_markersize,
                markeredgecolor=color,
                markeredgewidth=0.6,
            )

        # set time axis to logarithmic scale
        ax.set_xscale("log")

        # set plot range
        ax.set_xlim(x_lim[0], x_lim[1])
        ax.set_ylim(y_lim[0], y_lim[1])

        # set title
        if title is not None:
            ax.set_title(title)

        # save figure
        f.savefig(treekin_plot, bbox_inches="tight")
        plt.close(f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute RNA-RNA interaction formation dynamics by solving the master equation.",
        # allow_abbrev=False,
    )

    # Required arguments
    parser.add_argument(
        "-s",
        "--states",
        help="Path to input states file generated by RRIkinDP.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-i",
        "--initial",
        help="Initial state with 100% population in the simulation defined by index of first and last base pair e.g., '-i 2:4' (zero based indices).",
        type=lambda x: tuple(map(int, x.split(":"))),
        required=True,
    )
    parser.add_argument(
        "-r",
        "--rates",
        help="Path to save generated rates file , formatted for treekin.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-p",
        "--probs",
        help="Path to save treekin's state probabilities output.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-f",
        "--figure",
        help="Path to save the plot of state probabilities over time.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-o",
        "--output_summary",
        help="Path to save the kinetic features to (json).",
        type=str,
        required=True,
    )

    # Absorbing state settings
    parser.add_argument(
        "--absorbing_states",
        help="List of absorbing states in base pair format, i.e. represented by the zero based index of there first and last base pair. Example: '14:17 5:7, 0:6'.",
        type=lambda s: [tuple(map(int, state.split(":"))) for state in s.split()],
        default=[],
    )
    parser.add_argument(
        "--absorbing_full_interaction",
        help="Attach an absorbing state to the full interaction. (True/False).",
        action="store_true",
    )
    parser.add_argument(
        "--dissociation_at",
        help="Interaction length threshold for direct dissociation (None for no dissociated state).",
        type=int,
        default=2,
    )
    parser.add_argument(
        "--association_scaling",
        help="Prefactor^-1 for association rates.",
        type=int,
        default=1000000000000,
    )
    parser.add_argument(
        "--non_absorbing_dissociated_state",
        help="Attach an absorbing state to the dissociated state. (True/False).",
        action="store_true",
    )
    parser.add_argument(
        "--energy_penalty_absorbing_state",
        help="Energy penalty in kcal/mol for the absorbing state relative to connected state. Recommended values are <= 10 kcal/mol for non-binary outputs.",
        type=float,
        default=19,
    )

    # State names output
    parser.add_argument(
        "--state_names_file",
        help="Path to save JSON file with state names. Defaults to <rate_file>.json.",
        type=str,
        default=None,
    )

    # Treekin solver options
    parser.add_argument(
        "--treekin_executable",
        help="Path to the treekin executable (default: 'treekin').",
        type=str,
        default="treekin",
    )
    parser.add_argument(
        "--human_readable_rates",
        help="Output rate file as text file instead of binary file (lower precision). (True/False).",
        action="store_true",
    )

    parser.add_argument(
        "--treekin_verbose",
        help="Print additional details from treekin's execution. (True/False).",
        action="store_true",
    )

    # Summary settings
    parser.add_argument(
        "--target_states",
        help="Additional states that should be included in the kinetic features summary, in addition to the full interaction and dissaociated state. Example: '0:3,6:12'.",
        type=lambda x: list(map(int, x.split(","))),
        default=[],
    )

    # Plotting options
    parser.add_argument(
        "--plot_no_labels",
        help="Display labels for each state in the plot. (False/True).",
        action="store_false",
    )
    parser.add_argument(
        "--plot_label_cutoff",
        help="Minimum population fraction required to label a state in the plot.",
        type=float,
        default=0.1,
    )
    parser.add_argument(
        "--figsize",
        help="Size of the plot in inches, e.g., '--figsize 7 4'.",
        type=lambda x: tuple(map(float, x.split())),
        default=(7, 4),
    )
    parser.add_argument(
        "--plot_x_lim",
        help="Plot range on x-axis as a tuple. (default: None, None).",
        type=lambda x: tuple(
            map(
                lambda v: None if v == None else float(v),
                x.split(),
            )
        ),
        default=(None, None),
    )
    parser.add_argument(
        "--plot_y_lim",
        help="Plot range on y-axis as a tuple (default: -0.05 1.05).",
        type=lambda x: tuple(map(float, x.split())),
        default=(-0.05, 1.05),
    )

    parser.add_argument(
        "--plot_title",
        help="Title to be shown in popluation probability plot (default: None).",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--E_mean_plot_path",
        help="Path to save E mean plot to (default: None).",
        type=str,
        default=None,
    )

    args = parser.parse_args()

    # processing args
    ## use binary rate file?
    binary = not (args.human_readable_rates)

    ## make disscociated state absorbing?
    absorbing_dissociated_state = True
    if args.non_absorbing_dissociated_state:
        absorbing_dissociated_state = False

    # Catch to large absorbing state penalties when using non binary rate files
    if not binary and args.energy_penalty_absorbing_state > 10:
        raise Exception(
            "The set energy penalty for absorbing states  leads to small rates that"
            + " can not be represented in non binary rate file. "
            + "The current energy penalty "
            + f"is {args.energy_penalty_absorbing_state} kcal/mol, corresponding "
            + f"to a rate of {MarkovProcess.get_rate(0,args.energy_penalty_absorbing_state)}.\n"
            + "Recomended energy penalties for absorbing states when using "
            + "non binary rate files are <= 10kcal/mol."
        )

    # Catch to large absorbing state penalties when using non binary rate files
    if not binary and absorbing_dissociated_state:
        if args.association_scaling > 100000000:
            raise Exception(
                "The set association scaling factor leads to rates that are to small"
                + " to be represented in a non-binary (text format) rates file. Either\n"
                + " - switch to binary rate files (no --human_readable_rates flag),\n"
                + " - use a smaller association scaling factor (eg --association_scaling 1 or 1000) or\n"
                + " - make the dissociated state non absorbing (--non_absorbing_dissociated_state)."
            )
        else:
            warnings.warn(
                "The set association scaling factor might lead to rates that are to small"
                + " to be represented in a non-binary (text format) rates file. If warnings/errors about to"
                + " small association rates appear (rates from dissociated state to other state), either\n"
                + " - switch to binary rate files (no --human_readable_rates flag),\n"
                + " - use a smaller association scaling factor (eg --association_scaling 1 or 1000) or\n"
                + " - make the dissociated state non absorbing (--non_absorbing_dissociated_state)."
            )

    # Generate rate matrix
    (
        matrix,
        states,
    ) = MarkovProcess.generate_treekin_rates_file(
        states_file=args.states,
        rate_file=args.rates,
        absorbing_states=args.absorbing_states,
        absorbin_full_interaction=args.absorbing_full_interaction,
        dissociation_at=args.dissociation_at,
        absorbing_dissociated_state=absorbing_dissociated_state,
        association_scaling_factor=args.association_scaling,
        energy_penalty_absorbing_state=args.energy_penalty_absorbing_state,
        binary=binary,
        state_names_file=args.state_names_file,
    )

    # Preprocess treekin input
    ## interaction length
    interaction_length = MarkovProcess.get_interaction_length(args.states)
    ## index of initial state
    initial_k, initial_l = args.initial
    initial_state_index = MarkovProcess.two2oneD(
        initial_k, initial_l, interaction_length
    )

    # Run treekin
    MarkovProcess.run_treekin(
        rate_file=args.rates,
        initial_distribution=[
            [
                State(
                    index=initial_state_index,
                    base_pairs=(
                        initial_k,
                        initial_l,
                    ),
                    absorbing=False,
                    energy=None,
                ),
                1,
            ]
        ],
        binary=binary,
        treekin_executable=args.treekin_executable,
        write_treekin_output_files=True,
        treekin_output_file=args.probs,
        verbose=args.treekin_verbose,
    )

    # Summarize dynamic features
    features = MarkovProcess.get_treekin_features(
        treekin_out_file=args.probs,
        states=states,
        target_states=args.target_states,
        features_json=args.output_summary,
        eval_full_interaction=True,
        eval_dissociated_state=True,
        eval_sum_absorbing=True,
    )

    e_features = MarkovProcess.get_E_mean_features(args.probs, states)

    print("Features form Markov process simulation:")
    for key, value in features.items():
        print(f"{key}: {value}")
    for key, value in e_features.items():
        print(f"{key}: {value}")

    if args.E_mean_plot_path is not None:
        MarkovProcess.plot_E_mean(
            args.probs,
            states,
            figure_path=args.E_mean_plot_path,
            enable_tex_fonts=False,
        )

    # Plot state probabilities
    MarkovProcess.plot_treekin(
        treekin_output=args.probs,
        treekin_plot=args.figure,
        states=states,
        labels=args.plot_no_labels,
        label_cutoff_fraction=args.plot_label_cutoff,
        figsize=args.figsize,
        x_lim=args.plot_x_lim,
        y_lim=args.plot_y_lim,
        title=args.plot_title,
        one_based_state_names=True,
    )
