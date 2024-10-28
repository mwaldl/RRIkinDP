import argparse
import warnings
import subprocess
from io import StringIO
import math
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sn
import json

R = 1.98720425864083 * math.pow(
    10, -3
)  # gas contant in dagcal⋅K−1⋅mol−1
T = 273.15  # 0 Celsius in K
MIN_RATE = 10 ** (-14)  # TODO: unit
DISSOCIATED_STATE_ENERGY = (
    0  # in dagcal⋅K−1⋅mol−1
)


def format_rates(rate):
    """Format single rate entry according to rates file format for treekin."""
    if rate == 0:
        rate_string = "{:>10}".format(rate)
    else:
        rate_string = "{:10.8f}".format(rate)
    return rate_string


def get_rate(energy_i, energy_j, t=37):
    """Compute transition rate between two states from their free energies.
    energy_i: free energy of first state in dagcal⋅mol−1
    energy_j: free energy of second state in dagcal⋅mol−1
    t: temperature in Celsius
    """
    deltaE = max(energy_i, energy_j) - energy_i
    rate = math.exp(-deltaE / (R * (T + t)))
    return rate


def check_rate(rate, states, min_rate=0.00000001):
    """Print warning if rate to small."""
    if rate < min_rate:
        warnings.warn(
            f"Transition {states[0]} to {states[1]} hast to small rate .\n"
            + f"rate: {rate}"
        )
    return rate


def get_interaction_length(states_file):
    """Return interaction length from directpaths states output."""
    df = pd.read_csv(states_file, sep="\t")
    return df["k"].max() + 1


def two2oneD(k, l, interaction_length):
    """Get index of state that starts at bp k and ends at bp l."""
    i = (
        k + 1
    ) * interaction_length  # all states (0:k,0:interaction_length)
    i -= (
        (k + 1) * k / 2
    )  # minus all states where l<k
    i -= (
        interaction_length - l - 1
    )  # minus states (k,l:interaction_length)
    return i


def treekin_rates_from_RRIkinDP_states(
    states_file,
    rate_file,
    energy_type="E",
    absorbing_states=[],
    absorbin_full_interaction=False,
    dissociation_at=None,
    absorbing_dissociated_state=False,
    energy_penalty_absorbing_state=19,
    binary=True,
    state_names_file=None,
):
    """Create rates input file for treekin from direct paths states output.

    states_file: states file from RRIkinDP
    rate_file: output file formated according to treekin rates input format
    energy_type: which energy colum in the states file to use
    absorbing_states: list of states that an absorbing states sould
        be attached to; states identified by there index in the states file;
        eg [12,1,24]
    absorbing_full_interaction: boolean; attach an absorbing state to the full
        interaction without needing to know the index of the full interaction
    dissociation_at: None if no dissociated state should be included; integer to
        specify from what interaction length an interaction can directly dissociate
        (eg set to 2 to skip single base pair state as barrier towards dissociation)
    absorbing_dissociated_state: boolean; attach an absorbing state to the
        dissociated state
    energy_penalty_absorbing_state: energy difference in kcal/mol between the
        absorbing state and the connected state; choose such that (large enough
        it results in an out rate that is negelectible within the simulated time
        while being nummerical stable; if the rate file is not written in a
        binary any energy difference ~larger than 10kcal/mol leads to a rate
        that is rounded to zero in the output format and thus not usable. For
        long simulation times look into high precision support in treekin
        (--mlapack-method).
    binary: bollean; whether to write output rate file in binary format
    """

    # Warning on absorbing state rates
    if (not binary) and (
        energy_penalty_absorbing_state > 10
    ):
        warnings.warn(
            "Reset energy penalty for absorbing states to 10kcal/mol. "
            + str(energy_penalty_absorbing_state)
            + "kcal/mol would lead to smaller rates than"
            + " what can be written to non-binary rate file. "
            + "With 10kcal/mol out-rates of absorbing states are not "
            + "negligible after 10E3 treekin time units."
        )
        energy_penalty_absorbing_state = 10

    # set up data structures with states info
    states = []
    states_dict = {}
    interaction_length = 0
    with open(states_file, "r") as f:
        i = 0
        for line in f:
            if i == 0:
                lables = line.strip().split("\t")
                index_e = lables.index(
                    energy_type
                )
                index_i = lables.index("k")
                index_j = lables.index("l")
            else:
                data = line.strip().split("\t")
                states.append(
                    [
                        float(data[index_e])
                        / 100,
                        (
                            int(data[index_i]),
                            int(data[index_j]),
                        ),
                    ]
                )
                states_dict[
                    (
                        int(data[index_i]),
                        int(data[index_j]),
                    )
                ] = {
                    "energy": float(data[index_e])
                    / 100,
                    "index": i - 1,
                }
                if (
                    int(data[index_j])
                    > interaction_length
                ):
                    interaction_length = int(
                        data[index_j]
                    )
            i = i + 1
    interaction_length += 1
    energies = [state[0] for state in states]

    # set up absorbing state
    if absorbin_full_interaction:
        absorbing_states.append(
            states_dict[
                0, interaction_length - 1
            ]["index"]
        )
    absorbing_states = list(set(absorbing_states))

    # count states
    number_of_states = len(states) + len(
        absorbing_states
    )
    if dissociation_at is not None:
        number_of_states += 1
    if absorbing_dissociated_state is not None:
        number_of_states += 1

    # build rate matrix
    matrix = []
    for k in range(len(states)):
        row = [0.0] * (len(states))

        current_i = states[k][1][0]
        current_j = states[k][1][1]

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
            row[l] = check_rate(
                get_rate(
                    energies[k], energies[l]
                ),
                (k, l),
                min_rate=MIN_RATE,
            )

        if absorbing_dissociated_state:
            if current_j - current_i == 0:
                row.append(
                    check_rate(
                        get_rate(
                            energies[k],
                            DISSOCIATED_STATE_ENERGY,
                        ),
                        (k, "dissociated-state"),
                        min_rate=MIN_RATE,
                    )
                )
            else:
                row.append(0.0)

        if (
            (dissociation_at is not None)
            and current_j - current_i
            < dissociation_at
        ):
            row.append(
                check_rate(
                    get_rate(
                        energies[k],
                        energies[k]
                        - energy_penalty_absorbing_state,
                    ),
                    (
                        k,
                        "dissociated-absorbing-state",
                    ),
                    min_rate=MIN_RATE,
                )
            )
        elif dissociation_at is not None:
            row.append(0.0)

        if k in absorbing_states:
            rates_to_absorbing = [0.0] * len(
                absorbing_states
            )
            index_absorbing = (
                absorbing_states.index(k)
            )
            rates_to_absorbing[
                index_absorbing
            ] = check_rate(
                get_rate(
                    energies[k],
                    energies[k]
                    - energy_penalty_absorbing_state,
                ),
                (k, "absorbing-state"),
                min_rate=MIN_RATE,
            )
            row += rates_to_absorbing
        else:
            row += [0.0] * len(absorbing_states)

        matrix.append(row)

    if absorbing_dissociated_state:
        row = [0.0] * (len(states) + 1)
        if dissociation_at is not None:
            row.append(
                check_rate(
                    get_rate(
                        DISSOCIATED_STATE_ENERGY,
                        DISSOCIATED_STATE_ENERGY
                        - energy_penalty_absorbing_state,
                    ),
                    (
                        "dissociated-state",
                        "dissociated-absorbing-state",
                    ),
                    min_rate=MIN_RATE,
                )
            )
        row += [0.0] * len(absorbing_states)
        for i in range(interaction_length):
            k = states_dict[(i, i)]["index"]
            row[k] = check_rate(
                get_rate(
                    DISSOCIATED_STATE_ENERGY,
                    states[k][0],
                ),
                ("dissociated-state", k),
                min_rate=MIN_RATE,
            )
        matrix.append(row)

    # add rate entries for a dissociated state
    if dissociation_at is not None:
        states.append(
            [DISSOCIATED_STATE_ENERGY, ("d", "d")]
        )
        row = [0.0] * number_of_states
        for i in range(interaction_length):
            for l in range(dissociation_at):
                if i + l < interaction_length:
                    k = states_dict[(i, i + l)][
                        "index"
                    ]
                    row[k] = check_rate(
                        get_rate(
                            states[k][0]
                            - energy_penalty_absorbing_state,
                            states[k][0],
                        ),
                        (
                            "dissociated-absorbing-state",
                            k,
                        ),
                        min_rate=MIN_RATE,
                    )
        if absorbing_dissociated_state:
            states.append(
                [
                    DISSOCIATED_STATE_ENERGY
                    - energy_penalty_absorbing_state,
                    ("a", "d"),
                ]
            )
            row[len(states)-1] = check_rate(
                get_rate(
                    DISSOCIATED_STATE_ENERGY
                    - energy_penalty_absorbing_state,
                    DISSOCIATED_STATE_ENERGY,
                ),
                (
                    "dissociated-absorbing-state",
                    "disscociated-state",
                ),
                min_rate=MIN_RATE,
            )
        matrix.append(row)

    # add rate entries for absorbing states (except for absorbing dissociated state)
    for a in absorbing_states:
        states.append(
            [
                states[a][0]
                - energy_penalty_absorbing_state,
                (
                    "a",
                    f"{states[a][1][0]}:{states[a][1][1]}",
                ),
            ]
        )
        row = [0.0] * number_of_states
        row[a] = check_rate(
            get_rate(
                states[a][0]
                - energy_penalty_absorbing_state,
                states[a][0],
            ),
            ("absorbing_state", a),
            min_rate=MIN_RATE,
        )
        matrix.append(row)

    # write matrix to file
    if binary:
        # transpose matrix
        t_matrix = [
            [
                matrix[j][i]
                for j in range(len(matrix))
            ]
            for i in range(len(matrix[0]))
        ]
        out = open(rate_file, "w+b")
        out.write(
            struct.pack("<i", number_of_states)
        )
        for row in t_matrix:
            for e in row:
                out.write(struct.pack("<d", e))
        out.close()
    else:
        out = open(rate_file, "w")
        for row in matrix:
            for e in row:
                out.write(format_rates(e) + " ")
            out.write("\n")
        out.close()

    # save state names
    if state_names_file is None:
        state_names_file = rate_file + ".json"
    state_names = [
        f"{state[1][0]}:{state[1][1]}"
        for state in states
    ]
    with open(state_names_file, "w") as json_file:
        json.dump(state_names, json_file)

    return matrix, state_names


# call treekin
def run_treekin(
    rate_file,
    start_state,
    binary=True,
    treekin_executable="treekin",
    write_treekin_output_files=True,
    treekin_output_file=None,
    verbose = False,
):
    """
    Run the treekin executable to compute RNA-RNA interaction dynamics.

    Parameters
    ----------
    rate_file : str
        Path to the rates input file in the specified format (binary or plain text).
    start_state : int
        Index of the starting state (0-based) in the interaction network.
    binary : bool, optional
        If True, the rates input file is in binary format, which offers higher precision (default is True).
    treekin_executable : str, optional
        Path to the treekin executable. Default assumes 'treekin' is available in the system PATH.
    write_treekin_output_files : bool, optional
        If True, writes the treekin output to a file specified by 'treekin_output_file' (default is True).
    treekin_output_file : str or None, optional
        Path to save the raw output from treekin. Required if 'write_treekin_output_files' is True.
    verbose : bool, optional
        If True, prints detailed output and errors from treekin execution (default is False).

    Returns
    -------
    None
        Executes treekin and saves its output if specified.

    Notes
    -----
    Ensure that the treekin executable is compatible with the rate file format and the command-line flags in use.
    """

    treekin_args = [
        treekin_executable,
        "-m",
        "I",
        "--ratesfile",
        rate_file,
        "--t8",
        "1E10",
        "--p0",
        str(start_state + 1) + "=1.0",
        "--mlapack-method",
        "DD"
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
    treekin_output = StringIO(
        stdout_data.decode()
    )

    if write_treekin_output_files:
        with open(
            treekin_output_file, "w"
        ) as out_handle:
            out_handle.write(
                treekin_output.read()
            )

    if verbose:
        print(stderr_data)


# TODO: get features from treekin output


def plot_treekin(
    treekin_output,
    treekin_plot,
    state_names=None,
    labels=True,
    label_cutoff_fraction=0.1,
    figsize=(7, 4),
    x_lim=(None, None),
    y_lim=(-0.05, 1.05),
):
    """Plot treekin output to file.

    Parameters
    ----------
    treekin_output : string
                     Path to treekin output file that is plotted.
    treekin_plot : string
                   Path that plot is saved to.
                   Fileformat is derived from ending. example: '.png', '.pdf'
    states_names : List of names for each state in the treekin output file.
    labels : bool
             Wether to plot state lables.
    label_cutoff_fraction : float between 0 and 1
                            Popluation fraction that a state has to represent at
                            at least one time point to get labled.
    figsize : duple of floats
              Figure size in inches.
    x_lim : duple of floats
            Plot range on x-axis.
    y_lim = duple of floats
            Plot range on y-axis.
    Todo: states labels as text instead of marker
    Notes: lables could also be done with https://pypi.org/project/matplotlib-label-lines/
    """
    # read treekin output file
    if state_names is not None:
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

    df = df.iloc[
        :, :-1
    ]  # remove empty column (tailing spaces in input)

    # set figure size
    f, ax = plt.subplots(figsize=figsize)

    # set axis labels
    ax.set_ylabel("Population")
    ax.set_xlabel("Time (a.u.)")

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

        edge_col = "white"
        text = col

        # mark absorbing states with grey marker edge and remove leading "a:" in state name
        if col.startswith("a"):
            edge_col = "grey"
            text = ":".join(col.split(":")[1:])

        # set text size in state labels
        text_markersize = 12
        if len(text)<4:
            text_markersize = 10
            if len(text)<2:
                text_markersize = 7

        # plot state label background
        plt.plot(
            max_population_time,
            max_population,
            marker="o",
            color="white",
            markeredgecolor=edge_col,
            alpha=0.7,
            markersize=14,
        )

        # plot state label text
        plt.plot(
            max_population_time,
            max_population,
            marker="$%s$" % text,
            color=color,
            markersize=text_markersize,
            # markeredgecolor="black",
            # markeredgewidth=0.05,
        )

    # set time axis to logarithmic scale
    ax.set_xscale("log")

    # set plot range
    ax.set_xlim(x_lim[0], x_lim[1])
    ax.set_ylim(y_lim[0], y_lim[1])

    # save figure
    f.savefig(treekin_plot, bbox_inches="tight")
    plt.close(f)


# TODO: run treekin for all seeds and summarize features


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Compute RNA-RNA interaction formation dynamics by solving the master equation."
    )

    # Required arguments
    parser.add_argument(
        "-s", "--states",
        help="Path to input states file generated by RRIkinDP.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-i", "--initial",
        help="Index of the initial state in the interaction network.",
        type=int,
        required=True,
    )
    parser.add_argument(
        "-r", "--rates",
        help="Path to save generated rates file , formatted for treekin.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-p", "--probs",
        help="Path to save treekin's state probabilities output.",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-f", "--figure",
        help="Path to save the plot of state probabilities over time.",
        type=str,
        required=True,
    )

    # Rate calculation settings
    parser.add_argument(
        "--energy_type",
        help="Specifies the energy column to use from the states file (e.g., 'E', 'ED1').",
        type=str,
        default="E",
    )
    '''
    parser.add_argument(
        "--min_rate",
        help="Minimum allowed rate value; rates below this threshold will raise a warning.",
        type=float,
        default=1e-14,
    )
    '''

    # Absorbing state settings
    parser.add_argument(
        "--absorbing_states",
        help="List of indices representing states to which absorbing states should be attached. Example: '12,1,24'.",
        type=lambda x: list(map(int, x.split(","))),
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
        default=None,
    )
    parser.add_argument(
        "--absorbing_dissociated_state",
        help="Attach an absorbing state to the dissociated state. (True/False).",
        action="store_true",
    )
    parser.add_argument(
        "--energy_penalty_absorbing_state",
        help="Energy penalty in kcal/mol for the absorbing state relative to connected state. Recommended values are <= 10 kcal/mol for non-binary outputs.",
        type=float,
        default=19,
    )

    # Treekin solver options
    parser.add_argument(
        "--treekin_executable",
        help="Path to the treekin executable (default: 'treekin').",
        type=str,
        default="treekin",
    )
    parser.add_argument(
        "--binary_rate_file",
        help="Output rate file in binary format for higher precision. (True/False).",
        action="store_true",
    )

    parser.add_argument(
        "--treekin_verbose",
        help="Print additional details from treekin's execution. (True/False).",
        action="store_true",
    )

    # Plotting options
    parser.add_argument(
        "--state_names_file",
        help="Path to save JSON file with state names. Defaults to <rate_file>.json.",
        type=str,
        default=None,
    )

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
        help="Plot range on x-axis as a tuple. (None for automatic limits).",
        type=lambda x: tuple(map(lambda v: None if v == None else float(v), x.split())),
        default=(None, None),
    )
    parser.add_argument(
        "--plot_y_lim",
        help="Plot range on y-axis as a tuple (default: -0.05, 1.05).",
        type=lambda x: tuple(map(float, x.split())),
        default=(-0.05, 1.05),
    )

    args = parser.parse_args()

    # Generate rate matrix
    matrix, state_names = treekin_rates_from_RRIkinDP_states(
        states_file=args.states,
        rate_file=args.rates,
        energy_type=args.energy_type,
        absorbing_states=args.absorbing_states,
        absorbin_full_interaction=args.absorbing_full_interaction,
        dissociation_at=args.dissociation_at,
        absorbing_dissociated_state=args.absorbing_dissociated_state,
        energy_penalty_absorbing_state=args.energy_penalty_absorbing_state,
        binary=args.binary_rate_file,
        state_names_file=args.state_names_file,
    )

    # Run treekin
    run_treekin(
        rate_file=args.rates,
        start_state=args.initial,
        binary=args.binary_rate_file,
        treekin_executable=args.treekin_executable,
        write_treekin_output_files=True,
        treekin_output_file=args.probs,
        verbose=args.treekin_verbose,
    )

    # Plot state probabilities
    plot_treekin(
        treekin_output=args.probs,
        treekin_plot=args.figure,
        state_names=state_names,
        labels=args.plot_no_labels,
        label_cutoff_fraction=args.plot_label_cutoff,
        figsize=args.figsize,
        x_lim=args.plot_x_lim,
        y_lim=args.plot_y_lim,
    )
