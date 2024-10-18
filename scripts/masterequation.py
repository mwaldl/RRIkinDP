import argparse
import warnings
import subprocess
from io import StringIO
import math
import pandas as pd
import matplotlib.pyplot as plt
#import seaborn as sn
import json

R = 1.98720425864083 * math.pow(10, -3)  # gas contant in dagcal⋅K−1⋅mol−1
T = 273.15  # 0 Celsius in K
MIN_RATE = 10 ** (-14) # TODO: unit
DISSOCIATED_STATE_ENERGY = 0 # in dagcal⋅K−1⋅mol−1


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
    i = (k + 1) * interaction_length  # all states (0:k,0:interaction_length)
    i -= (k + 1) * k / 2  # minus all states where l<k
    i -= interaction_length - l - 1  # minus states (k,l:interaction_length)
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
    if (not binary) and (energy_penalty_absorbing_state > 10):
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
                index_e = lables.index(energy_type)
                index_i = lables.index("k")
                index_j = lables.index("l")
            else:
                data = line.strip().split("\t")
                states.append(
                    [
                        float(data[index_e]) / 100,
                        (int(data[index_i]), int(data[index_j])),
                    ]
                )
                states_dict[(int(data[index_i]), int(data[index_j]))] = {
                    "energy": float(data[index_e]) / 100,
                    "index": i - 1,
                }
                if int(data[index_j]) > interaction_length:
                    interaction_length = int(data[index_j])
            i = i + 1
    interaction_length += 1
    energies = [state[0] for state in states]

    # set up absorbing state
    if absorbin_full_interaction:
        absorbing_states.append(states_dict[0, interaction_length - 1]["index"])
    absorbing_states = list(set(absorbing_states))

    # count states
    number_of_states = len(states) + len(absorbing_states)
    if dissociation_at is not "None":
        number_of_states += 1
    if absorbing_dissociated_state is not "None":
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
                get_rate(energies[k], energies[l]), (k, l), min_rate=MIN_RATE
            )

        if absorbing_dissociated_state:
            if current_j - current_i == 0:
                row.append(
                    check_rate(
                        get_rate(energies[k], DISSOCIATED_STATE_ENERGY),
                        (k, "dissociated-state"),
                        min_rate=MIN_RATE,
                    )
                )
            else:
                row.append(0.0)

        if (dissociation_at is not "None") and current_j - current_i < dissociation_at:
            row.append(
                check_rate(
                    get_rate(
                        energies[k],
                        energies[k] - energy_penalty_absorbing_state,
                    ),
                    (k, "dissociated-absorbing-state"),
                    min_rate=MIN_RATE,
                )
            )
        elif dissociation_at is not "None":
            row.append(0.0)

        if k in absorbing_states:
            rates_to_absorbing = [0.0] * len(absorbing_states)
            index_absorbing = absorbing_states.index(k)
            rates_to_absorbing[index_absorbing] = check_rate(
                get_rate(energies[k], energies[k] - energy_penalty_absorbing_state),
                (k, "absorbing-state"),
                min_rate=MIN_RATE,
            )
            row += rates_to_absorbing
        else:
            row += [0.0] * len(absorbing_states)

        matrix.append(row)

    if absorbing_dissociated_state:
        row = [0.0] * (len(states) + 1)
        if dissociation_at is not "None":
            row.append(
                check_rate(
                    get_rate(
                        DISSOCIATED_STATE_ENERGY,
                        DISSOCIATED_STATE_ENERGY - energy_penalty_absorbing_state,
                    ),
                    ("dissociated-state", "dissociated-absorbing-state"),
                    min_rate=MIN_RATE,
                )
            )
        row += [0.0] * len(absorbing_states)
        for i in range(interaction_length):
            k = states_dict[(i, i)]["index"]
            row[k] = check_rate(
                get_rate(DISSOCIATED_STATE_ENERGY, states[k][0]),
                ("dissociated-state", k),
                min_rate=MIN_RATE,
            )
        matrix.append(row)

    # add rate entries for a dissociated state
    if dissociation_at is not "None":
        states.append([DISSOCIATED_STATE_ENERGY, ("d", "d")])
        row = [0.0] * number_of_states
        for i in range(interaction_length):
            for l in range(dissociation_at):
                if i + l < interaction_length:
                    k = states_dict[(i, i + l)]["index"]
                    row[k] = check_rate(
                        get_rate(
                            states[k][0] - energy_penalty_absorbing_state,
                            states[k][0],
                        ),
                        ("dissociated-absorbing-state", k),
                        min_rate=MIN_RATE,
                    )
        if absorbing_dissociated_state:
            states.append([DISSOCIATED_STATE_ENERGY- energy_penalty_absorbing_state, ("a", "d")])
            row[len(states)] = check_rate(
                get_rate(
                    DISSOCIATED_STATE_ENERGY - energy_penalty_absorbing_state,
                    DISSOCIATED_STATE_ENERGY,
                ),
                ("dissociated-absorbing-state", "disscociated-state"),
                min_rate=MIN_RATE,
            )
        matrix.append(row)

    # add rate entries for absorbing states (except for absorbing dissociated state)
    for a in absorbing_states:
        states.append([states[a][0] - energy_penalty_absorbing_state, ("a", f"{states[a][1][0]}:{states[a][1][1]}")])
        row = [0.0] * number_of_states
        row[a] = check_rate(
            get_rate(states[a][0] - energy_penalty_absorbing_state, states[a][0]),
            ("absorbing_state", a),
            min_rate=MIN_RATE,
        )
        matrix.append(row)

    # write matrix to file
    if binary:
        # transpose matrix
        t_matrix = [
            [matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))
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
                out.write(format_rates(e) + " ")
            out.write("\n")
        out.close()

    # save state names
    if state_names_file == None:
        state_names_file = rate_file + ".json"
    state_names = [f"{state[1][0]}:{state[1][1]}" for state in states]
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
):

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
    stdout_data, stderr_data = treekin_process.communicate()
    stderr_data = stderr_data.decode()
    treekin_output = StringIO(stdout_data.decode())

    if write_treekin_output_files:
        with open(treekin_output_file, "w") as out_handle:
            out_handle.write(treekin_output.read())

    print(stderr_data)


# TODO: get features from treekin output

def plot_treekin(treekin_output, treekin_plot):
    """Plot treekin output to file.

    Parameters
    ----------
    treekin_output : string
                     Path to treekin output file that is plotted.
    treekin_plot : string
                   Path that plot is saved to.
                   Fileformat is derived from ending. example: '.png', '.pdf'

    Todo: states lables, figure size
    """
    data = pd.read_csv(treekin_output, header=None, sep=" ", comment="#")
    data = data.iloc[:, :-1]  # remove empty column (tailing spaces in input)
    f = plt.figure()
    # f.set_size_inches(7, 5)
    data[len(data.columns)] = data[len(data.columns) - 1] + data[len(data.columns) - 2]
    plt.plot(data[0], data.loc[:, data.columns != 0])
    plt.ylabel("Population")
    plt.xlabel("Time (a.u.)")
    # plt.legend(
    #     title="states:",
    #     labels=data.columns.values.tolist()[1:],
    #     loc="upper right",
    # )
    plt.ylim(-0.1, 1.1)
    plt.xscale("log")
    f.savefig(treekin_plot, bbox_inches="tight")
    plt.close(f)


# TODO: run treekin for all seeds and summarize features


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Compute RNA-RNA interaction formation dynamics by solving the master equation.")
    parser.add_argument("states", help="filepath to input states file as generated by RRIkinDP", type=str)
    parser.add_argument("start", help="index of initial state", type=int)
    parser.add_argument("rates", help="filepath to save rate file to", type=str)
    parser.add_argument("out", help="filepath to save state probabilities to", type=str)
    parser.add_argument("figure", help="filepath to save state probabilities plot to", type=str)

    args = parser.parse_args()
    states_file = args.states
    rates_file = args.rates
    treekin_out = args.out
    treekin_plot = args.figure
    start_state = args.start

    treekin_rates_from_RRIkinDP_states(
        states_file,
        rates_file,
        energy_type="E",
        absorbing_states=[],
        absorbin_full_interaction=False,
        dissociation_at=2,
        absorbing_dissociated_state=False,
        energy_penalty_absorbing_state=19,
        binary=False,
    )
    run_treekin(
        rates_file,
        start_state,
        binary=False,
        write_treekin_output_files=True,
        treekin_output_file=treekin_out,
    )
    plot_treekin(treekin_out, treekin_plot)
