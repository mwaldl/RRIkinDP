# RRIkinDP

RRIkinDP models RNA-RNA interaction formation between two RNAs along direct paths from a
first base pair (or interaction seed) to a full input (candidate) interaction. It evaluates
thermodynamic and kinetic features derived from the underlying state space.

## Contents

- [What is included](#what-is-included)
- [Installation with conda](#installation-with-conda)
- [Installation and Usage with Docker](#installation-and-usage-with-docker)
- [Building yourself](#building-yourself)
- [Run RRIkinDP](#run-rrikindp)
- [Python interface](#python-interface)
- [Plot energy landscape](#plot-energy-landscape)
- [Citation](#citation)
- [License](#license)

## What is included

RRIkinDP consists of two components, which can be installed and used independently of each
other:

- **`RRIkinDP`, the command line tool**, computes the energies of all states of the state
  space and the minimal energy barrier of the best direct path from a given start
  interaction to the full input interaction. States, structures and barriers are written to
  tab separated files. It is written in C++ and does not require python; continue at
  [Run RRIkinDP](#run-rrikindp).
- **`rrikindp`, the python module**, gives access to the same state space and barrier
  computation from within python, provides utilities for generating its input with IntaRNA,
  plots the energy landscape of the state space and models the interaction formation as a
  markov process to obtain kinetic features; continue at
  [Python interface](#python-interface).

## Installation with conda

Provides the C++ tool and the python interface:

`conda install bioconda::rrikindp` 

## Installation and Usage with Docker

The following steps create local Docker image called `rrikindp` containing all dependencies, which allows the program to be called using the command line. The image provides the command line tool only; for the python module see [Installation with conda](#installation-with-conda) or [Building yourself](#building-yourself).

### 1. Docker image creation

1. Clone this repository
2. Navigate to the repository directory
3. Call `docker build -t rrikindp .`

The Docker container is then available on your computer globally.

### 2. Call program using Docker image in data folder

Docker images cannot access the full file system of the computer they are running on. If the file outputs of `RRIkinDP` are of interest, they need to be written to a location accessible to both the host and the Docker container. In this example, a folder for data results is called `dataresults` on the host filesystem and `/data` on the container.

1. Navigate to a folder where you'd like your result files to be located
2. Create a subfolder `dataresults`: `mkdir dataresults`
3. Call the program inside the Docker container, mapping the `dataresults` folder to the container filesystem using `-v`

   ```
   docker run \
     -v "`pwd`/dataresults":/data \
     -t rrikindp /RRIkinDP  \
     --id_b "ChiX" \
     --seq_b "acaccgucgcuuaaagugacggcauaauaauaaaaaaaugaaauuccucuuugacgggccaauagcgauauuggccauuuuuuu" \
     --id_a "b1737" \
     --seq_a "GUUUGUUACCCAACAAACCGGUUGAAGUAAUUGACUCGCUGCUUUAUGGCAAAGUCGAUGGUUUAGGCGUGCUUAAGGCUGCGGUUGCAGCGAUUAAAAAAGCCGCAGCAAAUUAAUUUAUUUUAAAUUUUCCCGUCAAAGAGUUAUUUCAUAAAUCAAUACCGCAAUAUUUAAAUUGCGGUUUUUAAGGGUAUUUUUCUAUGAGUAAUGUUAUUGCAUCGCUUGAAAAGGUACUCCUCCCUUUUGCAGUUAAAAUAGGAAAGCAGCCACACGUUAAUGCAAUCAAAAAUGGCUUUAUUC" \
     --interaction_bps "(134,56):(135,55):(136,54):(137,53):(138,52):(139,51):(140,50):(141,49):(142,48):(143,47):(146,44):(147,43):(148,42):(149,41):(150,40):(151,39):(152,38)" \
     --seed 4 \
     --fixed_intramolecular_structures \
     --str_b "......(((((...)))))............(((((((((.................(((((((....))))))))))))))))" \
     --write_states /data/states.tsv \
     --write_structures /data/test.fa \
     --write_all_barriers /data/barriers.tsv
     ```

4. You will see terminal output and result files in `dataresults`


## Building yourself

`Dockerfile` contains the relevant commands for building `src/rrikindp/RRIkinDP.cpp` in a Debian-like environment using `conda` as a package manager for installing most dependencies listed below.

The command line tool and the python module are built independently of each other.

### Build the command line tool

The `RRIkinDP` executable only requires the dependencies listed under `RRIkinDP command line tool`; neither python nor the dependencies of the python module are needed.

```
cd src/rrikindp
make
```

If the dependencies are not installed system wide, their include and library paths are passed via the usual variables:

```
make CXXFLAGS="-I/your/prefix/include" LDFLAGS="-L/your/prefix/lib"
```

The resulting `RRIkinDP` binary can then be copied to any directory in your `PATH`.

### Install the python module

The `rrikindp` python module additionally requires the dependencies listed under `rrikindp python module` and is installed from the repository root:

```
pip install .
```

This builds and installs the python module only; the command line tool is built with `make` as described above.

### Dependencies

The dependencies can be installed by any means, e.g. from source, via the system package manager or via conda.

#### RRIkinDP command line tool
- compiler supporting the C++17 standard and OpenMP, e.g. g++ (-fopenmp)
- boost C++ library version >= 1.50.0 (lboost)
- IntaRNA (lIntaRNA)
  - https://github.com/BackofenLab/IntaRNA/#install
- Vienna RNA package version >= 2.4.14 (lRNA)
  - https://www.tbi.univie.ac.at/RNA/#download
  - https://github.com/ViennaRNA/ViennaRNA
  - conda install -c bioconda viennarna
- Easylogging++ logging framework (leasylogging)
  - https://github.com/amrayn/easyloggingpp

#### rrikindp python module
- python3 >= 3.9
- pybind11 (build time only)
- matplotlib
- seaborn
- pandas
- treekin version >= 0.5.1 (only required for `MarkovProcess`)
  - https://github.com/ViennaRNA/Treekin
  - conda install -c bioconda treekin



## Run RRIkinDP


### example call

```
RRIkinDP \
--id_b "ChiX" \
--seq_b "acaccgucgcuuaaagugacggcauaauaauaaaaaaaugaaauuccucuuugacgggccaauagcgauauuggccauuuuuuu" \
--id_a "b1737" \
--seq_a "GUUUGUUACCCAACAAACCGGUUGAAGUAAUUGACUCGCUGCUUUAUGGCAAAGUCGAUGGUUUAGGCGUGCUUAAGGCUGCGGUUGCAGCGAUUAAAAAAGCCGCAGCAAAUUAAUUUAUUUUAAAUUUUCCCGUCAAAGAGUUAUUUCAUAAAUCAAUACCGCAAUAUUUAAAUUGCGGUUUUUAAGGGUAUUUUUCUAUGAGUAAUGUUAUUGCAUCGCUUGAAAAGGUACUCCUCCCUUUUGCAGUUAAAAUAGGAAAGCAGCCACACGUUAAUGCAAUCAAAAAUGGCUUUAUUC" \
--interaction_bps "(134,56):(135,55):(136,54):(137,53):(138,52):(139,51):(140,50):(141,49):(142,48):(143,47):(146,44):(147,43):(148,42):(149,41):(150,40):(151,39):(152,38)" \
--seed 4 \
--write_states "states.tsv" \
--write_structures "test.fa" \
--write_all_barriers "barriers.tsv" \
--fixed_intramolecular_structures \
--str_b "......(((((...)))))............(((((((((.................(((((((....))))))))))))))))"
```

### Parameters
```
RRIkinDP Usage:
  --help                             Display this help message
  --version                          Display the version number
  --id_a arg                         id of first sequence
  --seq_a arg                        first sequence
  --str_a arg                        intramolecular structure of first sequence
                                     in dotbraket notation
  --id_b arg                         id of second sequence
  --seq_b arg                        second sequence
  --str_b arg                        intramolecular structure of second
                                     sequence in dotbraket notation
  --interaction_bps arg              interaction base pair list as string (one
                                     based)
  --seed arg                         seed length
  --write_all_barriers arg           file path to write minimal barriers for
                                     all seeds to
  --write_states arg                 file path to write states and their
                                     energies to
  --compute_states_only              compute states and output their energies
                                     to path specified in 'write_states'
  --no_dangle                        turn off dangle contributions at
                                     interaction ends
  --fixed_intramolecular_structures  compute accessibilities based on fixed
                                     intramolecular structures instead of based
                                     on partition function
  --write_structures arg             file path to write intramolecular and
                                     fully extended intermolecular structures
                                     to
  --temperature arg (=37)            temperature in Celsius
```

## Python interface

The `rrikindp` module provides the state space and barrier computation of RRIkinDP as
a python class, utilities to prepare its input from IntaRNA, plotting of the energy
landscape and modeling of the interaction formation as markov process.

```
import rrikindp
```

### Energy landscape

`DPLandscape` gives access to the RRIkinDP state space without writing intermediate
files. It is set up from the two sequences and the interaction base pair list; base
pair indices are zero based. All energies are returned in kcal/mol.

```
mRNA = 'GUUUGUUACCCAACAAACCGGUUGAAGUAAUUGACUCGCUGCUUUAUGGCAAAGUCGAUGGUUUAGGCGUGCUUAAGGCUGCGGUUGCAGCGAUUAAAAAAGCCGCAGCAAAUUAAUUUAUUUUAAAUUUUCCCGUCAAAGAGUUAUUUCAUAAAUCAAUACCGCAAUAUUUAAAUUGCGGUUUUUAAGGGUAUUUUUCUAUGAGUAAUGUUAUUGCAUCGCUUGAAAAGGUACUCCUCCCUUUUGCAGUUAAAAUAGGAAAGCAGCCACACGUUAAUGCAAUCAAAAAUGGCUUUAUUC'
sRNA = 'acaccgucgcuuaaagugacggcauaauaauaaaaaaaugaaauuccucuuugacgggccaauagcgauauuggccauuuuuuu'
bp_list_intarna = '(134,56):(135,55):(136,54):(137,53):(138,52):(139,51):(140,50):(141,49):(142,48):(143,47):(146,44):(147,43):(148,42):(149,41):(150,40):(151,39):(152,38)'
seed_length = 4

bp_list = rrikindp.utilities.intarna_to_bplist(bp_list_intarna, zero_based=True)

el = rrikindp.DPLandscape(mRNA, sRNA, bp_list)

full_E = el.get_full_E()
seed_Es = el.get_seed_Es(seed_length)
min_barrier_energies = el.get_min_barrier_Es(seed_length)
seed, barrier_state = el.get_seed_barrier_state(seed_length)

el.save_states('states.tsv')
```

The optional arguments of `DPLandscape` correspond to the command line options:
`id_a`, `id_b`, `str_a`, `str_b`, `accessibility_from_pf`, `dangles` and `temperature`.

Further accessors are `get_full_hybridE()`, `get_full_ED()`, `get_full_ED1()`,
`get_full_ED2()`, `get_seed_EDs()`, `get_seed_ED1s()`, `get_seed_ED2s()` and
`get_states()`. The underlying C++ classes `EM`, `RnaSequence`, `Interaction` and
`BasePair` are exposed as well; they use decacalories per mol, as the files written by
the command line tool do.

### Utilities

`rrikindp.utilities` prepares the RRIkinDP input from IntaRNA predictions.

```
from rrikindp import utilities

# run IntaRNA; returns a pandas DataFrame for outMode 'C', the raw output otherwise
prediction = utilities.run_intarna(mRNA, sRNA, id1='b1737', id2='ChiX')
bp_list_intarna = prediction['bpList'][0]

# base pair list as string to list of tuples; zero based for DPLandscape
bp_list = utilities.intarna_to_bplist(bp_list_intarna, zero_based=True)

# three line representation of the interaction; expects one based base pairs
structure = utilities.get_RRI_string_representation(
    mRNA, sRNA, utilities.intarna_to_bplist(bp_list_intarna),
    id1='b1737', id2='ChiX'
)
```

`run_intarna()` takes sequences or paths to fasta files and forwards further options to
IntaRNA via `intarna_args`; `temperature`, `threads`, `outMode`, `outCsvCols`,
`out_file` and `intarna_executable` are available as keyword arguments.

### Interaction formation as markov process

`MarkovProcess` models the interaction formation on the state space of direct paths and
evaluates it with treekin. States are represented by `State` objects; absorbing states
and a dissociated state are added to the state space on request.

```
from rrikindp import MarkovProcess

# rate matrix for treekin; returns the matrix and the list of states
rates, states = MarkovProcess.generate_treekin_rates_file(
    'states.tsv', 'rates.bin', absorbing_full_interaction=True
)

# start the simulation in the dissociated state
start_state = [s for s in states if s.base_pairs == ('d', 'd')][0]
MarkovProcess.run_treekin(
    'rates.bin', [[start_state, 1.0]], treekin_output_file='probs.txt'
)

# kinetic features and plots
features = MarkovProcess.get_treekin_features(
    'probs.txt', states, eval_full_interaction=True, eval_dissociated_state=True
)
MarkovProcess.plot_treekin('probs.txt', 'probs.pdf', states=states)
MarkovProcess.plot_E_mean('probs.txt', states, figure_path='E_mean.pdf')
```

`generate_treekin_rates_file()` controls the state space via `absorbing_states`,
`absorbing_full_interaction`, `dissociation_at`, `absorbing_dissociated_state`,
`association_scaling_factor` and `energy_penalty_absorbing_state`; `run_treekin()`
exposes the treekin numerics via `precision`, `time_increment`, `sim_start_time`,
`sim_end_time` and `temperature`. Mean energies over time are available from
`get_E_mean()` and `get_E_mean_features()`.

## Plot energy landscape

The energylandscape of all interaction structures within the RRIkinDP state space can be plotted with the `landscape` module of the installed `rrikindp` package.


### minimal call

```
python -m rrikindp.landscape states.tsv
```

### with annotated interaction structure

Structure annotations are ploted when both sequences and the interaction base pairs are provided.

```
python -m rrikindp.landscape states.tsv \
--id1 b1737 \
--id2 ChiX \
--seq1 GUUUGUUACCCAACAAACCGGUUGAAGUAAUUGACUCGCUGCUUUAUGGCAAAGUCGAUGGUUUAGGCGUGCUUAAGGCUGCGGUUGCAGCGAUUAAAAAAGCCGCAGCAAAUUAAUUUAUUUUAAAUUUUCCCGUCAAAGAGUUAUUUCAUAAAUCAAUACCGCAAUAUUUAAAUUGCGGUUUUUAAGGGUAUUUUUCUAUGAGUAAUGUUAUUGCAUCGCUUGAAAAGGUACUCCUCCCUUUUGCAGUUAAAAUAGGAAAGCAGCCACACGUUAAUGCAAUCAAAAAUGGCUUUAUUC \
--seq2 acaccgucgcuuaaagugacggcauaauaauaaaaaaaugaaauuccucuuugacgggccaauagcgauauuggccauuuuuuu \
--bp_list '(134,56):(135,55):(136,54):(137,53):(138,52):(139,51):(140,50):(141,49):(142,48):(143,47):(146,44):(147,43):(148,42):(149,41):(150,40):(151,39):(152,38)' \
--out test.pdf
```


### additional parameters
```
usage: landscape.py [-h] [-o OUT] [--id1 ID1] [--id2 ID2] [--energy ENERGY] [--seq1 SEQ1] [--seq2 SEQ2] [--bp_list BP_LIST] [--e_min E_MIN]
                    [--e_max E_MAX] [--figuresize FIGURESIZE FIGURESIZE] [--annotate | --remove_annotations]
                    states

Plot energylandscape from states provided by RRIkinDP.

positional arguments:
  states                path to states file from RRIkinDP

optional arguments:
  -h, --help            show this help message and exit
  -o OUT, --out OUT     path to output figure
  --id1 ID1             name of first RNA
  --id2 ID2             name of second RNA
  --energy ENERGY       free energy to plot: E, ED1, ED2 or Ehybrid
  --seq1 SEQ1           sequence of first RNA
  --seq2 SEQ2           sequence of second RNA
  --bp_list BP_LIST     base pair list as string; example: (134,56):(135,55):(136,54):(137,53):(138,52):(139,51)
  --e_min E_MIN         path to output figure
  --e_max E_MAX         path to output figure
  --figuresize FIGURESIZE FIGURESIZE
                        figure dimensions in inches; provide two arguments seperated by space; width hight
  --annotate            annotate energy per interaction structure; default: do not annotate if more than 15 base pairs in interaction
  --remove_annotations  do not annotate energy per interaction structure; default: annotate if less than 16 base pairs in interaction
```

### example usage in a python script

```
import rrikindp


id1 = 'b1737'
id2 = 'ChiX'
seq1 = 'GUUUGUUACCCAACAAACCGGUUGAAGUAAUUGACUCGCUGCUUUAUGGCAAAGUCGAUGGUUUAGGCGUGCUUAAGGCUGCGGUUGCAGCGAUUAAAAAAGCCGCAGCAAAUUAAUUUAUUUUAAAUUUUCCCGUCAAAGAGUUAUUUCAUAAAUCAAUACCGCAAUAUUUAAAUUGCGGUUUUUAAGGGUAUUUUUCUAUGAGUAAUGUUAUUGCAUCGCUUGAAAAGGUACUCCUCCCUUUUGCAGUUAAAAUAGGAAAGCAGCCACACGUUAAUGCAAUCAAAAAUGGCUUUAUUC'
seq2 = 'acaccgucgcuuaaagugacggcauaauaauaaaaaaaugaaauuccucuuugacgggccaauagcgauauuggccauuuuuuu'
bp_list_intarna = '(134,56):(135,55):(136,54):(137,53):(138,52):(139,51):(140,50):(141,49):(142,48):(143,47):(146,44):(147,43):(148,42):(149,41):(150,40):(151,39):(152,38)'
states_file = 'states.tsv'

bp_list = rrikindp.utilities.intarna_to_bplist(bp_list_intarna)
structure = rrikindp.utilities.get_RRI_string_representation(
    seq1, seq2, bp_list, id1=id1, id2=id2
)

rrikindp.plot_landscape(states_file,
                        figure_path=None,
                        energy="E",
                        structure=structure,
                        seqId1=id1,
                        seqId2=id2,
                       )
```

## Citation

RRIkinDP is described in *RRIkinDP: Targeted RNA-RNA interaction kinetics*, currently
available as a preprint:

https://doi.org/10.1101/2023.07.28.548983

## License

RRIkinDP is distributed under the GNU General Public License v3.0; see `COPYING`.
