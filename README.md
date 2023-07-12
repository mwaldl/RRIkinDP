# RRIkinDP

### example call

```
./paths --id_b "ChiX" --seq_b "acaccgucgcuuaaagugacggcauaauaauaaaaaaaugaaauuccucuuugacgggccaauagcgauauuggccauuuuuuu" --id_a "b1737" --seq_a "GUUUGUUACCCAACAAACCGGUUGAAGUAAUUGACUCGCUGCUUUAUGGCAAAGUCGAUGGUUUAGGCGUGCUUAAGGCUGCGGUUGCAGCGAUUAAAAAAGCCGCAGCAAAUUAAUUUAUUUUAAAUUUUCCCGUCAAAGAGUUAUUUCAUAAAUCAAUACCGCAAUAUUUAAAUUGCGGUUUUUAAGGGUAUUUUUCUAUGAGUAAUGUUAUUGCAUCGCUUGAAAAGGUACUCCUCCCUUUUGCAGUUAAAAUAGGAAAGCAGCCACACGUUAAUGCAAUCAAAAAUGGCUUUAUUC" --interaction_bps "(134,56):(135,55):(136,54):(137,53):(138,52):(139,51):(140,50):(141,49):(142,48):(143,47):(146,44):(147,43):(148,42):(149,41):(150,40):(151,39):(152,38)" --seed 4  --write_states "states.tsv"  --write_structures test.fa --write_all_barriers "barriers.tsv" --fixed_intramolecular_structures --str_b "......(((((...)))))............(((((((((.................(((((((....))))))))))))))))"
```

## Recommended Usage: Docker

The following steps create local Docker image called `rrikindp` containing all dependencies, which allows the program to be called using the command line.

### 1. Docker image creation

1. Clone this repository
2. Navigate to the repository directory
3. Call `docker build -t rrikindp .`

The Docker container is then available on your computer globally.

### 2. Call program using Docker image in data folder

Docker images cannot access the full file system of the computer they are running on. If file outputs of `paths` are of interest, they need to be written to a location accessible to both the host and the Docker container. In this example, a folder for data results is called `dataresults` on the host filesystem and `/data` on the container.

1. Navigate to a folder where you'd like your result files to be located
2. Create a subfolder `dataresults`: `mkdir dataresults`
3. Call the program inside the Docker container, mapping the `dataresults` folder to the container filesystem using `-v`

   ```
   docker run \
     -v "`pwd`/dataresults":/data \
     -t rrikindp2 /paths  \
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

`Dockerfile` contains the relevant commands for building `src/directPaths.cpp` in a Debian-like environment using `miniconda3` as a package manager for installing most dependencies listed below.

### Dependencies


#### RRIkinDP
- compiler supporting C++11 standard and OpenMPg++ (fopenmp)
- boost C++ library version >= 1.50.0 (lboost)
- IntaRNA (lIntaRNA)
  - https://github.com/BackofenLab/IntaRNA/#install
- Vienna RNA package version >= 2.4.14 (lRNA)
  - https://www.tbi.univie.ac.at/RNA/#download
  - https://github.com/ViennaRNA/ViennaRNA
  - conda install -c bioconda viennarna
- Easylogging++ logging framework (leasylogging)
  - https://github.com/amrayn/easyloggingpp

#### landscapes.py
- python3
- matplotlib
- seaborn
- pandas


### Parameters

--help                             Display this help message  
--version                          Display the version number  
--id_a arg                         id of first sequence  
--seq_a arg                        first sequence  
--str_a arg                        intramolecular structure of first sequence
                                   in dotbraket notation  
--id_b arg                         id of second sequenc  
--seq_b arg                        second sequence  
--str_b arg                        intramolecular structure of second
                                   sequence in dotbraket notation  
--interaction_bps arg              interaction base pair list as string (one
                                   based)  
--seed arg                         seed length  
--write_all_barriers arg           file paths to write minmal barriers for
                                   all seeds to  
--write_states arg                 file paths to write states and their
                                   energies to  
--compute_states_only              compute states and output their energies
                                   to path specified in 'write_states'  
--no_dangle                        turn off dangle contributions at
                                   interaction ends  
--fixed_intramolecular_structures  compute accessibilities based on fixed
                                   intramolecular structures instead of based
                                   on partition function  
--write_structures arg             file paths to write intramolecular and
                                   fully extended intermolecular structures
                                   to  
--temperature arg (=37)            temperature in Celsius  
