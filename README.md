# RRIkinDP

### example call
``
./paths --id_b "ChiX" --seq_b "acaccgucgcuuaaagugacggcauaauaauaaaaaaaugaaauuccucuuugacgggccaauagcgauauuggccauuuuuuu" --id_a "b1737" --seq_a "GUUUGUUACCCAACAAACCGGUUGAAGUAAUUGACUCGCUGCUUUAUGGCAAAGUCGAUGGUUUAGGCGUGCUUAAGGCUGCGGUUGCAGCGAUUAAAAAAGCCGCAGCAAAUUAAUUUAUUUUAAAUUUUCCCGUCAAAGAGUUAUUUCAUAAAUCAAUACCGCAAUAUUUAAAUUGCGGUUUUUAAGGGUAUUUUUCUAUGAGUAAUGUUAUUGCAUCGCUUGAAAAGGUACUCCUCCCUUUUGCAGUUAAAAUAGGAAAGCAGCCACACGUUAAUGCAAUCAAAAAUGGCUUUAUUC" --interaction_bps "(134,56):(135,55):(136,54):(137,53):(138,52):(139,51):(140,50):(141,49):(142,48):(143,47):(146,44):(147,43):(148,42):(149,41):(150,40):(151,39):(152,38)" --seed 4  --write_states "states.tsv"  --write_structures test.fa --write_all_barriers "barriers.tsv" --fixed_intramolecular_structures --str_b "......(((((...)))))............(((((((((.................(((((((....))))))))))))))))"
``

### dependencies
- g++
- lboost
- lIntaRNA
- fopenmp
- lRNA
- leasylogging


### help

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
