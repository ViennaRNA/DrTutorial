# DrTutorial Data files and scripts

This repository contains data, scripts and other resources for the DrTransformer Tutorial.


## FASTA formatted RNA sequence files

The `sequences/` directory contains the following FASTA-formatted
sequences

- [SRPn.fa](sequences/SRPn.fa) ... The native E.coli SRP sequence as used in [1], [2] and [3]
- [SRPt.fa](sequences/SRPt.fa) ... The trapped E.coli SRP mutant (U21C) as used in [1] and [3]
- [SRPr.fa](sequences/SRPr.fa) ... The rescue E.coli SRP mutant (U21C/C22U/G93A) as used in [1]
- [SRPf.fa](sequences/SRPf.fa) ... The folding-path E.coli SRP mutant (U35C/U37C) as used in [3]


## SHAPE reactivity data

The `SHAPE/` directory contains co-transcriptional SHAPE data provided by [1] and [2] and downloaded from
the [RMDB database](https://rmdb.stanford.edu/).

- [SRPECLI_BZCN_0001.rdat](SHAPE/SRPECLI_BZCN_0001.rdat) ... Replicate 1 for native E.coli SRP sequence from [2]
- [SRPECLI_BZCN_0002.rdat](SHAPE/SRPECLI_BZCN_0002.rdat) ... Replicate 2 for native E.coli SRP sequence from [2]
- [SRPECLI_BZCN_0003.rdat](SHAPE/SRPECLI_BZCN_0003.rdat) ... Replicate 3 for native E.coli SRP sequence from [2]
- [SRPECLI_BZCN_0004.rdat](SHAPE/SRPECLI_BZCN_0004.rdat) ... Equilibrium refold for native E.coli SRP sequence from [2]
- [SRPU21C_BZCN_0001.rdat](SHAPE/SRPU21C_BZCN_0001.rdat) ... Replicate 1 for trapped E.coli SRP mutant (U21C) from [1]
- [SRPU21C_BZCN_0002.rdat](SHAPE/SRPU21C_BZCN_0002.rdat) ... Replicate 2 for trapped E.coli SRP mutant (U21C) from [1]
- [SRPU21C_BZCN_0003.rdat](SHAPE/SRPU21C_BZCN_0003.rdat) ... Replicate 3 for trapped E.coli SRP mutant (U21C) from [1]
- [SRPU21C_BZCN_0004.rdat](SHAPE/SRPU21C_BZCN_0004.rdat) ... Equilibrium refold 1 for trapped E.coli SRP mutant (U21C) from [1]
- [SRPU21C_BZCN_0005.rdat](SHAPE/SRPU21C_BZCN_0005.rdat) ... Equilibrium refold 2 for trapped E.coli SRP mutant (U21C) from [1]
- [SRP21CR_BZCN_0001.rdat](SHAPE/SRP21CR_BZCN_0001.rdat) ... Replicate 1 for rescue E.coli SRP mutant (U21C/C22U/G93A) from [1]
- [SRP21CR_BZCN_0002.rdat](SHAPE/SRP21CR_BZCN_0002.rdat) ... Replicate 2 for rescue E.coli SRP mutant (U21C/C22U/G93A) from [1]
- [SRP21CR_BZCN_0003.rdat](SHAPE/SRP21CR_BZCN_0003.rdat) ... Replicate 3 for rescue E.coli SRP mutant (U21C/C22U/G93A) from [1]

Each of these `.rdat` files is accompanied by a corresponding `.dat` file that
contains the same reactivity data in comma-separated value (CSV) format.


## Python and R scripts

the `scripts/` directory contains various prediction-, parser- and plotting scripts.

- [thermo_energy.py](scripts/thermo_energy.py) ... 
- [thermo_accessibility.py](scripts/thermo_accessibility.py) ... 
- [drf_parser.py](scripts/drf_parser.py) ... 
- [plot_accessibility.R](scripts/plot_accessibility.R) ... 


## References

- [1] Yu, A. M., Gasper, P. M., Cheng, L., Lai, L. B., Kaur, S., Gopalan, V.,
Chen, A. A., and Lucks, J. B. (2021). "Computationally reconstructing
cotranscriptional RNA folding from experimental data reveals rearrangement
of non-native folding intermediates.", Molecular cell , 81(4), 870-883.
- [2] Watters, K. E., Strobel, E. J., Yu, A. M., Lis, J. T., and Lucks, J. B.
(2016). "Cotranscriptional folding of a riboswitch at nucleotide resolution.",
Nature structural & molecular biology, 23(12), 1124.
- [3] Fukuda, S., Yan, S., Komi, Y., Sun, M., Gabizon, R., and Bustamante, C.
(2020). "The biogenesis of SRP RNA is modulated by an RNA folding intermediate
attained during transcription.", Molecular Cell , 77(2), 241-250.
