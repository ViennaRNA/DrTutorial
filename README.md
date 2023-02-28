# A guide to computational cotranscriptional folding using the SRP RNA - Data and Scripts

This repository contains data, scripts and other resources for the DrTutorial book chapter

**"A guide to computational cotranscriptional folding using the SRP RNA"**

Below is a list of resources collected within this repository together with their references and sources (if applicable).

## Installation

The resources can be installed globally via `pip`:

```
pip install .
```

or locally in the users home directory via:

```
pip install --user .
```


## FASTA formatted RNA sequence files

The `sequences/` directory contains the following FASTA-formatted
sequences:

| sequence id | local file | data set (sequence name) | reference |
| ----------- | ---------- | ------------------------ | --------- |
| SRPn | [SRPn.fa](sequences/SRPn.fa) | Native E.coli SRP sequence | [1], [2], [3], [4] |
| SRPt | [SRPt.fa](sequences/SRPt.fa) | Trapped E.coli SRP mutant (U21C) | [1], [3] |
| SRPr | [SRPr.fa](sequences/SRPr.fa) | Rescue E.coli SRP mutant (U21C/C22U/G93A) | [1] |
| SRPf | [SRPf.fa](sequences/SRPf.fa) | Folding-path E.coli SRP mutant (U35C/U37C) | [3] |

and additional structural motifs of transient helices in
FASTA format with dot-bracket structures as reported in
literature and within this tutorial:

| sequence id | local file | description | reference |
| ----------- | ---------- | ----------- | --------- |
| SRPn | [SRPn-H1.fa](sequences/SRPn-H1.fa) | Transient helix H1 | [1],[3] |
| SRPn | [SRPn-H1a.fa](sequences/SRPn-H1a.fa) | Transient helix H1a | [4] |
| SRPn | [SRPn-H1b.fa](sequences/SRPn-H1b.fa) | Transient helix H1b | [3] |
| SRPn | [SRPn-H1c.fa](sequences/SRPn-H1c.fa) | Transient helix H1c |  |
| SRPn | [SRPn-H1-H2-R1-H3.fa](sequences/SRPn-H1-H2-R1-H3.fa) | Transient helices H1, H2, R1, and H3 | [1],[3] |
| SRPn | [SRPn-H1c-H2a-H3a.fa](sequences/SRPn-H1c-H2a-H3a.fa) | Transient helices H1c, H2a, and H3a | [1] |
| SRPn | [SRPn-H2-R1-R2.fa](sequences/SRPn-H2-R1-R2.fa) | Functional conformation with helices H2, R1, and R2 | [1],[2],[3] |


## SHAPE reactivity data

The `SHAPE/` directory contains co-transcriptional SHAPE data provided by [1] and [2] and downloaded from
the [RMDB database](https://rmdb.stanford.edu/).

| sequence id | local file | data set (with link to RMDB) | reference |
| ----------- | ---------- | ---------------------------- | --------- |
| SRPn | [SRPECLI_BZCN_0001.rdat](SHAPE/SRPECLI_BZCN_0001.rdat)  |  [Replicate 1 for native E.coli SRP sequence](https://rmdb.stanford.edu/detail/SRPECLI_BZCN_0001) | [2] |
| SRPn | [SRPECLI_BZCN_0002.rdat](SHAPE/SRPECLI_BZCN_0002.rdat)  |  [Replicate 2 for native E.coli SRP sequence](https://rmdb.stanford.edu/detail/SRPECLI_BZCN_0002) | [2] |
| SRPn | [SRPECLI_BZCN_0003.rdat](SHAPE/SRPECLI_BZCN_0003.rdat)  |  [Replicate 3 for native E.coli SRP sequence](https://rmdb.stanford.edu/detail/SRPECLI_BZCN_0003) | [2] |
| SRPn | [SRPECLI_BZCN_0004.rdat](SHAPE/SRPECLI_BZCN_0004.rdat)  |  [Equilibrium refold for native E.coli SRP sequence](https://rmdb.stanford.edu/detail/SRPECLI_BZCN_0004) | [2] |
| SRPt | [SRPU21C_BZCN_0001.rdat](SHAPE/SRPU21C_BZCN_0001.rdat)  |  [Replicate 1 for trapped E.coli SRP mutant (U21C)](https://rmdb.stanford.edu/detail/SRPU21C_BZCN_0001) | [1] |
| SRPt | [SRPU21C_BZCN_0002.rdat](SHAPE/SRPU21C_BZCN_0002.rdat)  |  [Replicate 2 for trapped E.coli SRP mutant (U21C)](https://rmdb.stanford.edu/detail/SRPU21C_BZCN_0002) | [1] |
| SRPt | [SRPU21C_BZCN_0003.rdat](SHAPE/SRPU21C_BZCN_0003.rdat)  |  [Replicate 3 for trapped E.coli SRP mutant (U21C)](https://rmdb.stanford.edu/detail/SRPU21C_BZCN_0003) | [1] |
| SRPt | [SRPU21C_BZCN_0004.rdat](SHAPE/SRPU21C_BZCN_0004.rdat)  |  [Equilibrium refold 1 for trapped E.coli SRP mutant (U21C)](https://rmdb.stanford.edu/detail/SRPU21C_BZCN_0004) | [1] |
| SRPt | [SRPU21C_BZCN_0005.rdat](SHAPE/SRPU21C_BZCN_0005.rdat)  |  [Equilibrium refold 2 for trapped E.coli SRP mutant (U21C)](https://rmdb.stanford.edu/detail/SRPU21C_BZCN_0005) | [1] |
| SRPr | [SRP21CR_BZCN_0001.rdat](SHAPE/SRP21CR_BZCN_0001.rdat)  |  [Replicate 1 for rescue E.coli SRP mutant (U21C/C22U/G93A)](https://rmdb.stanford.edu/detail/SRP21CR_BZCN_0001) | [1] |
| SRPr | [SRP21CR_BZCN_0002.rdat](SHAPE/SRP21CR_BZCN_0002.rdat)  |  [Replicate 2 for rescue E.coli SRP mutant (U21C/C22U/G93A)](https://rmdb.stanford.edu/detail/SRP21CR_BZCN_0002) | [1] |
| SRPr | [SRP21CR_BZCN_0003.rdat](SHAPE/SRP21CR_BZCN_0003.rdat)  |  [Replicate 3 for rescue E.coli SRP mutant (U21C/C22U/G93A)](https://rmdb.stanford.edu/detail/SRP21CR_BZCN_0003) | [1] |


## Python and R scripts

The `scripts/` directory contains various prediction-, parser- and plotting scripts.

| script name | purpose |
| ----------- | ------- |
| [thermo_predict.py](scripts/thermo_predict.py) | Predict thermodynamic equilibrium profiles for energies and accessibilities of nascent transcripts |
| [drf_parser.py](scripts/drf_parser.py) | Parser for *.drf output as produced by the `DrTransformer` and `DrKinfold` programs |
| [convert_rdat.py](scripts/convert_rdat.py) | Convert cotranscriptional SHAPE reactivity data from RMDBs .rdat files into the CSV format produced by `drf_parser.py` and `thermo_predict.py` |
| [plot_energy_bands.R](scripts/plot_energy_bands.R) | Produce an energy distribution plot for cotranscriptionally formed structures |
| [plot_accessibility.R](scripts/plot_accessibility.R) | Plot accessibility profiles for nascent transcripts |
| [make_SRP_images.sh](scripts/make_SRP_images.sh) | Create annotated secondary structure plots of the transient helix motifs |

Additionally, in the `drconverters/` directory, this repository contains a snapshot of the
[`drconverters`](https://github.com/bad-ants-fleet/drconverters) script package. This package
will be automagically included in the installation process.

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
- [4] Wong, T. N., Sosnick, T. R., and Pan, T. (2007). "Folding of noncoding
RNAs during transcription facilitated by pausing-induced nonnative structures.",
Proceedings of the National Academy of Sciences, 104(46), 17995-18000.
