#!/usr/bin/python
#
import csv
import argparse
import re
import RNA
import pandas as pd


def get_sequence_line(filename):
    """
    Read a file and return the first line that
    consists of RNA/DNA sequence letters. If no
    such line is found or there is any problem
    with the input file, return None
    """
    sequence = None
    header   = None

    seq_pattern       = re.compile(r"([ACGUTNacgutn]+)")
    fasta_header_pat  = re.compile(r"^>(.*)$")

    with open(filename) as f:
        for line in f:
            m = fasta_header_pat.match(line)
            if m:
                # only process first FASTA entry in file
                if header and sequence:
                    return (sequence, header)

                header    = m.group(1)
                sequence  = ""
                continue

            m = seq_pattern.match(line)
            if m:
                # for input without FASTA header, each
                # sequence must be on a single line
                if not header:
                    return (m.group(1), None)
                elif sequence:
                    sequence += m.group(1)
                else:
                    sequence = m.group(1)

    return (sequence, header)


def get_SHAPE_data(filename, n, offset = 14):
    """
    Read csv separated SHAPE data from a file and return it
    as list of lists of SAHPE reactivities
    """
    SHAPE_data  = []

    with open(filename, "r") as f:
        reader = csv.reader(f)
        next(reader, None) # skip header
        for row in reader:
            i   = int(row[0])
            dat = [-999.] + [ float(d) if d != "NA" else -999. for d in row[1:] ]
            if i > offset:
                if i - offset - 1 > len(SHAPE_data):
                    SHAPE_data += [ [-999.0 for j in range(0, n + 1) ] for k in range(i - offset - len(SHAPE_data)) ]
                SHAPE_data.append(dat)

    return SHAPE_data



def fold_and_print(sequence, options, SHAPE_data, outfile):
    """
    Compute MFE and obtain Boltzmann samples from all nascent transcript lengths.
    Additionally, guide structure prediction by SHAPE data if available.
    """
    md          = RNA.md()
    md.uniq_ML  = 1

    n       = len(sequence)
    fc      = RNA.fold_compound(sequence, md)
    ss, mfe = fc.mfe()

    fc.exp_params_rescale(mfe)
    fc.pf()

    for i in range(options.start, n + 1):
        energies  = []
        subseq    = sequence[0 : i]
        (ss, mfe) = RNA.fold(subseq)

        if SHAPE_data:
            fc_sub = RNA.fold_compound(subseq, md)
            if i < len(SHAPE_data):
                fc_sub.sc_add_SHAPE_deigan(SHAPE_data[i], 1.1, -0.3)
            fc_sub.exp_params_rescale(mfe)
            fc_sub.pf()
        else:
            fc_sub = fc

        for s in fc_sub.pbacktrack5(options.samples, i):
            energies.append(RNA.eval_structure_simple(subseq, s))

        df  = pd.Series(energies)
        qt  = df.quantile([0.25,0.75])

        # print result for sampling approach
        line = [str(i), "sampling", options.sequence_id]
        line += ["{:6.2f}".format(d) for d in [ qt[0.25],
                                                qt[0.75],
                                                df.median(),
                                                df.mean(),
                                                df.min(),
                                                df.max() ] ]
        print("\t".join(line), file=outfile)

        if options.mfe:
            line = [str(i), "MFE", options.sequence_id]
            line += ["{:6.2f}".format(d) for d in [mfe for i in range(6)] ]

            print("\t".join(line), file=outfile)
        

def main():
    outfile     = None
    n           = 0
    SHAPE_data  = None
    sequence    = None
    parser      = argparse.ArgumentParser()
    group       = parser.add_mutually_exclusive_group()

    parser.add_argument("-i", "--input", type=str, help="input file")
    group.add_argument("-o", "--output",
                        type=str,
                        help="output file name. Defaults to print to stdout",
                        default = None)
    group.add_argument("-a", "--append-to",
                        type = str,
                        help="Append output to an existing file instead of overwriting it")
    parser.add_argument("-n", "--samples", type=int, help="Number of samples per subsequence", default = 1000)
    parser.add_argument("-m", "--mfe", help="Add MFE values", action="store_true")
    parser.add_argument("--header", help="Add header line", action="store_true")
    parser.add_argument("-s", "--sequence-id", type=str, help="overwrite sequence identifier", default = None)
    parser.add_argument("--SHAPE", type=str, help="SHAPE file")
    parser.add_argument("--start", type=int, help="Start length", default = 15)
    parser.add_argument("--offset", type=int, help="Offset of SHAPE data", default = 14)
    parser.add_argument("-P", "--params", type=str, help="Energy parameters")


    args = parser.parse_args()

    # read input sequence
    if args.input:
        sequence, seq_id  = get_sequence_line(args.input)
        n                 = len(sequence)

    # exit script if no sequence is available
    if not sequence:
        exit(1)

    if args.output:
        outfile = open(args.output, "w")
    elif args.append_to:
        outfile = open(args.append_to, "a")

    if not outfile:
        outfile = sys.stdout

    # parse further arguments
    if args.SHAPE:
        SHAPE_data = get_SHAPE_data(args.SHAPE, n, args.offset)

    if args.params:
        RNA.read_parameter_file(args.params)

    if args.header:
        header = ["step",
                  "method",
                  "sequence",
                  "quant25",
                  "quant75",
                  "median",
                  "mean",
                  "min",
                  "max"]

        print("\t".join(header), file=outfile)

    if not args.sequence_id:
        args.sequence_id = seq_id if seq_id else "sequence"


    fold_and_print(sequence, args, SHAPE_data, outfile)



if __name__ == '__main__':
    main()
