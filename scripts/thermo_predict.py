#!/usr/bin/python
#
import sys
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
    fasta_header_pat  = re.compile(r"^>\s*([^\s]+)")

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


def accessibility(args, sequence, outfile):
    """
    Predict accessibility profiles
    """
    # print header line
    if args.header:
        head_list = ["length", "method", "name"]
        head_list += [str(i) for i in range(1, n + 1) ]
        print(",".join(head_list), file=outfile)

    # loop over all nascent transcripts
    for l in range(1, len(sequence) + 1):
        # create fold_compound for subsequence
        fc  = RNA.fold_compound(sequence[0:l])
        # compute MFE
        (ss, mfe) = fc.mfe()
        # rescale Boltzmann factors
        fc.exp_params_rescale(mfe)
        # compute partition function and base pair probabilities
        fc.pf()
        # retrieve base pair probabilities
        bpp = fc.bpp()
        # initialize list for accessibilities
        q = [ 0 for i in range(0, l + 1) ]
        # sum-up probabilities to be paired
        for i in range(1, l + 1):
            for j in range(i, l + 1):
                q[i] += bpp[i][j]
                q[j] += bpp[i][j]
        # turn probabilities to be paired into actual accessibilities
        for i in range(1, l + 1):
            q[i] = 1 - q[i]
        # collect data for current line of accessibilities
        line_list = [str(l), "equilibrium", args.sequence_id]
        line_list += [ "{:g}".format(q[i]) for i in range(1, l + 1) ]
        line_list += [ "NA" for i in range(l + 1, len(sequence) + 1) ]
        # print accessibilities
        print(",".join(line_list), file=outfile)

    if outfile != sys.stdout:
        outfile.close()

def fold_and_print(args, sequence, outfile):
    """
    Compute MFE and obtain Boltzmann samples from all nascent transcript lengths.
    Additionally, guide structure prediction by SHAPE data if available.
    """
    SHAPE_data = None

    if args.SHAPE:
        SHAPE_data = get_SHAPE_data(args.SHAPE, n, args.offset)

    if args.header:
        header = ["length",
                  "method",
                  "name",
                  "Q25",
                  "Q75",
                  "Qmedian",
                  "Qmean",
                  "Qmin",
                  "Qmax"]

        print(",".join(header), file=outfile)

    md          = RNA.md()
    md.uniq_ML  = 1

    n       = len(sequence)
    fc      = RNA.fold_compound(sequence, md)
    ss, mfe = fc.mfe()

    fc.exp_params_rescale(mfe)
    fc.pf()

    for i in range(args.start, n + 1):
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

        for s in fc_sub.pbacktrack5(args.samples, i):
            energies.append(RNA.eval_structure_simple(subseq, s))

        df  = pd.Series(energies)
        qt  = df.quantile([0.25,0.75])

        # print result for sampling approach
        line = [str(i), "sampling", args.sequence_id]
        line += ["{:.2f}".format(d) for d in [ qt[0.25],
                                                qt[0.75],
                                                df.median(),
                                                df.mean(),
                                                df.min(),
                                                df.max() ] ]
        print(",".join(line), file=outfile)

        if args.mfe:
            line = [str(i), "MFE", args.sequence_id]
            line += ["{:.2f}".format(d) for d in [mfe for i in range(6)] ]

            print(",".join(line), file=outfile)

def main():
    outfile       = None
    n             = 0
    sequence      = None
    parser        = argparse.ArgumentParser()
    group_output  = parser.add_mutually_exclusive_group()
    group_header  = parser.add_mutually_exclusive_group()

    parser.add_argument("-i",
                        "--input",
                        type = str,
                        help = "Sequence input file, e.g. FASTA formatted.",
                        required = True)
    group_output.add_argument("-o", "--output",
                        type = str,
                        help = "Output file name. Defaults to print to stdout.",
                        default = None)
    group_output.add_argument("-a", "--append-to",
                        type = str,
                        help = "Append output to an existing file instead of overwriting it.")
    group_header.add_argument("--header",
                        action="store_true",
                        help="Add header line")
    group_header.add_argument("--no-header",
                        action = "store_true",
                        help = "Do not add header line if using -o/--output option or when printing to stdout.")
    parser.add_argument("-s", "--sequence-id",
                        type = str,
                        help = "Overwrite sequence identifier. Defaults to extract identifier from FASTA header.",
                        default = None)
    parser.add_argument("-P", "--params",
                        type = str,
                        help = "Load a different energy parameter set.")

    # create sub-parsers for the different modes of this script
    sub_parsers = parser.add_subparsers(title = 'subcommands',
                                        description = 'valid sub-commands',
                                        required = True)

    # options for the 'energy distribution' mode
    parser_en = sub_parsers.add_parser('energy',
                                       help='Energy distribution help')
    parser_en.add_argument("-n", "--samples",
                           type = int,
                           help = "Number of samples per subsequence.",
                           default = 1000)
    parser_en.add_argument("--mfe",
                           action = "store_true",
                           help = "Add MFE values.")
    parser_en.add_argument("--SHAPE",
                          type = str,
                          help = "cotranscriptional SHAPE reactivity file (csv formatted).")
    parser_en.add_argument("--start",
                           type = int,
                           help = "Start length",
                           default = 15)
    parser_en.add_argument("--offset",
                           type = int,
                           help = "Offset of SHAPE data",
                           default = 14)
    parser_en.set_defaults(func = fold_and_print)

    # options for the 'accessibility profile' mode
    parser_up = sub_parsers.add_parser('accessibility',
                                       help = 'Accessibility profile help')
    # no further options for this mode (yet)
    parser_up.set_defaults(func = accessibility)

    args = parser.parse_args()

    # read input sequence
    sequence, seq_id  = get_sequence_line(args.input)
    n                 = len(sequence)

    # exit script if no sequence is available
    if not sequence:
        print(f'Unable to parse any sequence data from file {args.input}')
        exit(1)

    # prepare output stream
    if args.output:
        outfile = open(args.output, "w")
        if not args.no_header:
            args.header = True
    elif args.append_to:
        outfile = open(args.append_to, "a")

    if not outfile:
        if not args.no_header:
            args.header = True
        outfile = sys.stdout

    # load energy parameters if necessary
    if args.params:
        RNA.read_parameter_file(args.params)

    # prepare sequence identifier
    if not args.sequence_id:
        args.sequence_id = seq_id if seq_id else "RNA"

    # call prediction mode function
    args.func(args, sequence, outfile)


if __name__ == '__main__':
    main()
