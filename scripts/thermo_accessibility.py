#!/usr/bin/env  python3
#
import argparse
import re
import sys
import RNA


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


def accessibility(sequence, outfile, identifier):
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
        line_list = [identifier, "equilibrium", str(l)]
        line_list += [ "{:g}".format(q[i]) for i in range(1, l + 1) ]
        line_list += [ "NA" for i in range(l + 1, len(sequence) + 1) ]

        # print accessibilities
        print(",".join(line_list), file=outfile)

    if outfile != sys.stdout:
        outfile.close()

def main():
    outfile   = None
    sequence  = None
    n         = 0
    parser    = argparse.ArgumentParser()
    group     = parser.add_mutually_exclusive_group()

    parser.add_argument("-i", "--input",
                        type=str,
                        help="input file (FASTA formatted)")
    group.add_argument("-o", "--output",
                        type=str,
                        help="output file name. Defaults to print to stdout",
                        default = None)
    group.add_argument("-a", "--append-to",
                        type = str,
                        help="Append output to an existing file instead of overwriting it")
    parser.add_argument("--header",
                        help="Add header line", action="store_true")
    parser.add_argument("--no-header",
                        action = "store_true",
                        help = "Do not add header line if using -o/--output option")
    parser.add_argument("-s", "--sequence-id",
                        type=str, help="Sequence identifier", default=None)

    args = parser.parse_args()

    if not args.input:
        exit(1)

    if args.output:
        outfile = open(args.output, "w")
        if not args.no_header:
            args.header = True
    elif args.append_to:
        outfile = open(args.append_to, "a")

    if not outfile:
        outfile = sys.stdout

    if args.input:
        sequence, seq_id  = get_sequence_line(args.input)
        n                 = len(sequence)

    # exit script if no sequence is available
    if not sequence:
        exit(1)

    # print header line
    if args.header:
        head_list = ["sequence", "method", "length"]
        head_list += [str(i) for i in range(1, n + 1) ]
        print(",".join(head_list), file=outfile)

    if not args.sequence_id:
        args.sequence_id = seq_id if seq_id else "RNA"

    accessibility(sequence, outfile, args.sequence_id)

if __name__ == '__main__':
    main()
