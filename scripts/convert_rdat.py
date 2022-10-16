#!/usr/bin/env python3
#

import sys
import re
import argparse


def rdat2csv(args, outfile):
    annot_pat = re.compile(r"^DATA_ANNOTATION:(\d+).*datatype:REACTIVITY.*ID:Length(\d+)")
    data_pat  = re.compile(r"^DATA:(\d+)\s+(.*)$")

    with open(args.input) as f:
        data = dict()
        min_l = 10000
        max_l = 0

        for line in f:
            m = annot_pat.match(line)
            if m:
                n = "DATA." + m.group(1)
                l = m.group(2)
                data[n] = {'length' : l, 'reactivities' : []}

                continue

            m = data_pat.match(line)
            if m:
                n = "DATA." + m.group(1)
                if n in data:
                    data[n]['reactivities'] = [float(rr) for rr in m.group(2).rstrip().split("\t") ]

                continue

    # re-order data
    reactivities = { dd['length'] : dd['reactivities'] for dd in data.values() }

    # determine the maximum length of the data
    max_l = max([int(n) for n in reactivities])

    # print header line
    if args.header:
        head_list = ["length", "method", "name"]
        head_list += [str(i) for i in range(1, max_l + 1) ]
        print(",".join(head_list), file=outfile)

    for l in range(1, max_l + 1):
        if str(l) in reactivities:
            remaining = max_l - l
            line_list = [str(l), "SHAPE", args.sequence_id]
            line_list += [ str(r) for r in reactivities[str(l)]]
            line_list += [ 'NA' for _ in range(remaining) ]
            print(",".join(line_list), file = outfile)
        else:
            print(",".join([str(l), "SHAPE", args.sequence_id] + ['NA' for _ in range(max_l) ]), file = outfile)


def main():
    outfile       = None
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
                        help = "Sequence identifier",
                        default = "RNA")

    args = parser.parse_args()

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

    # call prediction mode function
    rdat2csv(args, outfile)


if __name__ == '__main__':
    main()
