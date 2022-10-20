#!/usr/bin/env python

import re
import sys
import math
import pandas as pd
import argparse

def drtrafo_get_drforna_energies(f, bins):
    last_time = 0
    last_step = 0
    last_bin = {}

    for e, line in enumerate(f):
        if e == 0:
            assert line == "id time occupancy structure energy\n"
            continue

        [_, time, occ, ss, en] = line.split()
        occ = min(int(round(float(occ)*10000)), 10000)
        assert 0 <= occ <= 10000
        en = float(en)
        # assert en <= 0

        if occ == 0 : 
            continue

        step = len(ss)
        if time == last_time and step == last_step:
            # let's add occuancy to last-time bin and continue
            if en in last_bin:
                last_bin[en] += occ
            else:
                last_bin[en] = occ
            continue

        elif step > last_step:
            # yes, let's push the last bin
            bins[last_step] = last_bin.copy()

        # else, initialize a new bin
        last_time = time
        last_step = step
        last_bin = {}
        last_bin[en] = occ

    bins[last_step] = last_bin.copy()
    return

def get_uprobs(drf):
    uprobs = []
    with open(drf) as f:
        llen, ltime, lltime = 0, '0', '0'
        up = []
        for i, line in enumerate(f):
            if i == 0:
                assert line == "id time occupancy structure energy\n"
                continue
            [_, stime, occ, ss, en] = line.split()
            occu = float(occ)
            if ltime == stime and len(ss) == llen:
                for j, b in enumerate(ss):
                     if b == '.': 
                         up[j] += occu
            else:
                if len(ss) > llen:
                    uprobs.append(up)
                up = [0 for _ in range(len(ss))]
                for j, b in enumerate(ss):
                     if b == '.': 
                         up[j] += occu
            llen = len(ss)
            ltime = stime
        uprobs.append(up)
    return uprobs

def main():
    """ For example:
    ./drtrafo_parser.py -i SRP.drf --mode energyrange
    ./drtrafo_parser.py -i SRP.drf --mode accessibility
    """
    outfile = None
    parser  = argparse.ArgumentParser()
    group   = parser.add_mutually_exclusive_group()

    parser.add_argument("-i", "--input",
                        required = True,
                        type = str,
                        help = "Input file name")
    group.add_argument("-o", "--output",
                        type=str,
                        help="output file name. Defaults to print to stdout. Automatically adds header unless (--no-header) is specified",
                        default = None)
    group.add_argument("-a", "--append-to",
                        type = str,
                        help="Append output to an existing file instead of overwriting it")
    parser.add_argument("-n", "--length",
                        type = int,
                        help = ("Maximum number of nucleotides, i.e. number of transcription steps"),
                        default = 117)
    parser.add_argument("-m", "--mode",
                        choices = ('energyrange', 'energy', 'shape', 'lshape','up'),
                        help = "Output mode, i.e. output file format",
                        default = 'energyrange')
    parser.add_argument("--header",
                        action = "store_true",
                        help = "Add header line") 
    parser.add_argument("--no-header",
                        action = "store_true",
                        help = "Do not add header line if using -o/--output option")
    parser.add_argument("-s", "--sequence-id",
                        type = str,
                        help = "Sequence identifier",
                        default="RNA")
    parser.add_argument("-t", "--tool-id",
                        type = str,
                        help = "Tool identifier",
                        default="DrTrafo")

    args = parser.parse_args()

    if args.output:
        outfile = open(args.output, "w")
        if not args.no_header:
            args.header = True
    elif args.append_to:
        outfile = open(args.append_to, "a")

    if not outfile:
        outfile = sys.stdout

    length  = args.length
    mode    = args.mode
    outlen = length + 1

    if mode in ['shape', 'lshape','up']:
        uprobs = get_uprobs(args.input)
        if args.header:
            header_list = ["length", "method", "name"]
            header_list += [ map(str,range(1, outlen)) ]

            print(",".join(header_list), file=outfile)

        if mode == 'lshape':
            lprobs = uprobs[-1]

            data_list = [f'{len(lprobs)},{args.tool_id},{args.sequence_id}']
            data_list += [f'{round(p,2)}' for p in lprobs ]
            data_list += ['NA' for _ in range(length - len(lprobs))]
            print(','.join(data_list), file=outfile)

        else:
            for s in range(1, outlen):
                data_list = [f'{s},{args.tool_id},{args.sequence_id}']
                data_list += [str(round(p,2)) for p in uprobs[s]]
                data_list += ['NA' for _ in range(length - s)]
                print(','.join(data_list), file=outfile)

    elif mode in ['energyrange', 'energy']:
        f = open(args.input)
        # Initialize a dict for every transcription step.
        bins = [ dict() for i in range(outlen) ]
        drtrafo_get_drforna_energies(f, bins)

        if args.header:
            header_list = [ "length",
                            "method",
                            "name",
                            "Q25",
                            "Q75",
                            "Qmedian",
                            "Qmean",
                            "Qmin",
                            "Qmax"]
            print(",".join(header_list), file=outfile)

        for s in range(1, outlen):
            # Expand energy counts to list
            energies = []
            for key in bins[s]:
                energies = energies + [float(key) for _ in range(0, bins[s][key])]

            df  = pd.Series(energies)
            qt  = df.quantile([0.25,0.75])
            data_list = [f'{s:d}',
                         f'{args.tool_id}',
                         f'{args.sequence_id}',
                         f'{qt[0.25]:.2f}',
                         f'{qt[0.75]:.2f}',
                         f'{df.median():.2f}',
                         f'{df.mean():.2f}',
                         f'{df.min():.2f}',
                         f'{df.max():.2f}']
            print(",".join(data_list), file=outfile)

        f.close()

    return

if __name__ == '__main__':
    main()
