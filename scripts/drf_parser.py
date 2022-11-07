#!/usr/bin/env python

import re
import sys
import math
import pandas as pd
import argparse

def drtrafo_get_drforna_energies(drf):
    bins = []
    with open(drf) as f:
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
                bins.append(last_bin.copy())

            # else, initialize a new bin
            last_time = time
            last_step = step
            last_bin = {}
            last_bin[en] = occ
        # let's push the last bin
        bins.append(last_bin.copy())
    return bins

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

    group.add_argument("-o", "--output", type = str,
                        help = "Output file name. Defaults to STDOUT. Automatically adds header unless (--no-header) is specified")
    group.add_argument("-a", "--append", type = str,
                        help = "Append output to an existing file. Does not print the header, disables STDOUT.")
    parser.add_argument("-i", "--input", required = True, type = str,
                        help = "Input file name. Defaults to STDIN.")
    parser.add_argument("--energy", action = "store_true",
                        help = "Calculate accessibilities.")
    parser.add_argument("--accessibility", action = "store_true",
                        help = "Calculate accessibilities.")
    parser.add_argument("--access-by-index", type = int, default = 0,
                        help = "Calculate accessibilities.")
    parser.add_argument("-l", "--length", type = int,
                        help = "Specify transcript length (e.g. when data is missing).")
    parser.add_argument("-n", "--name", type = str, required = True,
                        help = "Sequence name")
    parser.add_argument("-m", "--method", type = str, required = True,
                        help = "Method name")
    parser.add_argument("--lshape", type = int,
                        help = "Shortcut to read results of a shapeseq simulation.")
    args = parser.parse_args()

    outfile = sys.stdout
    if args.output:
        outfile = open(args.output, "w")
    elif args.append:
        outfile = open(args.append, "a")

    if args.lshape:
        assert not args.length
        assert not args.access_by_index
        args.length = args.lshape
        args.access_by_index = -1

    if args.accessibility or args.access_by_index:
        uprobs = get_uprobs(args.input) # uprobs[0] = []
        datalen = args.length + 1 if args.length else len(uprobs)
        if args.output:
            print(f"length,method,name,{','.join(map(str, range(1, datalen)))}", file = outfile)
        if args.accessibility:
            for l in range(1, datalen):
                data = [str(round(p, 2)) for p in uprobs[l]] + ['NA' for _ in range(datalen - 1 - l)]
                print(f"{l} {args.method} {args.name} {','.join(data)}", file = outfile)
        if args.access_by_index:
            idx = args.access_by_index % len(uprobs)
            data = [str(round(p, 2)) for p in uprobs[idx]] + ['NA' for _ in range(datalen - 1 - idx)]
            print(f"{len(uprobs)-1} {args.method} {args.name} {','.join(data)}", file = outfile)

    if args.energy:
        eranges = drtrafo_get_drforna_energies(args.input)
        datalen = args.length + 1 if args.length else len(eranges)
        if args.output:
            header_list = ["length",
                           "method",
                           "name",
                           "Q25",
                           "Q75",
                           "Qmedian",
                           "Qmean",
                           "Qmin",
                           "Qmax"]
            print(",".join(header_list), file = outfile)

        for l in range(1, datalen):
            # Expand energy counts to list
            energies = []
            for key in eranges[l]:
                energies = energies + [float(key) for _ in range(0, eranges[l][key])]
            df  = pd.Series(energies)
            qt  = df.quantile([0.25,0.75])
            data_list = [f'{l:d}',
                         f'{args.method}',
                         f'{args.name}',
                         f'{qt[0.25]:.2f}',
                         f'{qt[0.75]:.2f}',
                         f'{df.median():.2f}',
                         f'{df.mean():.2f}',
                         f'{df.min():.2f}',
                         f'{df.max():.2f}']
            print(",".join(data_list), file=outfile)

    if args.output or args.append:
        outfile.close()

if __name__ == '__main__':
    main()
