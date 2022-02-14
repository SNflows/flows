# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import argparse
from flows.plots import plots_interactive
import flows

if __name__ == '__main__':
    # Parse command line arguments:
    parser = argparse.ArgumentParser(description='Run photometry pipeline.')
    parser.add_argument('-t', '--target', type=str, help='TIC identifier of target.', nargs='?', default=2)
    parser.add_argument('-s', '--site', type=int, help='TIC identifier of target.', nargs='?', default=None)
    parser.add_argument('-d', '--date', type=str, help='TIC identifier of target.', nargs='?', default=None)
    parser.add_argument('-o', '--output', type=str, help='TIC identifier of target.', nargs='?', default=None)
    args = parser.parse_args()

    if args.output is None:
        plots_interactive()

    flows.visibility(target=args.target, siteid=args.site, date=args.date, output=args.output)
