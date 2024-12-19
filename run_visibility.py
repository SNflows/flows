# -*- coding: utf-8 -*-
"""

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

import argparse
from flows.plots import plots_interactive
import flows
from flows.target import Target
import re
import os.path

if __name__ == '__main__':
    # Parse command line arguments:
    parser = argparse.ArgumentParser(description='Run photometry pipeline.')
    parser.add_argument('-t', '--target', type=str, help='TIC identifier of target.', nargs='?', default=2)
    parser.add_argument('-s', '--site', type=int, help='TIC identifier of target.', nargs='?', default=None)
    parser.add_argument('-d', '--date', type=str, help='TIC identifier of target.', nargs='?', default=None)
    parser.add_argument('-o', '--output', type=str, help='Base output directory for visibility plots.', nargs='?', default=None)
    args = parser.parse_args()

    if args.output is None:
        plots_interactive()

    t = Target.from_tname(args.target)

    output = args.output
    if output is not None:
        regex_year = re.compile(r'^(\d+)')
        m = regex_year.match(t.name)
        year = int(m.group(1))
        output = os.path.join(output, str(year), t.name)
        if not os.path.exists(output):
            os.makedirs(output)

    flows.visibility(target=t, siteid=args.site, date=args.date, output=output)
