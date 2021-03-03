import os, sys, re
import argparse, logging

from astropy.io import fits

from flows.fitscmd import COMMANDS

def _add():

    parser = argparse.ArgumentParser(description='Add FITS specific command.')
    parser.add_argument('command', choices=list(COMMANDS.keys()), help='command')
    parser.add_argument('parameters', help='command parameters')
    parser.add_argument('fitsfiles', nargs='+', help='FITS files')
    args = parser.parse_args()

    try:
        parameters = COMMANDS[args.command][1](*re.match(COMMANDS[args.command][0], args.parameters).groups())
    except:
        logging.critical('can not parse %s', args.parameters)
        quit(1)

    for f in args.fitsfiles:
        try:
            with fits.open(f, mode='update') as hdul:
                if not COMMANDS[args.command][2](*parameters, hdul):
                    logging.error('%s %s is not valid for %s', args.command, args.parameters, f)
                    continue
                hdul[0].header[''] = 'FLOWSCMD: %s = %s' % (args.command, args.parameters)
        except Exception as e:
            logging.error('could not open %s', f)
            continue

def _get():

    parser = argparse.ArgumentParser(description='Get FITS specific commands.')
    parser.add_argument('fitsfile', help='FITS file')
    args = parser.parse_args()

    with fits.open(args.fitsfile) as hdul:

        if not '' in hdul[0].header:
            return

        commands = [(i, c[10:]) for i, c in enumerate(hdul[0].header['']) if c.startswith('FLOWSCMD: ')]
        for i, c in commands:
            print('%%%dd) %%s' % (len(str(len(commands))),) % (i, c))

def _del():

    parser = argparse.ArgumentParser(description='Delete FITS specific command.')
    parser.add_argument('command', type=int, help='Command ID')
    parser.add_argument('fitsfile', help='FITS file')
    args = parser.parse_args()

    with fits.open(args.fitsfile, mode='update') as hdul:

        if not '' in hdul[0].header:
            logging.critical('no commands in fits')
            quit(1)

        commands = [i for i, c in enumerate(hdul[0].header['']) if c.startswith('FLOWSCMD: ')]
        if not args.command in commands:
            logging.critical('command id does not exist')
            quit(1)

        head = list(hdul[0].header[''])
        head.pop(args.command)
        del hdul[0].header['']

        for l in head:
            hdul[0].header[''] = l

if __name__ == '__main__':

    if len(sys.argv) > 1:

        method = sys.argv.pop(1)

        if not {
                'add' : _add,
                'get' : _get,
                'del' : _del,
        }.get(method, lambda *a: 1)():

            quit(0)

        sys.argv.insert(1, method)

    parser = argparse.ArgumentParser(description='Control FITS specific commands.')
    parser.add_argument('method', choices=('add', 'get', 'del'), help='method')
    args = parser.parse_args()
