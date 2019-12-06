__author__ = 'Connor Morgan-Lang'

"""
samsum command line.
"""
import sys
import argparse
import logging

from samsum.commands import (info, stats)

usage = """
samsum <command> [<args>]
** Commands include:
stats          Write the number of reads that mapped to each reference sequence
** Other commands:
info           Display samsum version and other information.
Use '-h' to get subcommand-specific help, e.g.
"""


def main():
    commands = {"stats": stats,
                "info": info}
    parser = argparse.ArgumentParser(description='Summarize read recruitments to reference sequences')
    parser.add_argument('command', nargs='?')
    args = parser.parse_args(sys.argv[1:2])

    if not args.command:
        sys.stderr.write(usage)
        sys.exit(1)

    if args.command not in commands:
        logging.error('Unrecognized command')
        sys.stderr.write(usage)
        sys.exit(1)

    cmd = commands.get(args.command)
    cmd(sys.argv[2:])
    logging.info("samsum has finished successfully.\n")


if __name__ == '__main__':
    main()
