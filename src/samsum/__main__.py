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


def main(cmd_args=None) -> int:
    """
    Main control function of samsum.
    Based on the subcommand received from the command-line it will call the corresponding function in samsum.commands.

    :return: None
    """
    commands = {"stats": stats,
                "info": info}
    parser = argparse.ArgumentParser(description='Summarize read recruitments to reference sequences')
    parser.add_argument('command', nargs='?')
    if not cmd_args:
        cmd_args = sys.argv
    args = parser.parse_args(cmd_args[1:2])

    if not args.command:
        sys.stderr.write(usage)
        return 1

    if args.command not in commands:
        logging.error('Unrecognized command')
        sys.stderr.write(usage)
        return 1

    cmd = commands.get(args.command)
    cmd(cmd_args[2:])
    logging.info("samsum has finished successfully.\n")
    return 0


if __name__ == '__main__':
    sys.exit(main())
