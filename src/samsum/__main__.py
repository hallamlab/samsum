__author__ = 'Connor Morgan-Lang'

"""
samsum command line.
"""
import sys
import argparse
import logging
import resource

from samsum.commands import (info, stats)

usage = """
samsum <command> [<args>]
** Commands include:
stats          Write the number of reads that mapped to each reference sequence
** Other commands:
info           Display samsum version and other information.
Use '-h' to get subcommand-specific help, e.g.
"""

def memory_limit():
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (get_memory() * 1024 / 2, hard))

def get_memory():
    with open('/proc/meminfo', 'r') as mem:
        free_memory = 0
        for i in mem:
            sline = i.split()
            if str(sline[0]) in ('MemFree:', 'Buffers:', 'Cached:'):
                free_memory += int(sline[1])
    return free_memory

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
    memory_limit() # Limitates maximun memory usage to 25% of RAM
    try:
        main()
    except MemoryError:
        sys.stderr.write('\n\nERROR: Memory Exception\n')
        sys.exit(1)
