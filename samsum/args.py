__author__ = 'Connor Morgan-Lang'

import argparse


class SAMSumArgumentParser(argparse.ArgumentParser):
    """
    A base argparse ArgumentParser for samsum with functions to furnish with common arguments.
    This standardizes the interface for a unified aesthetic across all sub-commands
    """
    def __init__(self, **kwargs):
        """
        Instantiate the argparse argument-parser and create three broad argument groups:
            reqs - for the required parameters
            optopt - for the optional parameters
            miscellany - for the miscellaneous parameters that are module agnostic,
            for example verbose, help, num_threads
        :param kwargs:
        """
        super(SAMSumArgumentParser, self).__init__(add_help=False, **kwargs)
        self.reqs = self.add_argument_group("Required parameters")
        self.seqops = self.add_argument_group("Sequence operation arguments")
        self.optopt = self.add_argument_group("Optional options")
        self.miscellany = self.add_argument_group("Miscellaneous options")

        self.miscellany.add_argument("-v", "--verbose", action="store_true", default=False,
                                     help="Prints a more verbose runtime log")
        self.miscellany.add_argument("-h", "--help",
                                     action="help",
                                     help="Show this help message and exit")

    def parse_args(self, args=None, namespace=None):
        args = super(SAMSumArgumentParser, self).parse_args(args=args, namespace=namespace)

        return args
