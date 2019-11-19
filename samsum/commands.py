__author__ = 'Connor Morgan-Lang'


import logging
import samsum
import numpy
import samsum.args as ss_args
import samsum.classy as ss_class

# import samsum.file_parsers as ss_fp
# import samsum.utilities as ss_utils


def info(sys_args):
    """
    Function for writing version information about samsum and python dependencies.
    Other related info (citation, executable versions, etc.) should also be written through this sub-command.
    Create a SAMSumBase object for the `info` sub-command

    :param sys_args: List of arguments parsed from the command-line.
    :return: None
    """
    parser = ss_args.SAMSumArgumentParser(description="Return package and executable information.")
    args = parser.parse_args(sys_args)
    ss_class.prep_logging()
    info_ss = ss_class.SAMSumBase("info")

    logging.info("samsum version " + samsum.version + ".\n")

    # Write the version of all python deps
    py_deps = {"numpy": numpy.__version__}

    logging.info("Python package dependency versions:\n\t" +
                 "\n\t".join([k + ": " + v for k, v in py_deps.items()]) + "\n")

    # TODO: Write path to installation directory

    # Write the version of executable deps
    # info_ss.furnish_with_arguments(args)
    # logging.info(ss_utils.executable_dependency_versions(info_ss.executables))

    if args.verbose:
        pass
        # logging.info(summary_str)

    return


def stats(sys_args):
    """

    :param sys_args: List of arguments parsed from the command-line.
    :return: None
    """
    parser = ss_args.SAMSumArgumentParser(description="Calculate read coverage stats over reference sequences.")
    args = parser.parse_args(sys_args)
    ss_class.prep_logging()
    info_ss = ss_class.SAMSumBase("info")

    return
