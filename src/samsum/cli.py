import os

import mappy
# import pysam
import pyfastx

from samsum import _version as ss_version
from samsum import args as ss_args
from samsum import classy as ss_class
from samsum import logger as ss_log
from samsum import file_parsers as ss_fp
from samsum import utilities as ss_utils


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
    ss_logger = ss_log.prep_logging()
    info_ss = ss_class.SAMSumBase("info")

    ss_logger.info("samsum version " + ss_version.__version__ + ".\n")

    # Write the version of all python deps
    py_deps = {"mappy": mappy.__version__,
               # "pysam": pysam.__version__,
               "pyfastx": pyfastx.version()}

    ss_logger.info("Python package dependency versions:\n\t" +
                   "\n\t".join([k + ": " + v for k, v in py_deps.items()]) + "\n")

    # TODO: Write path to installation directory

    # Write the version of executable deps
    ss_logger.info(ss_utils.executable_dependency_versions(info_ss.executables))

    if args.verbose:
        pass
        # ss_logger.info(summary_str)

    return 0


def stats(sys_args):
    """
    A user-facing sub-command to write an abundance table from provided SAM and FASTA files.

    :param sys_args: List of arguments parsed from the command-line.
    :return: None
    """
    parser = ss_args.SAMSumArgumentParser(description="Calculate read coverage stats over reference sequences.")
    parser.add_stats_args()
    args = parser.parse_args(sys_args)

    stats_ss = ss_class.SAMSumBase("stats")
    stats_ss.log = ss_log.prep_logging(log_file_name=os.path.dirname(args.output_table) + os.sep + "samsum_log.txt",
                                       verbosity=args.verbose)
    stats_ss.furnish_with_arguments(args)
    references = stats_ss.calculate_abundances()

    # Write the summary table with each of the above metrics as well as variance for each
    ss_fp.write_summary_table(references, args.output_table,
                              ss_utils.file_prefix(stats_ss.aln_file), stats_ss.unmapped, args.sep)

    return 0
