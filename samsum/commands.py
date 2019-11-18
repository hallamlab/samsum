__author__ = 'Connor Morgan-Lang'


import logging
import samsum
import numpy
import samsum.utilities as ss_utils
import samsum.args as ss_args
# import samsum.file_parsers as ss_fp
import samsum.classy as ss_class


def info(sys_args):
    """

    """
    parser = ss_args.samsum_ArgumentParser(description="Return package and executable information.")
    args = parser.parse_args(sys_args)
    ss_class.prep_logging()
    ts_info = ss_class.samsum_base("info")

    logging.info("TreeSAPP version " + samsum.version + ".\n")

    # Write the version of all python deps
    py_deps = {"numpy": numpy.__version__}

    logging.info("Python package dependency versions:\n\t" +
                 "\n\t".join([k + ": " + v for k, v in py_deps.items()]) + "\n")

    # Write the version of executable deps
    ts_info.furnish_with_arguments(args)
    logging.info(ss_utils.executable_dependency_versions(ts_info.executables))

    if args.verbose:
        pass
        # logging.info(summary_str)

    return


def stats(sys_args):

    return
