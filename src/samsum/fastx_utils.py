import logging
import os
import sys

from pyfastx import Fasta, Fastq
from pyfastxcli import fastx_format_check

LOGGER = logging.getLogger("samsum")


def gen_fastx(fastx: str):
    if not os.path.exists(fastx):
        LOGGER.debug("'{}' fastx file doesn't exist.\n".format(fastx))
        return

    try:
        fastx_type = fastx_format_check(fastx)
    except Exception:
        LOGGER.error("Unable to detect file type for '{}'\n".format(fastx))
        sys.exit(3)

    try:
        if fastx_type == 'fasta':
            return Fasta(fastx, build_index=False, full_name=True)
        elif fastx_type == 'fastq':
            return Fastq(fastx, build_index=False, full_name=True)
        else:
            LOGGER.error("Unknown fastx type: '{}'\n".format(fastx_type))
            sys.exit(3)
    except RuntimeError as error:
        LOGGER.debug(str(error) + "\n")
        return
