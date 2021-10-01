import logging
import os
import sys

from pyfastx import Fasta, Fastq
from pyfastxcli import fastx_format_check

LOGGER = logging.getLogger("samsum")


def read_fastq_to_dict(fastq_file: str, num_records=0) -> dict:
    fasta_dict = dict()
    num_parsed = 0

    try:
        py_fa = Fastq(fastq_file, build_index=False, full_name=True)
    except RuntimeError as error:
        LOGGER.debug(str(error) + "\n")
        return fasta_dict

    for name, seq, _qual in py_fa:  # type: (str, str)
        fasta_dict[name] = seq.upper()
        num_parsed += 1
        if 0 < num_records <= num_parsed:
            break

    return fasta_dict


def read_fasta_to_dict(fasta_file: str, num_records=0) -> dict:
    """
    Reads any fasta file using the pyfastx library

    :param fasta_file: Path to a FASTA file to be read into a dict
    :param num_records: The number of sequence records to parse and return (default is all)
    :return: Dict where headers/record names are keys and sequences are the values
    """
    fasta_dict = dict()
    num_parsed = 0

    if not os.path.exists(fasta_file):
        LOGGER.debug("'{}' fasta file doesn't exist.\n".format(fasta_file))
        return fasta_dict

    try:
        py_fa = Fasta(fasta_file, build_index=False, full_name=True)
    except RuntimeError as error:
        LOGGER.debug(str(error) + "\n")
        return fasta_dict

    for name, seq in py_fa:  # type: (str, str)
        fasta_dict[name] = seq.upper()
        num_parsed += 1
        if 0 < num_records <= num_parsed:
            break

    return fasta_dict


def read_fastx_to_dict(fastx: str, num_records=0) -> dict:
    if not os.path.exists(fastx):
        LOGGER.debug("'{}' fasta file doesn't exist.\n".format(fastx))
        return {}

    try:
        fastx_type = fastx_format_check(fastx)
    except Exception:
        LOGGER.error("Unable to detect file type for '{}'\n".format(fastx))
        sys.exit(3)

    if fastx_type == 'fasta':
        return read_fasta_to_dict(fastx, num_records)
    elif fastx_type == 'fastq':
        return read_fastq_to_dict(fastx, num_records)
    else:
        LOGGER.error("Unknown fastx type: '{}'\n".format(fastx_type))
        sys.exit(3)
