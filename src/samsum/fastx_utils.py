import logging
import os
import sys

from pyfastx import Fasta, Fastq
from pyfastxcli import fastx_format_check

from samsum import utilities

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


def split_fastx(fastx: str, output_dir: str, n_splits: int) -> list:
    """

    :param fastx:
    :return:
    """
    if n_splits == 1:
        return [fastx]

    utilities.make_sure_dir_exists(output_dir)
    split_file_prefix = output_dir + os.sep + os.path.basename(fastx)
    fx_size = 0
    buff_list = []

    # Load the FASTA records into memory
    fastx_iter = gen_fastx(fastx)
    for name, seq in fastx_iter:
        buff_list.append(">{}\n{}\n".format(name, seq))
        fx_size += len(buff_list[-1])

    max_split_size = fx_size / n_splits
    split_files = []
    file_size_acc = 0
    fastx_split_path = split_file_prefix + '.' + utilities.int_to_str(len(split_files),
                                                                      n_digits=len(str(n_splits)))
    out_f = open(fastx_split_path, 'w')
    while buff_list:
        record = buff_list.pop(0)
        out_f.write(record)
        file_size_acc += len(record)
        if file_size_acc >= max_split_size:
            out_f.close()
            split_files.append(fastx_split_path)
            fastx_split_path = split_file_prefix + '.' + utilities.int_to_str(len(split_files),
                                                                              n_digits=len(str(n_splits)))
            out_f = open(fastx_split_path, 'w')
            file_size_acc = 0

    out_f.close()
    split_files.append(fastx_split_path)

    assert len(split_files) == n_splits

    return split_files


def fasta_seq_lengths(fasta_file: str, min_seq_length=0) -> (dict, str):
    """
    Function for calculating the lengths of all sequences in a FASTA file.

    :param fasta_file: Path to a FASTA file to be parsed
    :param min_seq_length: The minimum length for a reference sequence to be included
    :return: A dictionary of sequence lengths indexed by their respective sequence names
    """
    seq_lengths_map = {}
    try:
        py_fa = Fasta(fasta_file, build_index=False, full_name=True)
    except RuntimeError as error:
        return seq_lengths_map, str(error)

    for name, seq in py_fa:  # type: (str, str)
        if len(seq) > min_seq_length:
            seq_lengths_map[name] = len(seq)

    return seq_lengths_map, ""
