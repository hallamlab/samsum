__author__ = 'Connor Morgan-Lang'


import subprocess
import logging
import sys
import os
import re
from samsum import classy


def launch_write_command(cmd_list, collect_all=True):
    """
    Wrapper function for opening subprocesses through subprocess.Popen()

    :param cmd_list: A list of strings forming a complete command call
    :param collect_all: A flag determining whether stdout and stderr are returned
    via stdout or just stderr is returned leaving stdout to be written to the screen
    :return: A string with stdout and/or stderr text and the returncode of the executable
    """
    stdout = ""
    if collect_all:
        proc = subprocess.Popen(' '.join(cmd_list),
                                shell=True,
                                preexec_fn=os.setsid,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        stdout = proc.communicate()[0].decode("utf-8")
    else:
        proc = subprocess.Popen(' '.join(cmd_list),
                                shell=True,
                                preexec_fn=os.setsid)
        proc.wait()

    # Ensure the command completed successfully
    if proc.returncode != 0:
        logging.error(cmd_list[0] + " did not complete successfully! Command used:\n" +
                      ' '.join(cmd_list) + "\nOutput:\n" + stdout)
        sys.exit(19)

    return stdout, proc.returncode


def executable_dependency_versions(exe_dict):
    """
    Function for retrieving the version numbers for each executable in exe_dict
    :param exe_dict: A dictionary mapping names of software to the path to their executable
    :return: A formatted string with the executable name and its respective version found
    """
    versions_dict = dict()
    versions_string = "Software versions used:\n"

    simple_v = ["prodigal"]
    no_params = ["bwa"]
    version_re = re.compile(r"[Vv]\d+.\d|version \d+.\d|\d\.\d\.\d")

    for exe in exe_dict:
        ##
        # Get the help/version statement for the software
        ##
        versions_dict[exe] = ""
        if exe in simple_v:
            stdout, returncode = launch_write_command([exe_dict[exe], "-v"])
        elif exe in no_params:
            stdout, returncode = launch_write_command([exe_dict[exe]])
        else:
            logging.warning("Unknown version command for " + exe + ".\n")
            continue
        ##
        # Identify the line with the version number (since often more than a single line is returned)
        ##
        for line in stdout.split("\n"):
            if version_re.search(line):
                # If a line was identified, try to get just the string with the version number
                for word in line.split(" "):
                    if re.search(r"\d\.\d", word):
                        versions_dict[exe] = re.sub(r"[,:()[\]]", '', word)
                        break
                break
            else:
                pass
        if not versions_dict[exe]:
            logging.debug("Unable to find version for " + exe + ".\n")

    ##
    # Format the string with the versions of all software
    ##
    for exe in sorted(versions_dict):
        n_spaces = 12-len(exe)
        versions_string += "\t" + exe + ' '*n_spaces + versions_dict[exe] + "\n"

    return versions_string


def load_references(refseq_lengths: dict) -> dict:
    """

    :param refseq_lengths:
    :return:
    """
    logging.debug("Loading the reference sequences into objects... ")
    references = {}
    for seq_name in refseq_lengths:  # type: str
        if seq_name in references:
            logging.error("Duplicate reference sequence names encountered: %s\n" % seq_name)
            sys.exit(3)
        ref_seq = classy.RefSequence(seq_name)
        references[seq_name] = ref_seq
        ref_seq.length = refseq_lengths[seq_name]
    logging.debug("done.\n")
    return references


def load_alignments(mapped_dict: dict) -> list:
    """

    :param mapped_dict:
    :return:
    """
    alignments = []
    for read_name in mapped_dict:
        map_stats = mapped_dict[read_name].split("\t")
        ref_seq = classy.AlignmentDat(read_name)
        ref_seq.load_sam(map_stats)
        alignments.append(ref_seq)
    return alignments


def load_reference_coverage(references: dict, alignments: list) -> None:
    """

    :param references:
    :param alignments:
    :return:
    """
    for aln in alignments:  # type: classy.AlignmentDat
        if aln.ref not in alignments:
            logging.error("Reference sequence from SAM file not found in FASTA: %s\n" % aln.ref)
            sys.exit(3)
        ref_seq = references[aln.ref]  # type: classy.RefSequence
        ref_seq.reads_mapped += 1
        ref_seq.weight_total += aln.weight
        if aln.start < ref_seq.leftmost:
            ref_seq.leftmost = aln.start
    return
