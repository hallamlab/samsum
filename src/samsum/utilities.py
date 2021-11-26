import subprocess
import logging
import sys
import os
import re

import multiprocessing
import tqdm

LOGGER = logging.getLogger("samsum")


def file_prefix(file_path: str) -> str:
    return os.path.basename('.'.join(file_path.split('.')[:-1]))


def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)


def which(program):
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path_element in os.environ["PATH"].split(os.pathsep):
            path_element = path_element.strip('"')
            exe_file = os.path.join(path_element, program)
            if is_exe(exe_file):
                return exe_file
    return None


def launch_write_command(cmd_list, just_do_it=False, collect_all=True):
    """
    Wrapper function for opening subprocesses through subprocess.Popen()

    :param cmd_list: A list of strings forming a complete command call
    :param just_do_it: Always return even if the returncode isn't 0
    :param collect_all: A flag determining whether stdout and stderr are returned
    via stdout or just stderr is returned leaving stdout to be written to the screen
    :return: A string with stdout and/or stderr text and the returncode of the executable
    """
    stdout = ""
    if collect_all:
        proc = subprocess.Popen(cmd_list,
                                shell=False,
                                preexec_fn=os.setsid,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT)
        stdout = proc.communicate()[0].decode("utf-8")
    else:
        proc = subprocess.Popen(cmd_list,
                                shell=False,
                                preexec_fn=os.setsid)
        proc.wait()

    # Ensure the command completed successfully
    if proc.returncode != 0 and not just_do_it:
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
            stdout, returncode = launch_write_command([exe_dict[exe], "-v"], True)
        elif exe in no_params:
            stdout, returncode = launch_write_command([exe_dict[exe]], True)
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


def tqdm_multiprocessing(func, arguments_list: list, num_processes: int, pbar_desc: str, disable=False) -> list:
    if len(arguments_list) == 0:
        return []
    pool = multiprocessing.Pool(processes=num_processes)

    jobs = []
    result_list_tqdm = []
    pbar = tqdm.tqdm(jobs, total=len(arguments_list), desc=pbar_desc, ncols=120, disable=disable)

    def update(*a):
        pbar.update()

    for args in arguments_list:
        jobs.append(pool.apply_async(func=func, args=(*args,), callback=update))
    pool.close()

    for job in pbar:
        result_list_tqdm.append(job.get())

    pbar.close()

    return result_list_tqdm


def make_sure_dir_exists(dir_path: str) -> None:
    """
    Check whether a directory exists and create it if not not found.

    :param dir_path:
    :return:
    """

    if not dir_path or os.path.isdir(dir_path):
        return

    try:
        os.mkdir(dir_path)
    except OSError:
        LOGGER.error("Unable to create directory for path '{}'\n".format(dir_path))
        sys.exit(1)

    return


def int_to_str(number: int, n_digits: int) -> str:
    x = str(number)
    if len(x) > n_digits:
        raise SystemExit("Unable to convert integer {} to string with {} digits.".format(number, n_digits))
    while len(x) < n_digits:
        x = '0' + x
    return x
