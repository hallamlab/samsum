import logging
import os
import sys


class SamsumFormatter(logging.Formatter):
    """
    Controls the formatting of the log through the logging package
    """
    error_fmt = "%(levelname)s - %(module)s, line %(lineno)d:\n%(message)s"
    warning_fmt = "%(levelname)s:\n%(message)s"
    debug_fmt = "%(asctime)s\n%(message)s"
    info_fmt = "%(message)s"

    def __init__(self):
        super().__init__(fmt="%(levelname)s: %(message)s",
                         datefmt="%d/%m %H:%M:%S")

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._style._fmt = SamsumFormatter.debug_fmt

        elif record.levelno == logging.INFO:
            self._style._fmt = SamsumFormatter.info_fmt

        elif record.levelno == logging.ERROR:
            self._style._fmt = SamsumFormatter.error_fmt

        elif record.levelno == logging.WARNING:
            self._style._fmt = SamsumFormatter.warning_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)

        # Restore the original format configured by the user
        self._style._fmt = format_orig

        return result


def prep_logging(log_file_name=None, verbosity=False, silent=False, stream=sys.stderr) -> logging.Logger:
    """
    Allows for multiple file handlers to be added to the root logger, but only a single stream handler.
    The new file handlers must be removed outside of this function explicitly

    :param log_file_name: Path to a file to write the TreeSAPP log
    :param verbosity: Whether debug-level information should be written (True) or not (False)
    :param stream: Which stream, sys.stdout or sys.stderr, should the console logger write to?
    :param silent: Suppresses writing to the output stream if True
    :return: A logging.Logger instance with the name "samsum"
    """
    if verbosity:
        logging_level = logging.DEBUG
    else:
        logging_level = logging.INFO

    # Detect whether a handlers are already present and return if true
    ss_logger = logging.getLogger("samsum")
    if silent:
        ss_logger.disabled = True
    ss_logger.setLevel(logging.DEBUG)
    if len(ss_logger.handlers):
        return ss_logger

    formatter = SamsumFormatter()
    # Set the console handler normally writing to stdout/stderr
    ch = logging.StreamHandler(stream=stream)
    ch.setLevel(logging_level)
    ch.terminator = ''
    ch.setFormatter(formatter)
    ss_logger.addHandler(ch)

    ss_logger.is_silent = False
    if silent:
        ss_logger.is_silent = True
        ch.setLevel(logging.ERROR)
        ss_logger.disabled = True

    if log_file_name:
        if not os.path.isabs(log_file_name):
            log_file_name = os.path.join(os.getcwd(), os.path.dirname(log_file_name), os.path.basename(log_file_name))
        log_dir = os.path.dirname(log_file_name)
        try:
            if log_dir and not os.path.isdir(log_dir):
                os.mkdir(log_dir)
        except (IOError, OSError):
            sys.stderr.write("ERROR: Unable to make directory '" + log_dir + "'.\n")
            sys.exit(3)
        ts_file_logger = logging.FileHandler(log_file_name, 'w')
        ts_file_logger.setFormatter(SamsumFormatter())
        ts_file_logger.setLevel(logging.DEBUG)
        ss_logger.addHandler(ts_file_logger)
        ss_logger.propagate = False
    return ss_logger
