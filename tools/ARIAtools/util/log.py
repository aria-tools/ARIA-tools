# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#
# Author: Rohan Weeden
# Copyright 2020, by the California Institute of Technology. ALL RIGHTS
# RESERVED. United States Government Sponsorship acknowledged.
#
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
"""
Global logging configuration
"""
import os
import sys
import logging

# Inspired by
# https://stackoverflow.com/questions/384076/how-can-i-color-python-logging-output


class UnixColorFormatter(logging.Formatter):
    YELLOW = "\x1b[33;21m"
    RED = "\x1b[31;21m"
    BOLD_RED = "\x1b[31;1m"
    RESET = "\x1b[0m"

    COLORS = {
        logging.WARNING: YELLOW,
        logging.ERROR: RED,
        logging.CRITICAL: BOLD_RED
    }

    def __init__(self, fmt=None, datefmt=None, style="%", use_color=True):
        """
        Init Formatting
        """
        super().__init__(fmt, datefmt, style)
        # Save the old function so we can call it later
        self.__formatMessage = self.formatMessage
        if use_color:
            self.formatMessage = self.formatMessageColor

    def formatMessageColor(self, record):
        message = self.__formatMessage(record)
        color = self.COLORS.get(record.levelno)
        if color:
            message = "".join([color, message, self.RESET])
        return message


class CustomFormatter(UnixColorFormatter):
    """
    Adds levelname prefixes to the message on warning or above.
    """

    def formatMessage(self, record):
        message = super().formatMessage(record)
        if record.levelno >= logging.WARNING:
            message = ": ".join((record.levelname, message))
        return message


# logger = logging.getLogger("ARIAtools")
# logger.setLevel(logging.INFO)
# 
# stdout_handler = logging.StreamHandler(sys.stdout)
# stdout_handler.setFormatter(CustomFormatter(use_color=os.name != "nt"))
# 
# errorfile_handler = logging.FileHandler("error.log")
# errorfile_handler.setFormatter(logging.Formatter(
#     "[{asctime}] {funcName:>20}:{lineno:<5} {levelname:<10} {message}",
#     style="{"
# ))
# errorfile_handler.setLevel(logging.WARNING)
# 
# logger.addHandler(stdout_handler)
# logger.addHandler(errorfile_handler)
