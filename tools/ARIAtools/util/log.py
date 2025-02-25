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
import logging

FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'


# Set log level to warning for some third party packages
for logger in ['botocore', 'urllib3', 'rasterio', 'asyncio']:
    logging.getLogger(logger).setLevel(logging.WARNING)

logging.getLogger("asf_search").setLevel("ERROR")

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
