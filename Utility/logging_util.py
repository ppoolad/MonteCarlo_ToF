''' 
A collection of utilities to help with logging
Feel free to add utilities as needed.

Dylan Vogel
2020-05-01

'''

# ============================================================================
# IMPORTS
# ============================================================================

import logging
import logging.config
import sys
import os

try:
    from ruamel.yaml import YAML
except ImportError:
    YAML = None

# ============================================================================
# GLOBALS
# ============================================================================

CRITICAL = logging.CRITICAL
ERROR = logging.ERROR
WARNING = logging.WARNING
INFO = logging.INFO
DEBUG = logging.DEBUG
NOTSET = logging.NOTSET

# ============================================================================
# FUNCTIONS
# ============================================================================

def setup_logging(name, debug_level=logging.INFO):
    ''' Sets up logging to stdout and returns a logger for use in your scripts

        Default logging level is INFO, but can be changed in the call
        For convenience I've re-defined the logging library levels in 
        this namespace, so the level can be supplied as `logging_util.DEBUG`
        The base library levels can be found here: 
        https://docs.python.org/3/library/logging.html#logging-levels

        Args:
            name (str):          name of the logger, easiest way is to pass '__main__'
            debug_level (int):   debug level
        
        Returns:
            logger(logging.logger):     returns a logging Logger with name and debug level set
        
        Example usage:
            logger = logsetup.SetupLogging(__name__, logging_util.DEBUG)
    '''

    logger = logging.getLogger(name)
    if logger.handlers:
        logger.debug("A logger has already been set up with the name {}".format(name))
    else: 
        consoleFormatter = logging.Formatter("%(levelname)-8s %(name)-8s: %(message)s")
        consoleHandler = logging.StreamHandler(sys.stdout)
        consoleHandler.setFormatter(consoleFormatter)
        logger.addHandler(consoleHandler)
        logger.setLevel(debug_level)
        logger.debug("Set up a logger for '{}'".format(name))
    
    return logger

def setup_logging_using_config(path):
    ''' Set a config file for logger configuration
    
        Args:
            config_file (str):  path to a YAML config file
    '''
    if YAML is None:
        logger.error('Could not import ruamel.yaml to read the configuration file. Please import it with `pip install ruamel.yaml`')

    if not os.path.isfile(path):
        raise FileNotFoundError('The path given does not point to a file: {}'.format(path))

    with open(path, 'r') as f:
        yaml = YAML(typ='safe')
        logging.config.dictConfig(yaml.load(f))