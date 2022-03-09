"""
A collection of utilities to help with config file parsing.
Feel free to add utilities as needed.

Author: Dylan Vogel
Last Modified: 2021-02-18
"""

# ============================================================================
# IMPORTS
# ============================================================================

# python modules
import os

# third-party modules
import numpy as np
from ruamel.yaml import YAML
from ruamel.yaml.main import parse  # import as `pip install ruamel.yaml`

yaml = YAML()

# ============================================================================
# FUNCTIONS
# ============================================================================


# empty class to hold attributes
class ARGS(object):
    pass


def parse_yaml_config(path, required=[]):
    """Parses a provided yaml file and returns a python object, similar to
    how parse_args() returns. This should allow for interoperability.

    Args:
        path (str): path to the yaml config file
    Returns
    """
    # create arguments abject
    args = ARGS()

    if not os.path.isfile(path):
        raise FileNotFoundError(
            "The path given does not point to a file: {}".format(path)
        )

    with open(path, "r") as file:
        # load file
        params = yaml.load(file)
        # check if we got anything
        if params == None:
            raise ValueError(f"Specified file does not appear to be valid yaml")
        # read keys
        keys = params.keys()

        # check if user has provided required args
        if required != []:
            # check if all required args are present
            for arg in required:
                if arg not in keys:
                    raise NameError(
                        f'Could not find required parameter "{arg}" in file {path}'
                    )

        # populate the args object
        for arg in keys:
            # set the object attribute
            setattr(args, arg, params[arg])

    return args


def generate_tdc_pattern(xdim, ydim, xtdc, ytdc, path=None):
    """Generate a tdc pattern for inputting into a config file
    Args:
        xdim (int): number of sipms in the x dimension
        ydim (int): number of sipms in the y dimension
        xtdc (int): size of tdc group in the x dimension
        ytdc (int): size of tdc group in the y dimension
        path (str): (optional) path of the config file. This function will
            lazily append to the end of the file (i.e. you might need to then
            go in and fix it up)

    Returns:
        None: Will just print to terminal or append to file
    """
    if path is not None and not os.path.isfile(path):
        raise FileNotFoundError(
            "The path given does not point to a file: {}".format(path)
        )

    if xdim % xtdc:
        raise ValueError(
            f"The number of SIPMs in the x-dimension ({xdim}) does not divide cleanly into the number of TDCs ({xtdc})"
        )

    if xdim % xtdc:
        raise ValueError(
            f"The number of SIPMs in the y-dimension ({ydim}) does not divide cleanly into the number of TDCs ({ytdc})"
        )

    tdc_pattern = np.zeros((ydim, xdim))

    tdc_ind = 0
    for i in range(xdim // xtdc):
        for j in range(ydim // ytdc):
            tdc_pattern[
                (j * ytdc) : (j * ytdc + ytdc), (i * xtdc) : (i * xtdc + xtdc)
            ] = tdc_ind
            # increment to the next TDC
            tdc_ind += 1
    
    # convert to standard python list for easy formatting
    tdc_pattern = [list(row.astype('int')) for row in tdc_pattern]

    if path is None:
        print("Paste this into your config:\n")
        for row in tdc_pattern:
            print(f" - {row}")
    else:
        with open(path, 'a+') as file:
            file.write('\ntdc_pattern:\n')
            for row in tdc_pattern:
                file.write(f' - {row}\n')


if __name__ == "__main__":
    # run a test of the module
    with open('test.yaml', 'w+') as file:
        file.write('---\n')
        file.write('foo: bar\n')
        file.write('bar: foo\n')

    
    # should open with no errors
    args = parse_yaml_config('test.yaml', required=[])
    assert(args.foo == 'bar')
    assert(args.bar == 'foo')

    # should raise filenotfounderror
    try:
        args = parse_yaml_config('test.lmay', required=[])
        raise AssertionError("False filename did not raise FileNotfoundError")
    except FileNotFoundError:
        pass

    # should raise nameerror
    try:
        args = parse_yaml_config('test.yaml', required=['pong'])
        raise AssertionError("False filename did not raise NameError")
    except NameError:
        pass

    # generate a tdc pattern and write it to the file
    generate_tdc_pattern(2, 2, 1, 1, 'test.yaml')

    # parse the config and check that it's correct
    args = parse_yaml_config('test.yaml', required=[])
    assert(args.tdc_pattern == [[0, 2], [1, 3]])


    # clean up our mess
    os.remove('test.yaml')

    print("Functional test of config_util.py passed with no errors")
