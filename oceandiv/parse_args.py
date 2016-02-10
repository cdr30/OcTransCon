"""
Routines to parse command line arguments.

"""

import argparse


def get_args():
    """
    Get arguments from command line.
    
    """
    parser = argparse.ArgumentParser(
        description='Calculate ocean transport convergences.')
    parser.add_argument(
        'namelist', type=str, help='Path to namelist.ini')
    args = parser.parse_args()

    return args
