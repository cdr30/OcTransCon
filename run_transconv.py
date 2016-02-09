#!/usr/bin/env python2.7
"""
Script to execute TransConv from the command line.

"""

import sys

from transconv.transconv import main


if __name__ == '__main__':
    
    try:
        main()
    except KeyboardInterrupt as err:
        print err
        sys.exit()
        

