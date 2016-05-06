"""
Useful functions
"""


import numpy as np


def print_message(config, msg):
    """ Print message to standard output """
    if config.getboolean('output', 'print_stdout'): 
        print '\n%s\n' % (msg) 


def print_progress(task_name, nmax, n, nbar=20):
    """ Print progress to standard out. """
    
    done = nbar * '|'
    todo = nbar * '.'
    flt_pct = 100. * np.float(n)/nmax
    progind = np.int(flt_pct)/(100/nbar)
    progbar = done[:progind] + todo[progind:]
    
    print ('\r%25s: %s %6.2f%%' %
          (task_name, progbar, flt_pct)),
    
    if np.int(flt_pct) == 100:
        print ''
