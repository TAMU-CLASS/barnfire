#! /usr/bin/env python

'''
Andrew Till
Winter 2016

Initialize directories and files needed for the barnfire framework
'''

# STDLIB
import sys
import os

#MINE
import directories as d

def initialize_scratch_space(callingFile):
    d.make_scratch_directories()
    d.copy_xnjoy()
    d.copy_common_group_structures()
    d.copy_xs_script(callingFile)

if __name__ == '__main__':
    callingFile = None
    if len(sys.argv) > 2:
        callingFile = os.path.join(sys.argv[1], sys.argv[2])
    initialize_scratch_space(callingFile)
