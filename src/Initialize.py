#! /usr/bin/env python

'''
Andrew Till
Winter 2016

Initialize directories and files needed for the barnfire framework
'''

# STDLIB
import sys
import os
import shutil

#MINE
from directories import make_scratch_directories, copy_xnjoy, copy_common_group_structures
import directories as dirs

def initialize_scratch_space():
    make_scratch_directories()
    copy_xnjoy()
    copy_common_group_structures()
    
if __name__ == '__main__':
    initialize_scratch_space()
