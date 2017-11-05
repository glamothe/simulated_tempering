# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 11:38:23 2016

@author: glamothe
"""

from ST_class import *
import sys
import os
from optparse import OptionParser

def check_file(name, comment=""):
    """ Checks if file exists in current directory. 
    """
    flag = os.path.exists(name)
    if not flag :
        msg =  "ERROR: file %s not found\n" %(name)
        if comment != "":
            msg += "ERROR: %s" %(comment)
        sys.exit(msg) 


###########################
##     MAIN ST PROGRAM    #
###########################

parser = OptionParser(usage="%prog --gro file.gro --mdp file.mdp --top file.top")
parser.add_option("--gro", action="store", type="string", dest="gro_filename", help="Gromacs structure file .gro")
parser.add_option("--mdp", action="store", type="string", dest="mdp_filename", help="Gromacs parameters file .mdp")
parser.add_option("--top", action="store", type="string", dest="top_filename", help="Gromacs topology file .top")
(options, args) = parser.parse_args()

if not options.gro_filename:
    parser.print_help()
    parser.error("option --gro is mandatory")
check_file(options.gro_filename)

if not options.mdp_filename:
    parser.print_help()
    parser.error("option --mdp is mandatory")
check_file(options.mdp_filename)

if not options.top_filename:
    parser.print_help()
    parser.error("option --top is mandatory")
check_file(options.top_filename)

# Create a simulated tempering object
st = ST(gro=options.gro_filename, mdp=options.mdp_filename, top=options.top_filename)
# Run the simulated tempering
st.start()