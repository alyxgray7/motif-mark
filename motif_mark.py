#!/usr/bin/env python

### Import methods ### 
import argparse
import gzip
import re

### ARGPARSE ###
# def get_args():
#     """
#     Will set a command line option to run arguments.
#     """
#     # Parser main
#     parser = argparse.ArgumentParser(description = "Remove PCR duplicates from a SAM file. \
#         SAM file must be previously sorted by ascending chromosome using tools such as \
#         SAMtools. Sorted SAM files should be unzipped and contain reads from a single-end \
#         sequencing experiment.")
    
#     # Parser arguments
#     parser.add_argument('-d', '--directory', help = "Desired directory for output files to \
#         write in. ", required = True)
#     parser.add_argument('-f','--SAMfile', help = "Absolute/path/to/<SAMfile.sam>. The file \
#         should be presorted by chromosome location using a program such as SAMtools.", \
#         required = True)
#     parser.add_argument('-e','--which_end', help = "Argument to specify if <SAMfile.sam> \
#         contains single or paired-end reads. Use 1 for single and 2 for paired. \
#         WARNING: this script does no account for paired-end reads.", required = True, \
#         type = int, nargs = 1)
#     parser.add_argument('-u', '--UMIfile', help = "Absolute/path/to/<UMIfile.txt> \
#         containing a list of UMIs used during the sequencing run.", required = True)
#     #parser.add_argument('-h','--help', help = "Print help statement here.")
    
#     return parser.parse_args()
# args = get_args()

# vIUPAC = {
#     "A":["A"            ],
#     "C":[    "C"        ],
#     "G":[        "G"    ],
#     "T":[            "T"],
#     "U":[            "U"],
#     "W":["A",        "T"],
#     "S":[    "C","G"    ],
#     "M":["A","C"        ],
#     "K":[        "G","T"],
#     "R":["A",    "G",   ],
#     "Y":[    "C",    "T"],
#     "B":[    "C","G","T"],
#     "D":["A",    "G","T"],
#     "H":["A","C",    "T"],
#     "V":["A","C","G",   ],
#     "N":["A","C","G","T"],
#     "Z":[               ],
# }