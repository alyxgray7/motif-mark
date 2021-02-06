#!/usr/bin/env python

### Import methods ### 
import argparse
import gzip
import re
import cairo

### FILES ###
fasta_file="/Users/agray11/bioinformatics/WINTER2021/BI625_ADVGEN/Figure_1.fasta"
motif_file="/Users/agray11/bioinformatics/WINTER2021/BI625_ADVGEN/Fig_1_motifs.txt"

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

IUPAC_dict = {
    "A":["A"            ],
    "C":[    "C"        ],
    "G":[        "G"    ],
    "T":[            "T","U"],
    "U":[            "U","T"],
    "W":["A",        "T"],
    "S":[    "C","G"    ],
    "M":["A","C"        ],
    "K":[        "G","T"],
    "R":["A",    "G",   ],
    "Y":[    "C",    "T"],
    "B":[    "C","G","T"],
    "D":["A",    "G","T"],
    "H":["A","C",    "T"],
    "V":["A","C","G",   ],
    "N":["A","C","G","T"],
    "Z":[               ],
}
#print(IUPAC_dict)

def save_motifs(motif_file):
    """
    Will read in the <motif.txt>, save the motif sequences, and create the possible 
    combinations of those sequences as a motif_dict. Will account for ambigous
    base changes. 

    The motif_dict will hold the motif sequence fed in by the <motif.txt> as keys
    and have values equal to it's possible nucleotide sequence combinations. 
    Example:

    {"ygcy":["tgct", "tgcc", "cgct", "cgcc", "ugcu", "ugcc", "cgcu"]}

    Needs to take into account specific and ambigous motifs-- if one sequence that 
    is specific (tgct) and ambiguous (ygcy), the specific one will go FIRST, then
    the ambiguous one will follow. Ambiguous key's can't have values of the specific 
    motif sequence. 
    Example:

    {"tgct":["tgct", "ugcu"], "ygcy":["tgcc", "cgct", "cgcc", "ugcc", "cgcu"]}
    """
    motif_dict = {}
    with open(motif_file, 'r') as motif_fh:
        LN = 0
        for motif in motif_fh:
            LN += 1
            motif = motif.strip()
            #motif_dict.setdefault(motif,[motif.lower(), motif.upper()])
            motif_dict.setdefault(motif,[motif, motif])

    return motif_dict

motif_dict = save_motifs(motif_file)
print(motif_dict) # NEED TO COME BACK TO THIS!!!!!

# motif_dict = {}
# with open(motif_file, 'r') as motif_fh:
#     LN = 0
#     for line in motif_fh:
#         LN += 1
#         motif = line.strip()

#         # Set all motif strings as keys
#         motif_dict.setdefault(motif,)

#         # Search for ambiguous bases in motif seqs
#         ambig = re.findall(r"[^agcAGC]", motif)
#         print(ambig)
        


# motifs = ["ygcy","GCAUG","catag","YYYYYYYYYY"]
# for motif in motifs:
#     motif2 = str(motif)
#     ambig = re.split(r"[^agcAGC]", motif2)
#     print(ambig)
#     # base = motif.split()
#     # print(base)

# for motif in motifs:
#     amb_motifs = re.search(r"\w[^agcAGC]", motifs)
#     print(amb_motifs.string())
    # if motif != "A" or "C" or "G":
    #     print(motif)

#for motif in motifs:
    #print(motif)
# txt = "The rain in Spain"
# x = re.split("\s", txt)
# print(x)

            # for i in motif_dict:
            #     amb_motifs = re.search(r"\w[^agcAGC]", motif_dict.keys())
            #     print(amb_motifs)


#print(motif_dict)
        #print(motif)
        # = umi.strip()
        # dictionary.setdefault(umi,0)
    #return dictionary

# def multiple_replace(dict, text):
#   # Create a regular expression  from the dictionary keys
#   regex = re.compile("(%s)" % "|".join(map(re.escape, dict.keys())))

#   # For each match, look-up corresponding value in dictionary
#   return regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text) 


# text = "Larry Wall is the creator of Perl"

# dict = {
# "Larry Wall" : "Guido van Rossum",
# "creator" : "Benevolent Dictator for Life",
# "Perl" : "Python",
# } 

# print(multiple_replace(dict, text))

