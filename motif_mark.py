#!/usr/bin/env python

### Import methods ### 
import argparse
import gzip
import re
import cairo
#import numpy as np

### FILES ###
#fasta_file="/Users/agray11/bioinformatics/WINTER2021/BI625_ADVGEN/test2.fasta"
fasta_file="/Users/agray11/bioinformatics/WINTER2021/BI625_ADVGEN/test.fasta"
#fasta_file="/Users/agray11/bioinformatics/WINTER2021/BI625_ADVGEN/Figure_1.fasta"
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
            motif_dict.setdefault(motif,)
            m = re.search(r"[^atgcATGC]", motif)
            if m:
                ambig =  m.group()
                #print(ambig)
                
    return motif_dict
motif_dict = save_motifs(motif_file)
#print(motif_dict) # NEED TO COME BACK TO THIS!!!!!

def store_info(fasta_file):
    """
    Takes in the original <file.fasta> and saves the headers and respective sequences
    in a dictionary.
    """
    fasta_dict = {}
    with open(fasta_file, 'r') as fasta_fh:
        LN = 0
        for line in fasta_fh:
            LN += 1
            # Add headers as dictionary keys
            if line.startswith(">") == True:
                cnt = 0
                seq_list = []
                header = line.strip()
                
            # Add respective sequences as values
            else:
                seq = line.strip()
                seq_list.append(seq)

                # join the seqs to make one long sequence
                long_seq = "".join(seq_list)

                # Sum the number of nucleotides in sequence for each key
                for seq in seq_list:
                    cnt += len(seq)
                    # Add all information to dictionary
                    fasta_dict[header] = [long_seq, cnt] # needed to subtract extra counts

    return fasta_dict
fasta_dict = store_info(fasta_file)
#print(fasta_dict,"\n\n")

def get_coordinates(dictionary):
    """
    Takes in fasta_dict and outputs coordinates as a dictionary for each 
    header. Specifically, takes the sequence portion and turns into position
    coordinates for intron and exon positions as a dictionary. Intron coordinates
    are first and exon coordinates are second. Example:
    
    Input dictionary:
    {'>MBNL chr3:152446461-152447003': ['ataatATGCGGGGAaaaatcgggccccATTTGGGGCaaaa', 40]}

    Output dictionary:
    {'>MBNL chr3:152446461-152447003': [
        [1, 2, 3, 4, 5, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 37, 38, 39, 40], 
        [6, 7, 8, 9, 10, 11, 12, 13, 14, 28, 29, 30, 31, 32, 33, 34, 35, 36]
        ]
    }
    """
    coordinates_dict = {}

    # Loop through each header contained in dict
    for header in fasta_dict:

        # Get sequence line from header value
        sequence = fasta_dict[header][0]

        # Regex to find introns / exons
        exons = re.finditer(r"[A-Z]", sequence)
        introns = re.finditer(r"[a-z]", sequence)

        # Loop through and save positions to lists
        intron_pos = []
        exon_pos = []
        for region in introns:
            pos = region.end()
            intron_pos.append(pos)
        for region in exons:
            pos = region.end()
            exon_pos.append(pos)

        # Add pos_coordinates to dictionary
        coordinates_dict[header] = [intron_pos, exon_pos]
            
    return coordinates_dict
coordinates_dict = get_coordinates(fasta_dict)
#print(coordinates_dict,"\n\n")


# def label_gene(ctx, pos, label):
#     """
#     Adds gene name to drawing. Gene name is taken from header information 
#     in <file.fasta> which is stored in coordinates_dict.
#     """

#     # Height of space where gene name and drawing will live in
#     gene_height = 90
#     # Coordinates of label
#     x0, y0 = pos

#     for read in coordinates_dict:
#         ctx.set_source_rgb(0,0,0) # black font
#         ctx.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
#         ctx.set_font_size(12)
#         ctx.move_to(X,Y)
#         ctx.show_text(read)


# def draw_gene():
#     """
#     Takes dictionary of coordinates and gene information and draws 
#     """

#     # Extract lengths for image dimensions 
#     lengths = []
#     for seq_length in fasta_dict.values():
#         lengths.append(seq_length[1])
#     #print(lengths)

#     # Set image dimensions based on fasta info
#     IMG_WIDTH = (max(lengths) / 4.5)
#     IMG_HEIGHT = (len(lengths) * 90)
#     x0,y0 = 10,10

#     # Initialize cairo surface
#     ps = cairo.SVGSurface("svgfile.svg", IMG_WIDTH, IMG_HEIGHT)
#     cr = cairo.Context(ps)

#     # Header font information
#     cr.set_source_rgb(0,0,0) # black text
#     cr.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
#     cr.set_font_size(12)

#     # Loop through and label headers
#     cr.move_to(x0,y0)
#     for i, header in enumerate(coordinates_dict):
#         cr.show_text(header)
        
#         # set new y coordinate
#         y1 = y0 + 90 + (90 * i)
#         cr.move_to(x0,y1)
    
#     # Loop through and draw intronic regions
#     y2 = y0 + (90/2)
#     cr.move_to(x0,y2)
#     for introns in coordinates_dict.values():

# draw_gene()
# if __name__ == "__main__":
#     main()

# for i, header in enumerate(coordinates_dict):
#     print(i) # position
#     print(header) # header name

for positions in coordinates_dict:
    print(coordinates_dict[positions][1])
    # introns = po
    # exons = positions[1]
    # print(exons)

# set width and height of image based on sequence lengths and #
# SCALE = 9
# WIDTH = (max(lengths) * 10)
# HEIGHT = (len(lengths) * 10) * SCALE
# X,Y = 10,10
# SPACE = 30
# INTRON_HEIGHT = 3.5
# EXON_HEIGHT = 9
# print(WIDTH, HEIGHT)

# for i, length in enumerate(lengths):
#     print(i)
#     print(length)

# # Create SVG surface
# with cairo.SVGSurface("test.svg", WIDTH, HEIGHT) as surface:
#     # start object for SVG surface
#     context = cairo.Context(surface)

#     # Loop through each read
#     for read in coordinates_dict:
#         context.set_source_rgb(0,0,0) # black
#         context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
#         context.set_font_size(12)
#         context.move_to(X,Y)
#         context.show_text(read)

#         # Draw line for length of gene (intron)
#         for i, length in enumerate(lengths):
#             context.rectangle(X, Y + (SPACE + (SPACE * i)), (length * SCALE), INTRON_HEIGHT)
#         #context.rectangle(X,Y+30,(max(lengths) * 9),3.5)
#             context.set_source_rgb(0.35,0,0)
#             context.fill()

#         # Draw exons

#         #context.rectangle()
#         # context.rectangle(50, 50, 50, 50)
#         # set color of context
#         # Fill color inside

#         context.stroke()
# print("Image Saved")



# lengths = {}
# for i, seq_length in enumerate(fasta_dict.values()):
#     lengths[i] = seq_length[1]

# sorted_lengths = {
#     sorted(lengths.items(),
#     key = lambda item: item[1],
#     reverse = True)
# }
# print(sorted_lengths)
    # lengths.append(i, length[1])

    # print(lengths)
# def create_img(dictionary):
#     """
    
#     """

# def append_value(dict_obj, key, value):
#     """
#     Appends a key value to an existing key in the dictionary.
#     """
#     # Check if key exist in dict or not
#     if key in dict_obj:
#         # Key exist in dict.
#         # Check if type of value of key is list or not
#         if not isinstance(dict_obj[key], list):
#             # If type is not list then make it list
#             dict_obj[key] = [dict_obj[key]]
#         # Append the value in list
#         dict_obj[key].append(value)
#     else:
#         # As key is not in dict,
#         # so, add key-value pair
#         dict_obj[key] = value
#     return 

# count the number of nucleotides for each read and append to end of dictionary value
#for i in header_dict.values():
    #print(i)

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


        # # Variables to save positions
        # intron_pos = []
        # exon_pos = []
        # for pos in introns:
        #     pos_end = pos.end()
        #     intron_pos.append(pos_end)
        # for pos in exons:


        #     #print(pos_end)

        # # Initialize lists to hold bases
        # intron_seq = []
        # exon_seq = []

        # # Loop through bases for each sequence
        # for base in sequence:
        #     base_counter += 1
        #     chopped_seqs = []

        #     # Add bases

        #     # Reached end of sequence, add chopped sequences to dictionary
        #     if base_counter > sequence_length:
        #         chopped_dict[i] = chopped_seqs
            
        #     # Continue to add bases 
        #     else:
        #         # Base is in intronic region
        #         if base.islower() == True:
        #             intron_seq.append(base)

        #         # Base is in exonic region
        #         else:
        #             full_intron = "".join(intron_seq)
        #             chopped_seqs.append(full_intron)
        #             exon_seq.append(base)
                

#print(chopped_dict)

# Get sequences and sequence length from dictionary
#chopped_seqs = set()
# chopped_dict = {}
# chopped_seqs = []
# for i in fasta_dict:
#     #print(i)
#     chopped_seqs = []
#     chopped_dict[i] = chopped_seqs
#     for sequences in fasta_dict.values():
#         sequence = sequences[0]
#         sequence_length = sequences[1]
#         base_counter = 0
#         intron_seq = []
#         exon_seq = []
#         # loop letter by letter in sequence
#         if base_counter <= sequence_length
#             for base in sequence:
#                 base_counter += 1
#                 #print(base)
            
#                 # base is in intronic region
#                 if base.islower() == True:
#                     # join the exon sequence, add to chopped seqs list, then clear
#                     full_exon = "".join(exon_seq)
#                     chopped_seqs.append(full_exon)
#                     #chopped_seqs.add(full_exon)
#                     exon_seq.clear()

#                     # begin to add intronic sequence 
#                     intron_seq.append(base)
#                     #print(full_intron)

#                 # base is in exonic region
#                 else:
#                     # join the intron sequence, add to chopped seqs, then clear
#                     full_intron = "".join(intron_seq)
#                     chopped_seqs.append(full_intron)
#                     #chopped_seqs.add(full_intron)
#                     intron_seq.clear()

#                     # begin to add exonic sequence
#                     exon_seq.append(base)
#                     #print(exon_seq)
#             #chopped_seqs = list(dict.fromkeys(chopped_seqs))
#         else:

# print(chopped_seqs)


    # for header in fasta_dict:
    #     chopped_dict[header] = chopped_seqs
    #     #chopped_seqs.clear()
#print(chopped_dict)
    #print(sequence)
    #seq_split1 = re.sub(r"[A-Z]*[a-z]*", sequence)
    #seq_split1 = re.sub(r"([A-Z])", r" \1", sequence).split()
    #seq_split1 = re.split(r"[A-Z]*[a-z]*", sequence)
    #print(seq_split1)

    # introns = re.compile(r"[a-z]+")
    # exons = re.compile(r"[A-Z]+")
    # chopped = introns.findall(sequence) + exons.findall(sequence)
    # print(chopped)
    
    # exons = re.split(r"[^A-Z]*[a-z]", sequence)
    # chopped_seqs = [None]*(len(introns) + len(exons))
    # chopped_seqs[::2] = introns
    # chopped_seqs[1::2] = exons
    # #chopped_seqs = list(zip(introns, exons))
    # print(chopped_seqs)
    # chopped = re.split(r"[^a-z][A-Z]*", sequence)
    # print(chopped)
    # chopped = [s for s in re.split(r"([A-Z][^A-Z]*)", sequence) if s]
    # print(chopped)
    # chopped.append(intron)
    # exon = re.findall(r"[^a-z][A-Z]*", sequence)
    # chopped.append(exon)
    # print(chopped)
    # sequence_length = sequences[1]

    # # Split sequence into intron and exon chunks
    # intron_seq = []
    # exon_seq = []
    # # Check base by base and increment counter
    # for base in sequence:
    #     base_counter += 1
        
    #     #if base_counter <= sequence_length:
    #     # introns have lowercase letters
    #     if base.islower() == True:
    #         intron_seq.append(base)
    #         whole_intron = "".join(intron_seq)
    #         #print(whole_intron)
    #     else: 
    #         exon_seq.append(base)
    #         whole_exon = "".join(exon_seq)
    #         #print(whole_exon)
        
        # else:
        #     whole_intron = "".join(intron_seq)
        #     chopped_sequence.append(whole_intron)

        #     whole_intron = "".join(intron_seq)
        #     chopped_sequence.append(whole_intron)
        #     # exons have UPPERcase letters
        #     else:
        #         exon_seq.append(base)
        #         whole_exon = "".join(exon_seq)
        #         chopped_sequence.append(whole_exon)

        # # Reset base counter and move to next sequence
        # else:
            
        #     print(base_counter)
        #     base_counter = 0
    
# print(chopped_sequence)
# print(len(chopped_sequence))

                # whole_intron = "".join(intron_seq)
                # chopped_sequence.append(whole_intron)
                #print(intron_seq)
            # #while base_counter <= sequence_length:
            #     # Build intron sequence
            #     if base.islower() == True:
            #         intron_seq.append(base)
            #     else:
            #         whole_intron = "".join(intron_seq)
            #         chopped_sequence.append(whole_intron)
            #         exon_seq.append(base)
            # else:
            #     whole_intron = "".join(intron_seq)
            #     whole_exon = "".join(exon_seq)

            # while base.islower() == True:
            #     intron_seq.append(base)
            # else:
            #     whole_intron = "".join(intron_seq)
            #     chopped_sequence.append(whole_intron)

            # while base.isupper() == True:
            #     exon_seq.append(base)
            # else:
            #     whole_exon = "".join(exon_seq)
            #     chopped_sequence.append(whole_exon)
            
            # print(base_counter)
            # break

            
        #     whole_intron = "".join(intron_seq)
        #     chopped_sequence
        #     whole_exon = "".join(exon_seq)
        #     else:
        #         whole_intron = "".join(intron_seq)
        #         chopped_sequence.append(whole_intron)
        #         exon_seq.append(base)
        #     else:
        #         exons.append(base)
        #         whole_exon = "".join(exons)
        #         exon_seqs.append(whole_exon)
        #     print(intron_seqs)
        #     print(exon_seqs)
        #         #print(whole_exon)
        # #         seq_split_dict[header] = [whole_intron]
        # #     else:
        # #         exons.append(base)
        # #         whole_exon = "".join(exons)
        # #         seq_split_dict[header] = []
                