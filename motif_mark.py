#!/usr/bin/env python


### Import methods ### 
import argparse
import itertools
import re
import cairo
import seaborn as sns 


def get_args():
    """
    Will set a command line option to run arguments.
    """
    # Parser main
    parser = argparse.ArgumentParser(description = "Creates a gene map image of motif binding sites for \
        a given sequence or set of sequences. One SVG image is created for each input <file.fasta>. \
        The <file.fasta> name will be used to name the output <image.svg>. Image will be saved \
        in the directory given.")
    
    # Parser arguments
    parser.add_argument('-f','--FASTAfile', help = "Absolute/path/to/<file.fasta>.", required = True)
    parser.add_argument('-m','--MOTIFfile', help = "Absolute/path/to/<motif.txt>.", required = True)

    return parser.parse_args()
args = get_args()


# Get filenames and directory
fasta_file = args.FASTAfile
motif_file = args.MOTIFfile


def generate_motifs(motif_file):
    """
    Will read in the <motif.txt>, save the motif sequences, and create the possible 
    combinations of those sequences as a motif_dict. Will account for ambigous
    base changes. 

    The motif_dict will hold the motif sequence fed in by the <motif.txt> as keys
    and have values equal to it's possible nucleotide sequence combinations. 
    Example:

    {"ygcy":["tgct", "tgcc", "cgct", "cgcc", "ugcu", "ugcc", "cgcu"]}

    ### Not sure if this part is needed????
    # Needs to take into account specific and ambigous motifs-- if one sequence that 
    # is specific (tgct) and ambiguous (ygcy), the specific one will go FIRST, then
    # the ambiguous one will follow. Ambiguous key's can't have values of the specific 
    # motif sequence. 
    # Example:

    # {"tgct":["tgct", "ugcu"], "ygcy":["tgcc", "cgct", "cgcc", "ugcc", "cgcu"]}
    """
    motif_dict = {}
    motif_list = []

    # Read in motif.txt
    with open(motif_file, 'r') as motif_fh:
        LN = 0

        # Loop through file to get motif sequences
        for motif in motif_fh:
            LN += 1

            # Make motif sequences lowercase
            motif = motif.strip().lower()

            # Add motifs to list
            motif_list.append(motif)

            # IUPAC dictionary
            IUPAC_dict = {
                # Use lowercase letters
                "a":["a"            ],
                "c":[    "c"        ],
                "g":[        "g"    ],
                "t":[        "u","t"],
                "u":[        "u","t"],
                "w":["a",        "t"],
                "s":[    "c","g"    ],
                "m":["a","c"        ],
                "k":[        "g","t"],
                "r":["a",    "g",   ],
                "y":[    "c",    "t"],
                "b":[    "c","g","t"],
                "d":["a",    "g","t"],
                "h":["a","c",    "t"],
                "v":["a","c","g",   ],
                "n":["a","c","g","t"],
                "z":[               ],
            }

            # Create itertools object for each motif
            groups = itertools.groupby(motif, lambda char:char not in IUPAC_dict)
            splits = []
            
            # Loop through each itertools group object
            for b, group in groups:

                # Motif contains ambiguous character
                if b:
                    # Add translated base to end of sequence
                    splits.extend([[g] for g in group])
                
                # No ambiguous character in motif
                else:

                    # Add nucleotide base to end of sequence
                    for nuc in group:
                        splits.append(IUPAC_dict[nuc])
            
            # Create dictionary and add translated sequences as values
            combos = ["".join(p) for p in itertools.product(*splits)]
            
            # Make all motifs lowercase
            motif_dict[motif] = combos
                
    return motif_dict, motif_list
motif_dict, motif_list = generate_motifs(motif_file)

# Sanity check
# for i, item in enumerate(motif_dict):
#     print(i, len(motif_dict[item]))
#print(motif_dict["yyyyyyyyyy"])
#print(motif_list)


def store_fasta_info(fasta_file):
    """
    Takes in the original <file.fasta> and saves the headers and respective sequences
    in a dictionary.
    """
    fasta_dict = {}

    # Read in file.fasta
    with open(fasta_file, 'r') as fasta_fh:
        LN = 0
        for line in fasta_fh:
            LN += 1

            # Add headers as dictionary keys
            if line.startswith(">") == True:
                
                # Initialize counter for total sequence length
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
                    fasta_dict[header] = [long_seq, cnt]
                    #print(cnt)

    return fasta_dict
fasta_dict = store_fasta_info(fasta_file)
#print(fasta_dict,"\n\n")


def get_int_ex_coordinates(fasta_dict):
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
    int_ex_coord_dict = {}

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
            pos = region.start() + 1 # adjust to make 1-based
            intron_pos.append(pos)
        for region in exons:
            pos = region.start() + 1 # adjust to make 1-based
            exon_pos.append(pos)

        # Add pos_coordinates to dictionary
        int_ex_coord_dict[header] = [intron_pos, exon_pos]
            
    return int_ex_coord_dict
int_ex_coord_dict = get_int_ex_coordinates(fasta_dict)
#print(int_ex_coord_dict,"\n\n")


def get_motif_coordinates(motif_dict):
    """
    Takes in the motif_dictionary and outputs coordinates as a dictionary.
    Keys = gene number (0 based)
    Values = dictionary of original motif sequence from <motif.txt> (keys) and possible motif 
        sequence combination coordinates (values). Example:

    Input dictionary:
    {'ygcy': ['cgcc', 'cgct', 'tgcc', 'tgct']}

    Output dictionary:
    {1: {'ycgy': [6, 11, 28, 33]}, 2: ...
    }
    """

    gene_motif_coord_dict = {}
    # Loop for each gene
    for i, item in enumerate(fasta_dict):
        gene_num = (i + 1)
        sequence = fasta_dict[item][0].lower()

        # For this gene, loop through each motif key
        motif_coord_dict = {}
        for motif in motif_dict:
            combos = motif_dict[motif]

            # Loop through each combo 
            #motif_pos = []
            for combo in combos:

                # Create regex pattern based on unique motif combo
                pattern = str("(?:" + combo + ")")

                # Use regex to find patterns in sequence
                motifs = re.finditer(pattern, sequence)

                # Get starting positions of motifs
                motif_pos = []
                for pos in motifs:
                    pos = pos.start() + 1 # adjust to make 1-based
                    motif_pos.append(pos)
            
            # Add to inner dictionary 
            motif_coord_dict[motif] = motif_pos
        
        # Add to outer dictionary
        gene_motif_coord_dict[gene_num] = motif_coord_dict

    return gene_motif_coord_dict
motif_coord_dict = get_motif_coordinates(motif_dict)
#print(motif_coord_dict)


def name_image(file_name):
    """
    Names image.svg from path given to args. 
    Example:

    args.FASTAfile = "full/path/to/the/file_name.fasta"

    name_image(args.FASTAfile) = "file_name.svg"
    """
    # Make image name the same as filename
    image_name = "%s" % fasta_file
    image_name = image_name.split("/")
    image_name = image_name[-1]
    image_name = image_name.split(".")
    image_name = str(image_name[0]) + ".svg"
    return image_name
image_name = name_image(fasta_file)
#print(image_name)


def get_dims(fasta_dict):
    """
    Sets dimensions of the image.
    """

    lengths = []
    scale_x = 7
    scale_y = 100
    for seq_length in fasta_dict.values():
        lengths.append(seq_length[1])
        width = float(max(lengths) // scale_x)
        height = float((len(lengths) * scale_y) + 40) # Add extra space for legend at bottom
    return width, height
WIDTH, HEIGHT = get_dims(fasta_dict)
#print(WIDTH, HEIGHT)


# Alignment coordinates
x0,y0 = 10,10

# Scaling factors
scale_y = 100 # each gene map should fit within this height

# Generate color palette (Pete approved)
num_colors = len(motif_dict)
#print(num_colors)
palette = sns.color_palette("colorblind", num_colors)

# Set motifs to color
motif_colors = {}
for count, motif in enumerate(motif_list):
    motif_colors[motif] = palette[count]
#print(motif_colors)


### Begin drawing image
with cairo.SVGSurface(image_name, WIDTH, HEIGHT) as surface:
    # start object for SVG surface
    cr = cairo.Context(surface)

    # Header font information
    cr.set_source_rgb(0,0,0) # black text
    cr.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    cr.set_font_size(12)

    # Loop through and label headers
    cr.move_to(x0,y0)
    for i, header in enumerate(int_ex_coord_dict):
        cr.show_text(header)

        # set new y coordinate
        y1 = y0 + scale_y + (scale_y * i)
        #print(y1)
        cr.move_to(x0,y1)
    
    # Set starting intron position for first gene
    y1 = y0 + (scale_y // 2.5)
    #cr.move_to(x0, y1)

    # Loop through and draw intronic regions
    intron_y = 5 # Set height of intronic region
    for i, intron in enumerate(int_ex_coord_dict):
        intron_coordinates = int_ex_coord_dict[intron][0]
        intron_x = max(intron_coordinates)
        #print(intron_x)
        cr.rectangle(x0, y1, intron_x, intron_y)
        cr.fill()
        cr.set_source_rgb(0,0,0)
        
        # set new y coordinate
        y1 = y1 + scale_y + (scale_y * i)
        print(y1)
        cr.move_to(x0, y1)

    # Draw exonic regions
    exon_y = intron_y * 7 # Set height of exonic region
    #print(exon_y)
    for i, exon in enumerate(int_ex_coord_dict):
        exon_coordinates = int_ex_coord_dict[exon][1]
        exon0 = min(exon_coordinates) # exon start
        exon1 = max(exon_coordinates) # exon end
        exon_width = exon1 - exon0 # length of rectangle
        
        # Get rectangle start coordinates
        x1 = x0 + exon0 
        y2 = y0 + 25 + (100 * i)
        
        # Set starting exon position for first gene
        cr.rectangle(x1, y2, exon_width, exon_y)
        cr.fill()
        cr.set_source_rgb(0,0,0)

        # set new y coordinate
        y2 = y2 + 100
    
    # Label legend
    cr.move_to(x0 + 50, HEIGHT - 10)
    for i, motif in enumerate(motif_dict):
        cr.show_text(motif)

        # set new x coordinate
        x2 = x0 + 150 + (100 * i)
        cr.move_to(x2, HEIGHT - 10)
    
    # Draw colored boxes (Legend)
    cr.move_to(x0 + 20, HEIGHT - 10)
    x3 = x0 + 20
    for i, colors in enumerate(motif_colors):
        cr.rectangle(x3, HEIGHT - 30, 20, 20)
        cr.set_source_rgba(motif_colors[colors][0], motif_colors[colors][1], motif_colors[colors][2], .8)
        cr.fill()

        # set new x coordinate
        x3 = x3 + 100 + (2 * i)
        #print(x3)
        cr.move_to(x3, HEIGHT - 10)
    
    # Map motifs on gene
    motif_y = exon_y

    # Loop through for each gene
    for gene_num, items in enumerate(motif_coord_dict):
        coord_dict = motif_coord_dict[items]
        #print(coord_dict)

        # Get motifs and colors
        for motif_num, motifs in enumerate(coord_dict):
            colors = motif_colors[motifs]
            coordinates = coord_dict[motifs]
            #print(motif_num, coord_dict[motifs], colors)

            # Set the color for each motif
            cr.set_source_rgba(colors[0], colors[1], colors[2], 0.8)

            # Mark motifs
            for coordinate in coordinates:

                # Set new x,y start position
                x4 = x0 + coordinate
                y3 = y0 + 25 + (100 * gene_num)
                cr.rectangle(x4, y3, 2.5, motif_y)
                cr.fill()
                #print(gene_num, motif_num, coordinate, (x4,y3), colors)

    cr.show_page()
print("Image Saved")







    # # Split the motif list into sublists
    # # Loop through dictionary to get motif coordinates
    # for items in motif_coord_dict:
    #     coordinates = motif_coord_dict[items]
    #     print(coordinates)
    #     check = coordinates[0]
        
    #     # Loop through coordinates
    #     for coordinate in coordinates:

    #         if coordinate >= check:
    #             print(check, coordinate)
            
    #         else:
    #             pass
    #         #print(coordinate)
    #     #print(check)

        # Loop through coordinates
        # check = 0
        # for coordinate in coordinates:
        #     check = coordinate

        #     checklist.append(coordinate)
        # print(coordinates)

    # for items in motif_coord_dict:
    #     coordinates = motif_coord_dict[items]
    #     print(coordinates)
        
        # if motif_coord_dict[i] > motif_coord_dict[i -1]:
        #     print(i)
    # for motifs in motif_coord_dict:
    #     motif_coordinates = (motif_coord_dict[motifs])
    #     for coordinate in motif_coordinates:
    #         print(motif_coordinates[coordinate])
        #print(motif_coordinates)
        # cutoff = motif_coordinates[-1]
        # print(cutoff)
        # for coordinate in motif_coordinates:
        #     if motif_coordinates[coordinate] < 
        #     counter = coordinate
        #     print(counter)
            # gene_counter = coordinate
            # if coordinate 
            # print(gene_counter, coordinate)
    


# for intron in int_ex_coord_dict:
#     intron_coordinates = int_ex_coord_dict[intron][0]
#     print(intron_coordinates)
    

# # Loop through and draw intronic regions
# y2 = y0 + (90/2)
# cr.move_to(x0,y2)
# for introns in int_ex_coord_dict.values():

# draw_gene()


# for i, header in enumerate(int_ex_coord_dict):
#     print(i) # position
#     print(header) # header name

# for positions in int_ex_coord_dict:
#     print(int_ex_coord_dict[positions][1])
#     # introns = po
#     # exons = positions[1]
#     # print(exons)

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
#     for read in int_ex_coord_dict:
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

                # # combo letters
                # # lowercase letters
                # "a":["a","A"            ],
                # "c":[    "c","C"        ],
                # "g":[        "g","G"    ],
                # "t":[            "t","T","u","U"],
                # "u":[            "u","U","t","T"],
                # "w":["a","A","t",        "T"],
                # "s":[    "c","C","g","G"    ],
                # "m":["a","A","c","C"        ],
                # "k":[        "g","G","t","T"],
                # "r":["a","A",    "g","G"   ],
                # "y":[    "c","C"    "t","T"],
                # "b":[    "c","C","g","G","t","T"],
                # "d":["a","A"    "g","G","t","T"],
                # "h":["a","A","c","C",    "t","T"],
                # "v":["a","A","c","C","g","G"   ],
                # "n":["a","A","c","C","g","G","t","T"],
                # "z":[               ],
                # # combo letters
                # # uppercase letters
                # "A":["a","A"            ],
                # "C":[    "c","C"        ],
                # "G":[        "g","G"    ],
                # "T":[            "t","T","u","U"],
                # "U":[            "u","U","t","T"],
                # "W":["a","A","T",        "t"],
                # "S":[    "c","C","g","G"    ],
                # "M":["a","A","c","C"        ],
                # "K":[        "g","G","t","T"],
                # "R":["a","A",    "g","G"   ],
                # "Y":[    "c","C"    "t","T"],
                # "B":[    "c","C","g","G","t","T"],
                # "D":["a","A"    "g","G","t","T"],
                # "H":["a","A","c","C",    "t","T"],
                # "V":["a","A","c","C","g","G"   ],
                # "N":["a","A","c","C","g","G","t","T"],
                # "Z":[               ]