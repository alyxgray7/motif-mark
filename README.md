# motif-mark
***An algorithm to visualize RNA-binding protein motifs on a transcript sequence.***

### The problem:
Alternative splicing events are important modes of genetic regulation and is a source of protein diversity from the genome [^1]. Typically, splice sites at the intron/exon junctions direct the excision of introns from pre-mRNA transcripts by spliceosomes, but some regulatory sequences (motifs) bind other RNA-binding proteins (RBPs) that allow for multiple protein isoforms to be encoded by the same gene. The exonic regions removed could be specific to different cell types or in response to cellular stresses / changes in environment that require a certain protein isoform to perform a desired function. After removing the intronic and desired exonic regions, the pre-mRNA transcript is processed to produce a mature mRNA and distinct protein isoform that is functionally different. Alternative splicing mechanisms are highly variable and splicing variants are implicated in many human diseases.

Cassette exon (or exon skipping) is one example of an alternative splicing mechanism where one exon is left out of the primary transcript. It is the most common splicing mechanism in mammalian tissue. 

![Image of exon skipping](../casset_exon.png) [^1]
[^1]: Image adapted from Wikipedia (https://en.wikipedia.org/wiki/Alternative_splicing)

***Recognizing motifs is important for understanding alternative splicing mechanisms.***

### How **motif-mark** works:
Utilizing Pycairo, a *to-scale transcript map* will be drawn for any given gene/genes to easily visualize the introns, exons, and motfis.
- Introns will be labeled with a thin, horizontal line.
- Exons will be drawn as a large, rectangular box.
- Motifs will be labeled with tick marks, marked at the beginning of the motif sequence. If the motif sequence is repeated consecutively in the sequence, a tick mark will be placed at the beginning of each open reading frame.
- Different motif sequences will be labeled with different colors. A legend will be printed on the top right of each outputted image.
- Each gene will be labeled with the header information contained in the `<file.fasta>`.

![Example output](/Users/agray11/bioinformatics/WINTER2021/BI625_ADVGEN/MOTIF_MARK/motif-mark/Figure_1.svg)

- For each `<file.fasta>`, one image will be drawn containing all genes contained 
in the file. *While any number of genes can be included in the file, the recommended max is 10.* 
- Image will be saved in `<file.svg>` format. 
- Gene map is drawn to scale.

***

### Requirements:
- A `<file.fasta>` containing gene sequences. 
    - The file can contain as many genes as wanted (the recommended max is 10).
    - Header information should be unique and contain information desired to identify the transcript. For example: `>INSR chr19:7150261-7150808 (reverse complement)`
    - Sequences should contain **lower-case** letters for *introns* and **UPPER-CASE** letters for *exons*.

- A `<file.txt>` containing motif sequences. 
    - All ambigous base changes are considered.
    - Any motif sequence in introns and exons will be identified, so upper or lower case letters don't need to be specified.
    - Each motif sequence should read line-by-line. Example `<file.txt>`:
        
        >ygcy
        >GCAUG
        >catag
        >YYYYYYYYYY
        