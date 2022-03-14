
"""
Counts the number of DOG-box motifs in a genes promotor region. The promotor regions gets selected
by moving the start of the gene (pulled from a GFF file) 500 bases upstream. Using samtools the genomic sequence of 500pb+TSS: TSS is selected.

"""


import io
import subprocess

from sys import argv

__author__ = "Joris van Steenbrugge"
__license__= "MIT"
__version__= "v1.0"


GFF_f = argv[1]
fasta = argv[2]
transcript_ids = argv[3:]

def shift_gene_start(start_pos, shift_width = 500):
    new_start = int(start_pos) - shift_width
    if new_start < 0:
        return False

    return str(new_start)

def shift_gene_end(end_pos, shift_width = 500):
    return str(int(end_pos) + shift_width)


def run_samtools(pos, fasta):
    pos_fmt = f"{pos[0]}:{pos[1]}-{pos[2]}" 
    strand = pos[3]
    
    command = ['samtools', 'faidx', fasta, pos_fmt]

    if strand == '-':
        command.append('-i')

    proc = subprocess.Popen(command, stdout=subprocess.PIPE)
    return [line.strip() for line in io.TextIOWrapper(proc.stdout, encoding='utf-8')]

def parse_gff(GFF_F, tid):
    with open(GFF_f) as input_gene_coords:
        for line in input_gene_coords:
            if line.startswith("#"):
                continue

            line = line.strip().split("\t")

            if line[2] != 'mRNA':
                continue

            gff_name = line[8].split(";")[0].replace("ID=", "")
            if gff_name != tid:
                continue

            start = shift_gene_start(line[3], 500)
            end = shift_gene_end(line[4], 500)

            if not start:
                continue


            
            strand = line[6]
            pos = (line[0], start, end, strand)
            yield pos

def Identify_dogbox(sequence, promotor_boundary = 500, MOTIF = 'ATGCCA'):
    promotor_sequence = sequence[0:promotor_boundary].upper()
    return promotor_sequence.count(MOTIF)

if __name__ == "__main__":
    print("Transcript_id\tDogBox")
    for transcript in transcript_ids:
    
        for position in parse_gff(GFF_f, transcript):
        
            fasta_entry = run_samtools(position, fasta)
            sequence = "".join(fasta_entry[1:])
        
            dogbox_count = Identify_dogbox(sequence)
            print(f"{transcript}\t{dogbox_count}")
