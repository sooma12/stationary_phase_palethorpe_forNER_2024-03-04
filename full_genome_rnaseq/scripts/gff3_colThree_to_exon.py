# gff3_colThree_to_exon.py
"""
Change cds ONLY to exon in gff3 genome file for STAR genomeGenerate

gff3 file has the following terms in the 3rd column (index 2):
{'region': 1, 'gene': 3628, 'CDS': 3591, 'rRNA': 18, 'exon': 94, 'tRNA': 72,
'pseudogene': 55, 'sequence_feature': 1, 'riboswitch': 7, 'SRP_RNA': 1,
'ncRNA': 1, 'tmRNA': 1, 'RNase_P_RNA': 1}

USAGE:
python3 gff3_colThree_to_exon.py <infile path> <outfile path> <3rd column string to change>
"""

GFF_FILE_IN = '/work/geisingerlab/Mark/REF_GENOMES/17978-mff/NZ_CP012004.gff3'
GFF_EXON_FILE_OUT = '/work/geisingerlab/Mark/rnaSeq/2024-01_rnaseq_pbpGlpsB/ref/NZ_CP012004_col3ex.gff3'

import argparse

def get_args():
    """Return parsed CL arguments"""

    parser = argparse.ArgumentParser(
        description="Convert 'CDS' in 3rd column of infile to 'exon'",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument(
        'infile',  # access with args.infile
        help="Provide filepath of input file, e.g. a gff3 or gtf file;"
             "defaults to NZ_CP012004.gff3",
        type=str,
        default=GFF_FILE_IN)

    parser.add_argument(
        'outfile',  #access with args.outfile
        help="Provide filepath for output file, e.g. a gtf file;"
             "defaults to rnaSeq/2024-01_rnaseq_pbpGlpsB/ref/NZ_CP012004_col3ex.gff3",
        type=str,
        default=GFF_EXON_FILE_OUT)

    parser.add_argument(
        'term_2_exon',  #access with args.term_2_exon
        help="Provide the string in the 3rd column to convert to exon;"
             "defaults to 'CDS'",
        type=str,
        default='CDS')

    return(parser.parse_args())

def convert_third_column(args):
    with open(args.infile, 'r') as infh, open(args.outfile, 'w') as out_fh:
        # TODO_temp_coldict = {}
        for line in infh:
            split_line = line.rstrip().split(sep='\t')
            # if split_line[0] == 'NZ_CP012004.1':
            #     item = split_line[2]
            #     if item in TODO_temp_coldict:
            #         TODO_temp_coldict[item] += 1
            #     else:
            #         TODO_temp_coldict[item] = 1
            new_line = split_line
            if new_line[0] == 'NZ_CP012004.1':
                if new_line[2] == args.term_2_exon:
                    new_line[2] = 'exon'
            line_to_write = '\t'.join(new_line)
            out_fh.write(f'{line_to_write}\n')

if __name__ == "__main__":
    convert_third_column(get_args())
