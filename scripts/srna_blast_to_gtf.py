# srna_blast_to_gtf.py
# MWS 3/21/2024
"""
Start with BLASTn output derived from Kroger's sRNA Table S3 as query vs. 17961 chromosome .fa
*Assumes Blast outfmt 6
Read in sRNA lengths also from Kroger Table S3
Option: filter for genes in EG's list
As argument, take sRNA % match threshold and % of length match threshold

Based on thresholding, determine genomic positions of sRNAs in 17961
Output a GTF file for RNA-seq alignment containing sRNAs in 17961.

"""


import argparse
import csv
import pandas as pd


KROGER_SRNA_FILE = "/Users/mws/Documents/geisinger_lab_research/bioinformatics_in_acinetobacter/sRNAs_17978/Kroger_2017_Table_S3_sRNAs.xlsx"
WORKING_DIR = '/Users/mws/Documents/geisinger_lab_research/bioinformatics_in_acinetobacter/sRNAs_17978/stationary_phase_palethorpe_forNER_2024-03-04'
SRNA_BLAST_IN = WORKING_DIR + '/data/srna_blast_output_header.txt'
SRNA_TARGETS = ['sRNA21', 'sRNA85', 'sRNA77', 'sRNA103', 'sRNA40', 'sRNA102', 'sRNA76', 'sRNA20', 'sRNA29', 'sRNA53',
                'sRNA87', 'sRNA30', 'sRNA35', 'sRNA54', 'sRNA84']

def get_args():
    """Return parsed command line arguments."""

    parser = argparse.ArgumentParser(
        description="Write GTF file of sRNAs in 17961",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--filter',  # access with args.filter (i.e. `if args.filter:`
                        help='Pass --filter to filter for an sRNA set (hardcoded).',
                        action='store_true')

    parser.add_argument('-l', '--length_cutoff',  # access with args.length_cutoff
                        metavar="LENGTH_CUTOFF",
                        help='Provide % cutoff for BLAST query sRNA length vs. Kroger 2017.',
                        type=float,
                        default=100)

    parser.add_argument('-i', '--id_cutoff',  # access with args.id_cutoff
                        metavar="ID_CUTOFF",
                        help='Provide % cutoff for BLAST query sRNA percent identity vs. Kroger 2017.',
                        type=float,
                        default=100)

    parser.add_argument('-o', '--output',  # access with args.output
                        metavar='FASTA',
                        help='Provide path and filename for output FASTA file.',
                        type=str,
                        required=True,
                        default=WORKING_DIR+'17961_sRNAs.gtf')

    return (parser.parse_args())

def parse_blast(BLAST_FILEPATH):
    """
    Parse a tab-separated outfmt 6 BLAST file to a dictionary
    Returns a dictionary of dictionary.  keys = query seqid; values = dictionary of each line
    Within subdict, keys = column names from blast outfmt 6
    ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
    """

    blast_dict = {}

    with open(BLAST_FILEPATH, 'r') as blast_infile:
        reader = csv.DictReader(blast_infile, delimiter='\t')
        for line in reader:
            if line['qseqid'] not in blast_dict.keys():
                blast_dict[line['qseqid']] = line
            elif float(line['bitscore']) > float(blast_dict[line['qseqid']]['bitscore']):
                # Check bitscore of current line against the existing data for this srna and overwrite if higher
                blast_dict[line['qseqid']] = line

    return blast_dict


def main():
    # args = get_args()
    temp_filter = True # TODO remove me
    blast_srna_dict = parse_blast(SRNA_BLAST_IN)

    # Code to filter sRNAs to designed list
    keys_to_remove = []
    if temp_filter:  # TODO replace with if args.filter:
        for srna, dict in blast_srna_dict.items():
            if srna not in SRNA_TARGETS:  # SRNA_TARGETS is set at top of file
                keys_to_remove.append(srna)
        for srna in keys_to_remove:
            blast_srna_dict.pop(srna, None)
    kroger_srna_df = pd.read_excel(KROGER_SRNA_FILE, skiprows=1)
    kroger_srna_df_sRNAnum_index = kroger_srna_df.set_index('Name')

    # Filter srnas against provided thresholds
    srnas_passed_threshold_list = []
    for srna, dict in blast_srna_dict.items():
        if float(blast_srna_dict[srna]['length']) > 0.99 * float(kroger_srna_df_sRNAnum_index.loc[srna]['Size']):
            if float(blast_srna_dict[srna]['pident']) > 80:
                srnas_passed_threshold_list.append(srna)
    # print(f'Applied cutoffs <{args.length_cutoff}> for % of length and <{args.id_cutoff}> for % identity')
    print(f'Out of {len(SRNA_TARGETS)} targets specified, {len(srnas_passed_threshold_list)} passed cutoff')  # TODO remove this line.
    # TODO: change the number thresholds here to arguments from args

    # Make lists containing attributes for GTF file
    seqnames = []
    sources = []
    features = []
    starts = []
    ends = []
    scores = []
    strands = []
    frames = []
    attributes = []
    for srna in srnas_passed_threshold_list:
        # compose attribute string
        # make a dictionary containing important attributes
        temp_attr_dict = {
            'transcript_id': srna,
            'gene_id': srna,
            'gene_name': srna,
            'status': 'verified by N blot' if srna in ['sRNA17', 'sRNA37', 'sRNA77', 'sRNA84', 'sRNA100'] else 'predicted',
            'genomic_location': kroger_srna_df_sRNAnum_index.loc[srna]['Genomic location'] if kroger_srna_df_sRNAnum_index.loc[srna]['Genomic location'] != "3'UTR" else "3-UTR",
            'adjacent_genes': f"{kroger_srna_df_sRNAnum_index.loc[srna]['Upstream coding gene']}/{kroger_srna_df_sRNAnum_index.loc[srna]['Downstream coding gene']}",
            'BLASTn_percent_identity': blast_srna_dict[srna]['pident'],
            'BLASTn_evalue': blast_srna_dict[srna]['evalue'],
        }

        # print dictionary to a string, formatted for gtf attributes (key "value"; ...)
        temp_attrs = []
        for attr, value in temp_attr_dict.items():
            temp_attrs.append(f'{attr} "{value}"')
        attr_string = "; ".join(temp_attrs)

        seqnames.append('CP065432.1')
        sources.append('mws2024-03-22')
        features.append('sRNA')
        starts.append(blast_srna_dict[srna]['sstart'] if int(blast_srna_dict[srna]['sstart']) < int(blast_srna_dict[srna]['send']) else blast_srna_dict[srna]['send'])
        ends.append(blast_srna_dict[srna]['send'] if int(blast_srna_dict[srna]['send']) > int(blast_srna_dict[srna]['sstart']) else blast_srna_dict[srna]['sstart'])
        scores.append('.')
        strands.append(kroger_srna_df_sRNAnum_index.loc[srna]['Strand'])
        frames.append('.')
        attributes.append(attr_string)

    # Write GTF file.
    # make dictionary from field lists
    dict_to_df = {
        'seqname': seqnames,
        'source': sources,
        'feature': features,
        'start': starts,
        'end': ends,
        'score': scores,
        'strand': strands,
        'frame': frames,
        'attribute': attributes
    }

    # Write GTF file using pandas
    outfile_path = WORKING_DIR + '/srnas_17961_24-03-22.gtf'
    # quoting/quotechar/escapechar below are needed to ensure the attributes field is properly double-quoted
    pd.DataFrame(dict_to_df).to_csv(outfile_path, sep='\t', header=False, index=False, quoting=csv.QUOTE_NONE, quotechar="",  escapechar="\\")
    print(f'Wrote file to: {outfile_path}')

if __name__ == '__main__':
    main()
