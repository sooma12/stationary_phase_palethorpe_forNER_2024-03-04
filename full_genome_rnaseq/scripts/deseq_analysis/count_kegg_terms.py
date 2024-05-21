# count_kegg_terms.py
# mws 2024-04-12
"""
Count occurrences of KEGG terms in a DESeq output dataset (first column of data = ACX60 locus tags)
Output table format:
KEGG_ID   KEGG_term_name   counts
Usage:
"""



import argparse, csv, sys
import pandas as pd


BASE_DIR = '//'
FUNAGE_ANNOT_FILE = BASE_DIR + 'REFERENCE/Acinetobacter_baumannii_ATCC_17978-mff_ASM107767v1_genome2d.GO'
FXN_ANNOT_FILE = BASE_DIR + 'REFERENCE/17978-mff_functional_annotations.csv'
GO_TERM2NAME_FILE = BASE_DIR + 'REFERENCE/go_term2name.csv'

def main():
    args = get_args()  # infile, condition, padj, outfile
    condition_call = args.condition.lower()
    if condition_call not in ['up', 'down']:
        print('Error: Condition must be either "up" or "down".')
        sys.exit(1)
    functional_anno_dict = parse_functional_annos()

    locus_tag_list = get_locus_list(deseq_file=args.infile, padj_threshold=args.padj, up_or_down=condition_call, force_bool=args.force)

    kegg_count_dict = {}
    kegg_not_found = []
    kegg_success_counter = 0
    for locus in locus_tag_list:
        try:
            kegg_term_list = [functional_anno_dict[locus]['Category1'], functional_anno_dict[locus]['Category2'], functional_anno_dict[locus]['Category3']]
            kegg_success_counter += 1
            for term in kegg_term_list:
                kegg_count_dict[term] = kegg_count_dict.get(term, 0) + 1
        except KeyError as e:
            kegg_not_found.append(e.args[0])

    # List of lists for Pandas
    found_terms = kegg_count_dict.keys()
    lol = []
    for term in found_terms:
        current_list = []
        if term:
            current_list.append(term)
            try:
                current_list.append(int(kegg_count_dict[term]))
            except KeyError:
                pass

        lol.append(current_list)

    # Use Pandas to make dataframe from lists and write as Excel
    df = pd.DataFrame(lol, columns=['KEGG term', 'count'])
    df = df.sort_values(by=['count'], ascending=False)
    df = df.dropna()
    df.to_excel(args.outfile)

    print(f'FUNAGE annotations: found {kegg_success_counter} loci but {len(kegg_not_found)} missing')


def get_args():
    """Return parsed CL arguments"""

    parser = argparse.ArgumentParser(
        description="Count KEGG annotations using EG functional annotations spreadsheet",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-infile',  # access with args.infile
        help="Provide filepath of input file (DESeq2 output, with 1st column locus tags, 3rd column log2FoldChange, and 6th column adjusted p value",
        type=str)

    parser.add_argument('-condition',  # access with args.condition
        help="Provide regulation condition to analyze. Pass `-condition up` for upregulated or `-condition down` for downregulated",
        type=str)

    parser.add_argument('--padj',  # access with args.padj
        help="Provide adjusted p value threshold to include",
        type=float,
        default=0.1)

    parser.add_argument('--outfile',  #access with args.outfile
        help="Provide name for output file (include .xlsx file extension, e.g. `out.xlsx`)",
        type=str,
        default='out.xlsx')

    parser.add_argument(
        '--force',  #access with args.force
        help="Pass True to force program to read 3rd column as log2FoldChange and 5th column as padj despite column headers",
        type=bool,
        default=False)

    return parser.parse_args()


def parse_functional_annos():
    """
    Open and parse EG's 17978-mff functional annotations
    Returns a dictionary of dictionaries with ACX60 locus tags as keys, and sub-dictionaries containing all other columns from this spreadsheet
    Access a specific piece of information with, e.g., `functional_anno_dict['ACX60_RS15100']['Gene ontology IDs']
    Columns: ['', 'Gene', 'Protein ID', 'K Number', 'Global Gene', 'Global Product', 'Tag1', 'Tag2', 'Tag3', 'Category1', 'Category2', 'Category3', 'Global Path', 'Entry', 'query', 'Entry name', 'Protein names', 'Gene names', 'Organism', 'Length', 'Gene names  (ordered locus )', 'Gene names  (primary )', 'EC number', 'Gene ontology (biological process)', 'Gene ontology (cellular component)', 'Gene ontology (GO)', 'Gene ontology (molecular function)', 'Gene ontology IDs', 'Organism ID']
    """

    functional_anno_dict = {}
    with open(FXN_ANNOT_FILE, 'r') as fxn_annos_infh:
        anno_dictreader = csv.DictReader(fxn_annos_infh)
        for line in anno_dictreader:
            functional_anno_dict[line['Gene']] = line

    return functional_anno_dict


def get_locus_list(deseq_file, padj_threshold, up_or_down, force_bool) -> list:
    """ Parse DESeq output and return a list of locus tags matching log2FoldChange and adjusted p value threshold
    """

    locus_list_to_count = []

    with open(deseq_file, 'r') as deseq_infh:
        deseq_reader = csv.reader(deseq_infh)
        # Validate DESeq2 format
        first_line = next(deseq_reader)
        if first_line[2] != 'log2FoldChange' or first_line[5] != 'padj' and force_bool == False:
            print('Error: File does not seem to match DESeq2 output format.')
            print('Make sure 3rd column is log2FoldChange and 6th column is adjusted p value, then pass --force True')
            sys.exit(1)

        for line in deseq_reader:
            if line[2] == 'NA' or line[5] == 'NA':
                continue  # skip NAs
            if up_or_down == 'up':
                if (float(line[5]) < padj_threshold) & (float(line[2]) > 0):
                    locus_list_to_count.append(line[0])
            if up_or_down == 'down':
                if (float(line[5]) < padj_threshold) & (float(line[2]) < 0):
                    locus_list_to_count.append(line[0])

        return locus_list_to_count




if __name__ == '__main__':
    main()
