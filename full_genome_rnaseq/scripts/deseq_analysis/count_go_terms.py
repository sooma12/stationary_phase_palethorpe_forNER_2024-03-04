# count_go_terms.py
# mws 2024-04-12
"""
Count occurrences of GO terms in a DESeq output dataset (first column of data = ACX60 locus tags)
Show counts using FUNAGE-Pro annotations for 17978-mff as well as EG functional annotation spreadsheet
Output table format:
GO_ID    GO_term_name   counts_funage   gene_list_funage    counts_fxn_annots   gene_list_fxn_annots
Usage: e.g. `python scripts/count_go_terms.py -infile ./DES_lpsB_2024-04-04.csv --outfile goterms_dlpsB_up.xlsx -condition up`
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
    funage_tag_to_id, funage_id_to_name = parse_funage_annos()
    functional_anno_dict = parse_functional_annos()
    term2name = _get_term2name_dict()

    locus_tag_list = get_locus_list(deseq_file=args.infile, padj_threshold=args.padj, up_or_down=condition_call, force_bool=args.force)

    funage_go_count_dict = {}
    fxn_annos_go_count_dict = {}
    funage_not_found = []
    fxn_not_found = []
    funage_success_counter = 0
    fxn_success_counter = 0
    for locus in locus_tag_list:
        try:
            funage_go_ids = funage_tag_to_id[locus]
            funage_success_counter += 1
            for id in funage_go_ids.split(','):
                funage_go_count_dict[id] = funage_go_count_dict.get(id, 0) + 1
        except KeyError as e:
            funage_not_found.append(e.args[0])

        try:
            fxn_annos_go_ids = functional_anno_dict[locus]['Gene ontology IDs'].split('; ')
            fxn_success_counter += 1
            for id in fxn_annos_go_ids:
                fxn_annos_go_count_dict[id] = fxn_annos_go_count_dict.get(id, 0) + 1
        except KeyError as e:
            fxn_not_found.append(e.args[0])

    # List of lists for Pandas
    found_ids = set(list(fxn_annos_go_count_dict.keys()) + list(funage_go_count_dict.keys()))
    lol = []
    for id in found_ids:
        current_list = []
        if id != None:
            current_list.append(id)
            try:
                current_list.append(term2name[id])
            except KeyError:
                try:
                    current_list.append(funage_id_to_name[id])
                except KeyError:
                    pass
            try:
                current_list.append(int(funage_go_count_dict[id]))
            except KeyError:
                pass
            try:
                current_list.append(int(fxn_annos_go_count_dict[id]))
            except KeyError:
                pass
        lol.append(current_list)

    # Use Pandas to make dataframe from lists and write as Excel
    df = pd.DataFrame(lol, columns=['GO ID', 'GO term name', 'FUNAGE count', 'EG functional annotations count'])
    df = df.sort_values(by=['FUNAGE count'], ascending=False)
    df = df.dropna(subset=['FUNAGE count', 'EG functional annotations count'], how='all')  # Drop rows with NaNs in both count columns
    df.to_excel(args.outfile)

    print(f'FUNAGE annotations: found {funage_success_counter} loci but {len(funage_not_found)} missing')
    print(f'EG functional annotations: found {fxn_success_counter} loci but {len(fxn_not_found)} missing')


def get_args():
    """Return parsed CL arguments"""

    parser = argparse.ArgumentParser(
        description="Count GO annotations using FUNAGE-Pro and EG functional annotations spreadsheets",
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


def parse_funage_annos():
    """
    Open and parse FUNAGE-Pro 17978-mff annotation file
    Return two dictionaries as a tuple:
    1. dictionary of 'ACX60 locus tag':'GO_IDs' pairs; may be >1 GO_ID separated by commas
    2. dictionary of 'GO_ID':'GO_term_name' pairs
    """

    with open(FUNAGE_ANNOT_FILE, 'r') as funage_infh:
        funage_reader = csv.reader(funage_infh, delimiter='\t')

        funage_tag_id_dict = {}
        funage_id_name_dict = {}

        for line in funage_reader:  # line[0] = ACX60 tag; line[1] = GO ID; line[2] = GO term name
            # Build locus_tag:go_id dictionary
            if line[0] not in funage_tag_id_dict:
                funage_tag_id_dict[line[0]] = line[1]
            else:
                current_id = funage_tag_id_dict[line[0]]
                updated_id = current_id + ',' + line[1]
                funage_tag_id_dict[line[0]] = updated_id

            # build dictionary of GO id:name pairs
            if line[1] not in funage_id_name_dict:
                funage_id_name_dict[line[1]] = line[2]

        return funage_tag_id_dict, funage_id_name_dict


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


def _get_term2name_dict():
    """ Return a dictionary of go ID:term name pairs.  Uses the whole GO database """
    term2name_dict = {}
    with open(GO_TERM2NAME_FILE, 'r') as infile:
        reader = csv.reader(infile)

        for line in reader:
            term2name_dict[line[0]] = line[1]

    return term2name_dict


if __name__ == '__main__':
    main()
