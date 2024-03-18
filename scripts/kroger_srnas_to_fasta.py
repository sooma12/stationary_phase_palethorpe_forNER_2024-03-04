# kroger_srnas_to_fasta.py
# MWS 3/18/24
"""
Take Kroger 2017 sRNA table; output a fasta in the format:
>sRNA #
<sequence>
"""


import pandas as pd


KROGER_SRNA_FILE = "/Users/mws/Documents/geisinger_lab_research/bioinformatics_in_acinetobacter/sRNAs_17978/Kroger_2017_Table_S3_sRNAs.xlsx"
OUTFILE = "/Users/mws/Documents/geisinger_lab_research/bioinformatics_in_acinetobacter/sRNAs_17978/stationary_phase_palethorpe_forNER_2024-03-04/17978-srnas.fasta"


kroger_srna_df = pd.read_excel(KROGER_SRNA_FILE, skiprows=1)
kroger_srna_df_sRNAnum_index = kroger_srna_df.set_index('Name')

sRNA_seq_dict = {}
for srna_id in kroger_srna_df.to_dict("list")['Name']:
    sRNA_seq_dict[srna_id] = kroger_srna_df_sRNAnum_index.loc[srna_id, "Sequence 5' -> 3'"]

fasta_entries = []
for srna, seq in sRNA_seq_dict.items():
    fasta_entry = f">{srna}\n{seq}\n"
    fasta_entries.append(fasta_entry)

with open(OUTFILE, 'w') as outfile:
    for entry in fasta_entries:
        outfile.write(entry)
