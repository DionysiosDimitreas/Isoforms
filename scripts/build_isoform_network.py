import pandas as pd
from collections import defaultdict
from itertools import product
import gzip

print("Starting the transcript network build process...")

# --- 1. Define Input Filenames ---
# Make sure these filenames match what you have in your directory.
GENE_PAIRS_FILE = "cellchat_rl_pairs.tsv"
ISOFORM_DOMAINS_FILE = "isoform_pfam_domains.txt" # Your pfam_scan.pl output
DOMAIN_INTERACTIONS_FILE = "3did_pfamInteractions_1804.tsv.gz" # Can read .gz directly
BIOMART_FILE = "biomart_export.tsv" # Your file with ENSG, ENST, ENSP, Gene Symbol

# --- 2. Load and Process All Data ---

print("Loading and processing input files...")

# Load gene-level interaction pairs
gene_pairs_df = pd.read_csv(GENE_PAIRS_FILE, sep="\t")

# Load the domain-domain interaction rulebook into a set for fast lookups
domain_interactions = set()
with gzip.open(DOMAIN_INTERACTIONS_FILE, "rt") as f:
    for line in f:
        parts = line.strip().split() 
        if len(parts) >= 2:
            d1, d2 = parts[0], parts[1]
            domain_interactions.add(tuple(sorted((d1, d2))))

# Load and parse the isoform domain data from pfam_scan.pl output
try:
    isoform_domains_raw_df = pd.read_csv(
        ISOFORM_DOMAINS_FILE,
        delim_whitespace=True,
        comment="#",
        header=None,
        usecols=[1, 3],
        names=["pfam_accession", "isoform_id"],
    )
    # Clean version numbers from IDs
    isoform_domains_raw_df["pfam_accession"] = isoform_domains_raw_df["pfam_accession"].str.split('.').str[0]
    isoform_domains_raw_df["isoform_id"] = isoform_domains_raw_df["isoform_id"].str.split('.').str[0]
except (FileNotFoundError, pd.errors.EmptyDataError) as e:
    print(f"Error loading {ISOFORM_DOMAINS_FILE}: {e}")
    exit()

# Load the BioMart ID mapping file
biomart_df = pd.read_csv(BIOMART_FILE, sep="\t")
biomart_df.rename(columns={
    'Gene stable ID': 'ENSG',
    'Transcript stable ID': 'ENST',
    'Protein stable ID': 'ENSP',
    'HGNC symbol': 'gene_symbol'
}, inplace=True)
# Drop rows where essential IDs are missing
biomart_df.dropna(subset=['ENSP', 'ENST', 'gene_symbol'], inplace=True)

print("Stripping version numbers from BioMart file...")
biomart_df['ENSP'] = biomart_df['ENSP'].str.split('.').str[0]
biomart_df['ENST'] = biomart_df['ENST'].str.split('.').str[0]

# --- 3. Create Efficient Lookup Dictionaries ---

print("Creating lookup maps for fast processing...")

# Create a dictionary: {gene_symbol: [isoform1 (ENSP), isoform2 (ENSP), ...]}
gene_to_isoforms = biomart_df.groupby('gene_symbol')['ENSP'].apply(list).to_dict()

# Create a dictionary: {isoform_id (ENSP): {domain1, domain2, ...}}
isoform_to_domains = isoform_domains_raw_df.groupby('isoform_id')['pfam_accession'].apply(set).to_dict()

# *** NEW: Create a map from Protein ID (ENSP) to Transcript ID (ENST) ***
ensp_to_enst_map = pd.Series(biomart_df.ENST.values, index=biomart_df.ENSP).to_dict()

# --- 4. Find Transcript-Specific Interactions ---

print("Identifying transcript-level interactions... This may take a while.")
final_transcript_interactions = []

# Iterate through each gene-level interaction
for _, row in gene_pairs_df.iterrows():
    ligand_gene = row["ligand"]
    receptor_gene = row["receptor"]

    ligand_isoforms = gene_to_isoforms.get(ligand_gene, [])
    receptor_isoforms = gene_to_isoforms.get(receptor_gene, [])

    if not ligand_isoforms or not receptor_isoforms:
        continue

    # This loop will now check ALL combinations
    for lig_iso, rec_iso in product(ligand_isoforms, receptor_isoforms):
        lig_domains = isoform_to_domains.get(lig_iso, set())
        rec_domains = isoform_to_domains.get(rec_iso, set())

        if not lig_domains or not rec_domains:
            continue

        for dom1, dom2 in product(lig_domains, rec_domains):
            if tuple(sorted((dom1, dom2))) in domain_interactions:
                lig_enst = ensp_to_enst_map.get(lig_iso, 'N/A')
                rec_enst = ensp_to_enst_map.get(rec_iso, 'N/A')

                final_transcript_interactions.append({
                    "ligand_gene": ligand_gene,
                    "receptor_gene": receptor_gene,
                    "ligand_transcript": lig_enst,
                    "receptor_transcript": rec_enst,
                })
                # Found an interaction, break the domain-checking loop
                # and move to the next isoform pair.
                break


# --- 5. Save the Final Results ---

print("Saving the final network...")
output_df = pd.DataFrame(final_transcript_interactions).drop_duplicates()
output_df.to_csv("transcript_specific_interactions.tsv", sep="\t", index=False)

print(f"Done! Found {len(output_df)} transcript-level interactions.")
print("Output saved to 'transcript_specific_interactions.tsv'")
