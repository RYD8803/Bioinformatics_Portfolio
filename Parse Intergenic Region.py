from Bio import Entrez, SeqIO

# Function to fetch a GenBank record
def fetch_genbank_record(accession):
    Entrez.email = "rafaelangeloyudhistira@example.com"  # Always provide your email
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record

# Function to extract intergenic regions between specified genes
def extract_intergenic_regions(record, gene_list):
    intergenic_regions = []
    gene_positions = []

    # Collect the positions of the specified genes
    for feature in record.features:
        if feature.type == "gene":
            gene_name = None
            for qualifier in ["gene", "locus_tag"]:
                if qualifier in feature.qualifiers:
                    gene_name = feature.qualifiers[qualifier][0]
                    break
            if gene_name in gene_list:
                gene_positions.append((feature.location.start, feature.location.end, gene_name))

    # Sort the gene positions by their start positions
    gene_positions.sort()

    # Extract intergenic regions between the specified genes
    for i in range(len(gene_positions) - 1):
        start = gene_positions[i][1]
        end = gene_positions[i + 1][0]
        if start < end:
            intergenic_region = record.seq[start:end]
            intergenic_regions.append((gene_positions[i][2], gene_positions[i + 1][2], intergenic_region))

    return intergenic_regions

# List of accession numbers for different species
accessions = ["NC_063512.1", "NC_065014.1", "NC_070310.1", "NC_073119.1", "NC_073117.1", "NC_056110.1", "NC_070315.1", "NC_061410.1", "NC_067030.1", "NC_073120.1", "NC_046385.1", "NC_063513.1", "NC_065245.1", "NC_070321.1", "NC_068748.1", "NC_073118.1", "NC_073122.1", "NC_071962.1", "NC_045096.1", "NC_070330.1", "NC_073123.1", "NC_073121.1", "NC_047450.1"]  # Add your accession numbers here

# Corresponding species names for the accession numbers
species_names = {
    "NC_063512.1": "B. arachnoidea",
    "NC_065014.1": "B. asteropyrifolia",
    "NC_070310.1": "B. bracteata",
    "NC_073119.1": "B. cathayana",
    "NC_073117.1": "B. cavalerei",
    "NC_056110.1": "B. coptidifolia",
    "NC_070315.1": "B. dipetala",
    "NC_061410.1": "B. emeiensis",
    "NC_067030.1": "B. ferox",
    "NC_073120.1": "B. grandis",
    "NC_046385.1": "B. guangxiensis",
    "NC_063513.1": "B. gulongshanensis",
    "NC_065245.1": "B. handelii",
    "NC_070321.1": "B. henryi",
    "NC_068748.1": "B. jinyuensis",
    "NC_073118.1": "B. leprosa",
    "NC_073122.1": "B. obsolescens",
    "NC_071962.1": "B. pedatifida",
    "NC_045096.1": "B. pulchrifolia",
    "NC_070330.1": "B. samhaensis",
    "NC_073123.1": "B. smithiana",
    "NC_073121.1": "B. umbraculifolia",
    "NC_047450.1": "B. versicolor"
}

# List of genes of interest
genes_of_interest = ["rpoC2", "rpoC1"]  # Replace with actual gene names

# Fetch records and extract intergenic regions
for accession in accessions:
    print(f"Processing accession: {accession}")
    record = fetch_genbank_record(accession)
    intergenic_regions = extract_intergenic_regions(record, genes_of_interest)
    species_name = species_names.get(accession, "Unknown species")
    print(f"Intergenic regions for {species_name} ({accession}):")
    if not intergenic_regions:
        print("No intergenic regions found between the specified genes.")
    for idx, (gene1, gene2, region) in enumerate(intergenic_regions):
        print(f">{species_name}|{accession}|region between {gene1} - {gene2} \n{region}\n")

import os

# Folder to store the output files
output_folder = "rpoC2-rpoC1_output"
os.makedirs(output_folder, exist_ok=True)

# Fetch records and extract intergenic regions, writing each to a separate file
for accession in accessions:
    print(f"Processing accession: {accession}")
    record = fetch_genbank_record(accession)
    intergenic_regions = extract_intergenic_regions(record, genes_of_interest)

    species_name = species_names.get(accession, "Unknown species")  # Move this line inside the loop

    output_file_path = os.path.join(output_folder, f"{species_name} {accession}_rpoC2-rpoC1.txt")
    with open(output_file_path, "w") as output_file:
        if not intergenic_regions:
            output_file.write("No intergenic regions found between the specified genes.\n")
        else:
            for idx, (gene1, gene2, region) in enumerate(intergenic_regions):
                output_file.write(f">{species_name}|{accession}|region between {gene1} - {gene2} \n{region}\n")

    print(f"Output written to {output_file_path}")

# Create a zip file
import subprocess
subprocess.run(["zip", "-r", "rpoC2-rpoC1_output.zip", "rpoC2-rpoC1_output/"])