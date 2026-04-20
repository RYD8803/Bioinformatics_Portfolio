#Note: This were an updated code of the previous code for parsing gene
#This function can be called to extract many genes (using lists) as shown below
#This code needs BioPython, so it should be installed first

#Install Biopython
pip install biopython

#Check the version installed
import Bio
print(Bio.__version__)

def extract_many_adj_sequence(accession_list, output_path='', gene=None):
    if gene is None:
        gene = []
    import time
    from Bio import Entrez, SeqIO
    import os

    Entrez.email = "rafaelangeloyudhistira@gmail.com"
    
    while True:
        if not output_path:
            output_path = input("Enter output FASTA file path: ")
            continue

        directory = os.path.dirname(output_path)

        if directory == "" or os.path.exists(directory):
            break
        else:
            print("Invalid directory, try again.")
            output_path = input("Enter output FASTA file path: ")

    success = 0
    failures = 0
    
    with open(output_path, "w") as file:
        for acc in accession_list:
            gene_positions = []
            time.sleep(0.5)

            try:
                stream = Entrez.efetch(
                    db="nucleotide",
                    id=acc,
                    rettype="gb",
                    retmode="text"
                )

                record = SeqIO.read(stream, "genbank")
                stream.close()

                # Default for full sequence
                seq = record.seq
                description = record.description
                header = acc

                # If gene specified → try extract
                regions = []
                if gene:
                    for feature in record.features:
                        if feature.type == "gene":
                            gene_name = None
                            for qualifier in ["gene", "locus_tag"]:
                                if qualifier in feature.qualifiers: 
                                    gene_name = feature.qualifiers[qualifier][0]
                                    break
                            if gene_name:
                                gene_positions.append((
                                    int(feature.location.start),
                                    int(feature.location.end),
                                    gene_name
                                ))
                    gene_positions.sort(key=lambda x: x[0])
                    gene_lower = [g.lower() for g in gene]
                    # Extract intergenic regions between the specified genes
                    for i in range(len(gene_positions)):
                        current_gene = gene_positions[i][2]
                        if current_gene.lower() in gene_lower:

                        # --- previous → current ---
                            if i > 0:
                                prev_gene = gene_positions[i - 1]
                                start = prev_gene[1]
                                end = gene_positions[i][0]
                                if start != end:
                                    region_seq = record.seq[start:end]
                                    regions.append((prev_gene[2], current_gene, region_seq))

                        # --- current → next ---
                            if i < len(gene_positions) - 1:
                                next_gene = gene_positions[i + 1]
                                start = gene_positions[i][1]
                                end = next_gene[0]
                                if start < end:
                                    region_seq = record.seq[start:end]
                                    regions.append((current_gene, next_gene[2], region_seq))

                # Writing the sequence
                description = record.description.replace(" ", "_")
                if "," in description:
                    description = description.split(",")[0]
                
                for idx, (gene1, gene2, region) in enumerate(regions):
                    header = f"{acc}_{gene1}-{gene2}"
                    file.write(f">{header}|{description}|region between {gene1} - {gene2} \n{region}\n")
                    print(f"Record: {header}: {description}_region between {gene1} - {gene2}")
                success += 1

            except Exception as e:
                print(f"Failed for {acc}: {e}")
                failures += 1
                continue

    print("=" * 50 + "\n")
    print(f"Finished fetching:\nSuccess: {success}\nFailures: {failures}")

Accessions = ["NC_063512.1", "NC_065014.1", "NC_070310.1", "NC_073119.1", "NC_073117.1",
"NC_056110.1", "NC_070315.1", "NC_061410.1", "NC_067030.1", "PQ675783.1", "NC_073120.1",
"NC_046385.1", "NC_063513.1", "NC_065245.1", "NC_070321.1", "PQ572754.1", "NC_068748.1",
"NC_073118.1", "PQ619426.1", "PQ675781.1", "OP618127.1", "NC_073122.1", "OR288087.1", "NC_045096.1",
"NC_088496.1", "NC_070330.1", "NC_073123.1", "PQ619425.1", "NC_073121.1", "NC_047450.1", "PQ675782.1",
"MZ580429.1", "NC_070304.1", "NC_070305.1", "NC_070306.1", "NC_070311.1", "NC_070312.1", "NC_070313.1",
"NC_070314.1", "NC_070317.1", "NC_070318.1", "NC_070319.1", "NC_070322.1", "NC_070323.1", "NC_070324.1",
"NC_070327.1", "NC_070329.1", "NC_070331.1", "NC_070332.1", "NC_070333.1", "NC_070334.1", "NC_070307.1",
"NC_070308.1", "NC_070309.1", "NC_070316.1", "NC_070320.1", "NC_070325.1", "NC_070326.1", "NC_070328.1",
"NC_000932.1", "NC_007144.1"]

genes = ["rpoC1", "rpoB", "psbE", "psbK", "petA", "rpl2", "rpl22", "rps12"]

# Create folder
import os
folder_name = "C:\\Users\\Regina Maria\\OneDrive\\Desktop\\fasta_files_sequences"
os.makedirs(folder_name, exist_ok=True)

for gene_name in genes:
    output = os.path.join(folder_name, f"Begonia_{gene_name}.fa")
    print('\033[1m' + gene_name.upper() + '\033[0m')
    
    extract_many_adj_sequence(
        Accessions,
        output,
        [gene_name]   
    )