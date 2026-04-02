#Note: This were an updated code of the previous code for parsing gene
#This function can be called to extract many genes (using lists) as shown below
#This code needs BioPython, so it should be installed first

#Do mind that this code extract CDS. 
#If you were intrested in parsing the full gene, you can modify the code below
#line 65 specified the need for CDS

#Install Biopython
pip install biopython

#Check the version installed
import Bio
print(Bio.__version__)

#Extract many sequences
def extract_many_sequence(accession_list, output_path='', gene=''):
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

    print("Fetching " + gene + "...")

    success = 0
    failures = 0

    with open(output_path, "w") as file:
        for acc in accession_list:
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
                if gene:
                    found = False

                    for feature in record.features:
                        if feature.type == "CDS":
                            if "gene" in feature.qualifiers:
                                if gene in feature.qualifiers["gene"]:
                                    seq = feature.extract(record.seq)
                                    header = f"{acc}_{gene}"
                                    found = True
                                    break

                    if not found:
                        print(f"{acc}: gene '{gene}' not found")
                        failures += 1
                        continue

                # Writing the sequence
                description = record.description.replace(" ", "_")
                comma = description.find(",")
                description = description[:comma]
                
                gene_label = gene if gene else "full"
                file.write(f">{header}|{description}|{gene_label}\n{seq}\n")

                print(f"Record: {header}: {description}")
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
    print('\033[1m' + gene_name.upper() + '\033[0m') #header in output file
    extract_many_sequence(Accessions, output, gene_name)