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
