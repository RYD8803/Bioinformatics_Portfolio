from Bio import SeqIO, Entrez
from Bio.SeqFeature import SeqFeature, FeatureLocation
import os
import zipfile

# Set your email for NCBI Entrez
Entrez.email = "rafaelangeloyudhistira@gmail.com"

# Mapping accession numbers to species names
species_names = {
    "NC_063512.1": "Begonia arachnoidea",
    "NC_065014.1": "Begonia asteropyrifolia",
    "NC_070310.1": "Begonia bracteata voucher Peng23521",
    "NC_073119.1": "Begonia cathayana",
    "NC_073117.1": "Begonia cavaleriei",
    "NC_056110.1": "Begonia coptidifolia",
    "NC_070315.1": "Begonia dipetala voucher Peng22521",
    "NC_061410.1": "Begonia emeiensis",
    "NC_067030.1": "Begonia ferox",
    "PQ675783.1": "Begonia filiformis",
    "NC_073120.1": "Begonia grandis",
    "NC_046385.1": "Begonia guangxiensis plastid",
    "NC_063513.1": "Begonia gulongshanensis",
    "NC_065245.1": "Begonia handelii",
    "NC_070321.1": "Begonia henryi voucher RBGE20141517",
    "PQ572754.1": "Begonia jingxiensis",
    "NC_068748.1": "Begonia jinyunensis",
    "NC_073118.1": "Begonia leprosa",
    "PQ619426.1": "Begonia liuyanii",
    "PQ675781.1": "Begonia luochengensis",
    "OP618127.1": "Begonia masoniana",
    "NC_073122.1": "Begonia obsolescens",
    "OR288087.1": "Begonia pedatifida",
    "NC_045096.1": "Begonia pulchrifolia",
    "NC_088496.1": "Begonia retinervia",
    "NC_070330.1": "Begonia samhaensis voucher RBGE19990398",
    "NC_073123.1": "Begonia smithiana",
    "PQ619425.1": "Begonia subcoriacea",
    "NC_073121.1": "Begonia umbraculifolia",
    "NC_047450.1": "Begonia versicolor",
    "PQ675782.1": "Begonia wangii",
    "MZ580429.1": "Begonia sp. RBGE:20160139 voucher RBGE20160139",
    "NC_070304.1" : "Begonia aconitifolia voucher Peng22224",
    "NC_070305.1" : "Begonia amoeboides voucher RBGE20180924",
    "NC_070306.1" : "Begonia anemoniflora voucher RBGE20160123",
    "NC_070311.1" : "Begonia buddleiifolia voucher RBGE20160126",
    "NC_070312.1" : "Begonia convolvulacea voucher Peng21267",
    "NC_070313.1" : "Begonia cubensis voucher Peng21285",
    "NC_070314.1" : "Begonia depauperata voucher Peng24271",
    "NC_070317.1" : "Begonia egregia voucher Peng23327",
    "NC_070318.1" : "Begonia fissistyla voucher Peng21417",
    "NC_070319.1" : "Begonia foliosa voucher Peng21254",
    "NC_070322.1" : "Begonia heydei voucher RBGE20131992",
    "NC_070323.1" : "Begonia karwinskyana voucher Peng20880",
    "NC_070324.1" : "Begonia ludwigii voucher Peng22333",
    "NC_070327.1" : "Begonia oaxacana voucher RBGKew (s.n.)",
    "NC_070329.1" : "Begonia rossmanniae voucher RBGE20151093",
    "NC_070331.1" : "Begonia sanguinea voucher Peng21284",
    "NC_070332.1" : "Begonia santos-limae voucher Peng21320",
    "NC_070333.1" : "Begonia ulmifolia voucher RBGE20030607",
    "NC_070334.1" : "Begonia undulata voucher Peng21275",
    "NC_070307.1" : "Begonia anisosepala voucher Peng24894",
    "NC_070308.1" : "Begonia baccata voucher K023612",
    "NC_070309.1" : "Begonia bogneri voucher Peng22541",
    "NC_070316.1" : "Begonia dregei voucher K19503",
    "NC_070320.1" : "Begonia henrilaportei voucher RBGE20160414",
    "NC_070325.1" : "Begonia meyeri-johannis voucher RBGE2013122",
    "NC_070326.1" : "Begonia microsperma voucher Peng20259",
    "NC_070328.1" : "Begonia oxyloba voucher RBGE19982761"
}

# Function to fetch GenBank record
def fetch_genbank_record(accession):
    try:
        with Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text") as handle:
            return SeqIO.read(handle, "genbank")
    except Exception as e:
        print(f"Error fetching {accession}: {e}")
        return None

# Function to extract rpoC1 CDS
def extract_rpoc1_cds(record, species_name, accession):
    if record is None:
        return None

    for feature in record.features:
        if feature.type == "CDS" and "gene" in feature.qualifiers:
            if "rpoc1" in feature.qualifiers["gene"][0].lower():  # Ensure case-insensitive match
                cds_sequence = feature.extract(record.seq)

                # Define output path
                output_file = os.path.join(output_folder, f"{species_name.replace(' ', '_')}_{accession}|rpoC1_CDS.txt")

                # Write sequence to file
                with open(output_file, "w") as out_handle:
                    out_handle.write(f">{species_name}|{accession}|rpoC1_CDS\n{cds_sequence}\n")

                return output_file
    print(f"rpoC1 CDS not found in {accession}")
    return None

# Folder for output files
output_folder = "rpoC1_CDS_output"
os.makedirs(output_folder, exist_ok=True)

# Fetch and extract rpoC1 CDS
for accession, species_name in species_names.items():
    print(f"Processing: {accession} ({species_name})")
    record = fetch_genbank_record(accession)

    output_file_path = extract_rpoc1_cds(record, species_name, accession)

    if output_file_path:
        print(f"Extracted rpoC1 CDS: {output_file_path}")
    else:
        print(f"No rpoC1 CDS found for {species_name} ({accession})")

# Create a zip file for all CDS sequences
zip_filename = "rpoC1_CDS_output.zip"
with zipfile.ZipFile(zip_filename, "w") as zipf:
    for file in os.listdir(output_folder):
        zipf.write(os.path.join(output_folder, file), file)

print(f"All extracted CDS sequences zipped into {zip_filename}")
