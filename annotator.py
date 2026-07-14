from Bio import Entrez
import vcf
import yaml, os, json
from dotenv import load_dotenv
import pandas as pd

load_dotenv()  # Load environment variables from .env file
Entrez.email = "liliannliuu@gmail.com" # module-level attributes for verification purposes
Entrez.api_key = os.getenv("NCBI_API_KEY")

def name_standardize(name): # helper function to standardize names for comparison
    stand_name = name.lower().replace("_", " ").replace("-", " ").strip()
    if "/" in stand_name:
        stand_name = stand_name.split("/")[0].strip()  # take the first part before the slash
    return stand_name

def parse(vcf_path): # reads vcf file, storing the variant information as a list of Record objects
    vcf_reader = vcf.Reader(open(vcf_path, 'r'))
    variant_info = []
    for record in vcf_reader:
        variant_info.append({
            "chrom": record.CHROM,
            "pos": record.POS,
            "ref": record.REF,
            "alt": record.ALT,
            "id": record.ID})
    return variant_info

def annotate(variant_info): # retrieves clinical significance and molecular consequence from ClinVar using NCBI Entrez API
    for variant in variant_info:
        # fetching ncbi ids
        handle = Entrez.esearch(db="clinvar", term=variant["id"], retmax=1) # esearch returns a list of ncbi ids matching the search term
        record = Entrez.read(handle) # record is a biopython object that is essentially a dictionary
        handle.close()
        ncbi_id = record["IdList"][0] if record["IdList"] else None # assigns the first ncbi id in the list to ncbi_id, or None if the list is empty
        variant['ncbi_id'] = ncbi_id

        # getting clinvar data using ncbi ids 
        handle = Entrez.esummary(db="clinvar", id=variant["ncbi_id"], retmode="json") # esummary returns info on the variant in JSON format
        data = json.loads(handle.read().decode()) # loads the JSON data from the response into a Python dictionary
        handle.close()
        
        # access & store clinical significance
        if 'result' not in data or ncbi_id not in data['result']: # accounts for variants without ClinVar entry
            variant['clin_sig'] = "Unknown"
            variant['mole_conseq'] = ["Unknown"] # this is a list to account for molecular cons being in list format in the ClinVar API
            continue

        elif 'germline_classification' in data['result'][ncbi_id]:
            variant['clin_sig'] = data['result'][ncbi_id]['germline_classification']['description']
        else:
            variant['clin_sig'] = "Unknown"

        # access & store functional consequence
        if 'molecular_consequence_list' in data['result'][ncbi_id]:
            variant['mole_conseq'] = data['result'][ncbi_id]['molecular_consequence_list']
        else:
            variant['mole_conseq'] = "Unknown"
    return variant_info
        
def classify(variant_info, priority_levels): # classifies variants based on info from annotate() to sort them by priority level
    # mapping of significance and molecular consequence to numerical values --> lower = more severe
    clin_dict = {"pathogenic":0, "likely pathogenic":1, "conflicting classifications of pathogenicity":2, "uncertain significance":3, "likely benign":4, "benign":5, "drug response":6, "not provided":7, "unknown":8}
    mole_dict = {"nonsense variant":0, "nonsense":0, "frameshift variant":1, "splice donor variant":2, "stop lost":3, "missense variant":4, "inframe deletion":5, "in frame deletion":5, "inframe insertion":5, "in frame insertion":5, "5 prime UTR variant":6, "3 prime UTR variant":7, "non-coding transcript variant":7, "intron variant":8, "synonymous variant":9, "no sequence alteration":10}

    for variant in variant_info:
        clin_sig = variant.get('clin_sig', '')
        clin_sig = name_standardize(clin_sig)  # standardize clinical significance
        variant['clin_sig'] = clin_sig  # update the variant dictionary with the standardized clinical significance
        clin_sig_num = clin_dict.get(clin_sig, 9)  # default to 9 if not found

        priority = "Unknown"
        variant['priority_level'] = priority
        for level in priority_levels: # associates clinical significance with priority level based on config.yaml
            if clin_sig in priority_levels[level]:
                priority = level
                variant['priority_level'] = priority
                break

        mole_conseq = variant.get('mole_conseq', '')
        mole_conseq_list = [] # init list to store molecular consequence values --> could be multiple for 1 variant
        for item in mole_conseq:
            item = name_standardize(item)  # standardize molecular consequence
            if item in mole_dict:
                mole_conseq_list.append(mole_dict[item])
        mole_conseq_num = min(mole_conseq_list) if mole_conseq_list else 11  # selects lowest number, default to 11 if not found
        class_tuple = (clin_sig_num, mole_conseq_num)
        variant['classification'] = class_tuple

    sorted_variant_info = sorted(variant_info, key=lambda x: x['classification']) # anon function to sort variants by classification tuple
    return sorted_variant_info

def report(sorted_variant_info, report_file): # generates a report of the sorted variants in a tab-separated format
    report_df = pd.DataFrame(sorted_variant_info, columns = ["chrom", "pos", "ref", "alt", "clin_sig", "mole_conseq", "priority_level"])
    report_df['alt'] = report_df['alt'].apply(lambda x: ','.join(str(i) for i in x)) # converts list of alt alleles to a string for report
    report_df['mole_conseq'] = report_df['mole_conseq'].apply(lambda x: ','.join(str(i) for i in x)) # converts list of molecular consequences to a string for report
    with open(report_file, 'w') as f:
        report_tsv = report_df.to_csv(f, sep="\t",index=False)



if __name__ == "__main__":
    with open('config.yaml', 'r') as f:
        config = yaml.safe_load(f) # loads the config file and stores it in a variable called config as a dictionary

    # stores config info in variables
    vcf_path = config['input']['vcf_file']
    priority_levels = config['priority_levels']
    report_path = config['report_file']

    # call the functions in order
    variant_info = parse(vcf_path)
    variant_info = annotate(variant_info)
    sorted_variant_info = classify(variant_info, priority_levels)
    report(sorted_variant_info, report_path)



