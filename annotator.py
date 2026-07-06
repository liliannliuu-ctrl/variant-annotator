from Bio import Entrez
from numpy import record
import vcf
import requests, yaml, os, json
from dotenv import load_dotenv
import pandas as pd

load_dotenv()  # Load environment variables from .env file
Entrez.email = "liliannliuu@gmail.com" # module-level attributes for verification purposes
Entrez.api_key = os.getenv("NCBI_API_KEY")

def parse(vcf_path):
    # reads vcf file, storing as list of Record objects?...with each row/variant being one record object
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

def annotate(variant_info):
    for variant in variant_info:
        # fetching ncbi ids
        handle = Entrez.esearch(db="clinvar", term=variant["id"], retmax=1)
        record = Entrez.read(handle) 
        handle.close()
        ncbi_id = record["IdList"][0] if record["IdList"] else None
        variant['ncbi_id'] = ncbi_id

        # getting clinvar data using ncbi ids
        handle = Entrez.esummary(db="clinvar", id=variant["ncbi_id"], retmode="json")
        data = json.loads(handle.read().decode())
        handle.close()

        # access & store clinical significance
        if 'germline_classification' in data['result'][ncbi_id]:
            variant['clin_sig'] = data['result'][ncbi_id]['germline_classification']['description']
        else:
            variant['clin_sig'] = "Unknown"
        print(variant['clin_sig'])

        # access & store functional consequence
        if 'molecular_consequence_list' in data['result'][ncbi_id]:
            variant['mole_conseq'] = data['result'][ncbi_id]['molecular_consequence_list']
        else:
            variant['mole_conseq'] = "Unknown"

    # write handle info into a file to see structure
    with open("debug_output.json", "w") as f:
        f.write(json.dumps(data))
    return variant_info
        
def classify(variant_info):
    # mapping of significance and molecular consequence to numerical values
    clin_dict = {"pathogenic":0, "likely pathogenic":1, "uncertain significance":2, "likely benign":3, "benign":4, "drug response":5}
    mole_dict = {"frameshift variant":0, "nonsense variant": 1, "inframe deletion":2, "in frame deletion":2, "inframe insertion":2, "in frame insertion":2,"missense variant":3, "synonymous variant":4}
    for variant in variant_info:
        clin_sig = variant.get('clin_sig', '')
        clin_sig = name_normalize(clin_sig)  # normalize clinical significance
        clin_sig_num = clin_dict.get(clin_sig, 6)  # default to 6 if not found
        print(f"Clinical significance: {clin_sig}, Numerical value: {clin_sig_num}")
        print()
        print()

        mole_conseq = variant.get('mole_conseq', '')
        mole_conseq_list = [] # init list to store molecular consequence values --> could be multiple for 1 variant
        for item in mole_conseq:
            item = name_normalize(item)  # normalize molecular consequence
            if item in mole_dict:
                mole_conseq_list.append(mole_dict[item])
        mole_conseq_num = min(mole_conseq_list) if mole_conseq_list else 5  # default to 5 if not found
        class_tuple = (clin_sig_num, mole_conseq_num)
        print(f"Molecular consequence: {item}, Numerical value: {mole_conseq_num}")
        print()
        print()
        variant['classification'] = class_tuple

    sorted_variant_info = sorted(variant_info, key=lambda x: x['classification']) # anon function to sort by classification tuple
    for variant in sorted_variant_info:
        print(f"Variant: {variant['id']}, Classification: {variant['classification']}")
    return sorted_variant_info

def report(sorted_variant_info, report_file):
    report_df = pd.DataFrame(sorted_variant_info, columns = ['chrom', 'pos', 'ref', 'alt', 'clin_sig', 'mole_conseq'])
    report_df['alt'] = report_df['alt'].apply(lambda x: ','.join(str(i) for i in x))
    report_df['mole_conseq'] = report_df['mole_conseq'].apply(lambda x: ','.join(str(i) for i in x))
    with open(report_file, 'w') as f:
        report_tsv = report_df.to_csv(f, sep="\t",index=False)

def name_normalize(name):
    # normalize names to lowercase and remove spaces
    norm_name = name.lower().replace("_", " ").replace("-", " ").strip()
    return norm_name



if __name__ == "__main__":
    with open('config.yaml', 'r') as f:
        config = yaml.safe_load(f) # loads the config file and stores it in a variable called config as a dictionary

    vcf_path = config['input']['vcf_file']
    # calls parse function
    variant_info = parse(vcf_path)
    variant_info = annotate(variant_info)
    sorted_variant_info = classify(variant_info)
    report_tsv = report(sorted_variant_info, config['output']['report_file'])



