from Bio import Entrez
from numpy import record
import vcf
import requests
import yaml

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
    clin_sigs = []
    return clin_sigs

def classify(clin_sigs):
    classifications = []
    return classifications

def report(classifications):
    printed_report = []
    return printed_report

if __name__ == "__main__":
    print("This module is intended to be imported and used as a library, not run directly.")
    with open('config.yaml', 'r') as f:
        config = yaml.safe_load(f) # loads the config file and stores it in a variable called config as a dictionary

    vcf_path = config['input']['vcf_file']
    # calls parse function
    variant_info = parse(vcf_path)
    print(variant_info)
