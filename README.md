# VARIANT ANNOTATOR
A vcf (variant call format) file is a standardized, tab-separated file commonly used in bioinformatics and genomics to store variations in a DNA sequence (such as insertions or deletions) compared to a reference genome. This variant annotator algorithm parses and analyzes an input vcf file, extracting the known clinical significance and functional consequence associated with that variant from public NCBI databases. Essentially, this algorithm retrieves information on whether the variant is disease-causing and the type of mutation it constitutes. Additionally, this algorithm also classifies the variant's severity level depending on its clinical significance and functional consequence. Finally, this algorithm writes a tsv file containing information about the variant for easy access. 

## TECHNOLOGIES
- Python
- PyVCF3
- Biopython (Bio.Entrez)
- pandas
- PyYAML
- python-dotenv

## SETUP
1. Clone the repo
2. Create and activate a virtual environment
3. Install dependencies (pip install -r requirements.txt)
4. Create a .env file with personal NCBI API key
    ie. NCBI_API_KEY=your_key_here
5. Add a VCF file to the data/directory
6. Update config.yaml with VCF filename

## USAGE
python annotator.py

## OUTPUT
Running this program outputs a file called 'report.tsv', which outputs information about each variant separated by tabs: chromosome, position on the genome, reference and alternate nucleotides, clinical significance (ie. pathogenic), molecular consequence (ie. missense variant), and priority level. Priority level was assigned by this algorithm depending on the retrieved NCBI annotations and flags clinically actionable variants. 
