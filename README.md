# Notes
- ambiguous bases are excluded during 'collect_homologs'
- sed command in translate function removes frame-info suffix (_1) from ids
- where is the correct step to remove duplicates, beginning and homologs?
- write script to regenerate duplicates
    - where should duplicates be reintroduced?
- multiple substitutions at the same codon lead to weird behavior, investigate
- How can I make ambiguous not count towards indentity? Should be filtered
- Write script to impute msa for amiguity handling. Fasta -> Transpose -> Transpose -> Fasta

# UpDownSelect
## Installation Process
```bash
git clone https://github.com/StefanFrankBio/UpDownSelect.git
micromamba env create -f environments/UpDownSelect.yml
micromamba activate UpDownSelect
```
## Input Data
1. fasta file containing nucleotide sequences to be analyzed
2. reference fasta file containing nucleotide sequences for genes of interest
3. tsv file containing metadata
## How to run the pipeline
```bash
UpDownSelect sequences.fasta references.fasta metadata.tsv output_directory threads
```
## Workflow
### Step1: Extract Orfs
