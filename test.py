import argparse
import pandas as pd
import os

def process_tsv(args):
    annotation_data = pd.read_csv("matrix_merged.csv")
    annotation_dict = dict(zip(annotation_data['Gene'], annotation_data['Annotation']))
    
    input_data = pd.read_csv(args.input_tsv, sep='\t', header=None, names=['filepath', 'no_sequences'])
    
    results = []
    
    for _, row in input_data.iterrows():
        filepath = row['filepath']
        df = pd.read_csv(filepath, sep='\t')
        count_positive = len(df[(df['N_dev'] > 0)])
        mean_n_dev_pos = df[df['N_dev'] > 0]['N_dev'].mean()
        stdev_n_dev_pos = df[df['N_dev'] > 0]['N_dev'].std()
        max_n_dev = df['N_dev'].max()  # Compute the maximum value of N_dev column
        aln_len = len(df)
        pos_per_codon = count_positive / aln_len
        base_name = os.path.basename(filepath)
        base_name_without_ext = os.path.splitext(base_name)[0]
        annotation = annotation_dict.get(base_name_without_ext, "N/A")
        
        results.append({
            "gene": base_name_without_ext,
            "annotation": annotation,
            "no_sequences": row['no_sequences'],
            "aln_len": aln_len,
            "pos_count": count_positive,
            "pos_per_codon": round(pos_per_codon, 3),
            "mean_n_dev_pos": round(mean_n_dev_pos, 3),
            "stdev_n_dev_pos": round(stdev_n_dev_pos, 3),
            "max_n_dev": round(max_n_dev, 3)  # Add the max_n_dev to the results
        })
    
    result_df = pd.DataFrame(results)
    result_df = result_df.sort_values(by="max_n_dev", ascending=False)
    
    return result_df

def main():
    parser = argparse.ArgumentParser(description="Process a TSV with file paths and integers.")
    parser.add_argument("input_tsv", help="Path to the input TSV file.")
    args = parser.parse_args()

    result_df = process_tsv(args)
    result_df.to_csv('ranking_table.tsv', sep='\t', index=False)

if __name__ == "__main__":
    main()
