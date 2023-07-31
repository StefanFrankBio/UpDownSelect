import argparse
import json
import pandas as pd

def parse_args():
    parser = argparse.ArgumentParser()
    # parser.add_argument("-m", "--meme" )
    parser.add_argument("-f", "--fel" )
    parser.add_argument("-s", "--slac" )
    # parser.add_argument("-c", "--codeml" )
    parser.add_argument("-d", "--dNdS" )
    # parser.add_argument("-n", "--dNdS_with_dups" )
    parser.add_argument("-o", "--output")
    return parser.parse_args()

def main():
    args = parse_args()
    # with open(args.meme) as json_file:
    #     data = json.load(json_file)
    #     content = data['MLE']['content']
    #     meme_p_values = [row[6] for row in content['0']]

    with open(args.fel) as json_file:
        data = json.load(json_file)
        content = data['MLE']['content']
        fel_p_values = [row[4] for row in content['0']]
        fel_alpha = [row[0] for row in content['0']]
        fel_beta = [row[1] for row in content['0']]
        fel_classification = []
        for i,j in zip(fel_alpha, fel_beta):
            if i > j:
                fel_classification.append('Neg')
            elif i < j:
                fel_classification.append('Pos')
            else:
                fel_classification.append('None')
    
    with open(args.slac) as json_file:
        data = json.load(json_file)
        content = data['MLE']['content']['0']['by-site']
        slac_p_values = [(row[-3], row[-2]) for row in content['AVERAGED']]
        slac_classification = []
        slac_p_values_filtered = []
        for i, j in slac_p_values:
            if i < j:
                slac_classification.append('Pos')
                slac_p_values_filtered.append(i)
            elif i > j:
                slac_classification.append('Neg')
                slac_p_values_filtered.append(j)
            else:
                slac_classification.append('None')
                slac_p_values_filtered.append(i)
    
    
    # with open(args.codeml) as file:
    #     beb_flag = False
    #     codeml_probabilities = []
    #     for line in file:
    #         if "Bayes Empirical Bayes" in line:
    #             beb_flag = True
    #         elif beb_flag and line.strip() != '' and line.strip()[0].isdigit():
    #             parts = line.split()
    #             if len(parts) > 11:
    #                 codeml_probabilities.append(parts[12])

    df1 = pd.read_csv(args.dNdS, sep="\t")
    # df2 = pd.read_csv(args.dNdS_with_dups, sep="\t")
    # df1['dNdS_p_values_with_dups'] = df2.iloc[:, -1].tolist()
    # df1['meme_p_value'] = meme_p_values
    df1['fel_p_value'] = fel_p_values
    df1['fel_classification'] = fel_classification
    df1['slac_p_value'] = slac_p_values_filtered
    df1['slac_classification'] = slac_classification
    # df1['codeml_probabilities'] = codeml_probabilities
    df1.to_csv(args.output, sep='\t', index=False)

if __name__ == '__main__':
    main()