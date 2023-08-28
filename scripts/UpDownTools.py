import argparse
import hashlib
import json
import re
import pandas as pd
import numpy as np
from collections import Counter
from contextlib import ExitStack

MIN_ORF_LENGTH = 300
FASTA_LINE_LENGTH = 60

def read_fasta(file_path, max_sequences=None):
    with open(file_path, 'r') as file:
        identifier = None
        sequence = []
        count = 0
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if identifier is not None:
                    yield identifier, ''.join(sequence)
                    count += 1
                    if max_sequences is not None and count >= max_sequences:
                        break
                identifier = line
                sequence = []
            else:
                sequence.append(line)
        if identifier is not None and (max_sequences is None or count < max_sequences):
            yield identifier, ''.join(sequence)
            count += 1

def remove_duplicate_sequences(filepath, count, fasta_out, json_out):
    hash_dict = {}
    duplicate_count = 0
    with open(fasta_out, 'w') as file:
        for identifier, sequence in read_fasta(filepath, count):
            hash = hashlib.sha256(sequence.encode()).hexdigest()
            if hash not in hash_dict:
                print(identifier, file=file)
                formatted_sequence = '\n'.join([sequence[i:i+FASTA_LINE_LENGTH] for i in range(0, len(sequence), FASTA_LINE_LENGTH)])
                print(formatted_sequence, file=file)
                hash_dict[hash] = [identifier]
            else:
                hash_dict[hash].append(identifier)
                duplicate_count += 1
    
    hash_dict = {v[0]: v[1:] for k, v in hash_dict.items()}
    with open(json_out, 'w') as file:
        json.dump(hash_dict, file)
    print(f"Removed {duplicate_count} duplicate sequences.")

def extract_orfs(fasta_in, count, fasta_out, prefix):
    with open(fasta_out, 'w') as file:
        for identifier, sequence in read_fasta(fasta_in, count):
            ss_orfs = []
            stop_positions = find_substring_indices(sequence)
            stop_positions_by_frame = split_by_modulo(stop_positions)
            stops = match_stop_positions(stop_positions_by_frame)
            ss_orfs.extend([sequence[x+3:y] for x,y in stops if y-x > MIN_ORF_LENGTH])
            ss_orfs_filtered = filter_amiguity(ss_orfs)
            
            for i, s in enumerate(ss_orfs_filtered):
                print(f'>{prefix}{identifier[1:]}_{i}', file=file)
                formatted_sequence = '\n'.join([s[i:i+FASTA_LINE_LENGTH] for i in range(0, len(s), FASTA_LINE_LENGTH)])
                print(formatted_sequence, file=file)

def find_substring_indices(seq):
    substrings = "TAA|TAG|TGA"
    
    matches = re.finditer(substrings, seq)
    indices = [match.start() for match in matches]
    
    return sorted(indices)

def split_by_modulo(lst):
    outlist = [[] for i in range(3)]
    for i in lst:
        outlist[i % 3].append(i)
    return outlist

def match_stop_positions(lst):
    results = []
    for i in lst:
        results.extend(list(zip(i, i[1:])))
    return results

def filter_amiguity(sequences):
    filtered = []
    for seq in sequences:
        non_atgc_count = sum(1 for char in seq if char not in 'ATGC')
        if non_atgc_count / len(seq) <= 0.3:
            filtered.append(seq)
    return filtered

def collect_homologs(fasta_in, blast_in, count, fasta_out):
    df = pd.read_csv(blast_in, sep="\t", names=["sseqid", "qlen", "length", "pident"])
    df["len_ident"] = df["length"] / df["qlen"] * df["pident"]
    df.sort_values("len_ident", inplace=True, ascending=False, ignore_index=True)
    df = df[df['len_ident'] >= 70]
    seq_ids = df['sseqid'].to_list()
    
    with open(fasta_out, 'a') as file:
        for identifier, sequence in read_fasta(fasta_in, count):
            if identifier[1:] in seq_ids:
                if all(base in 'ATGC' for base in sequence):
                    print(identifier.rsplit('_', 1)[0], file=file)
                    formatted_sequence = '\n'.join([sequence[i:i+FASTA_LINE_LENGTH] for i in range(0, len(sequence), FASTA_LINE_LENGTH)])
                    print(formatted_sequence, file=file)

def ava_homologs(fasta_in, blast_in, count, fasta_out):
    df = pd.read_csv(blast_in, sep="\t", names=["qseqid", "sseqid", "qlen", "length", "pident"])
    slen_dict = {(row['qseqid'], row['sseqid']): row['qlen'] for index, row in df.iterrows()}
    df['slen'] = [slen_dict.get((row['sseqid'], row['qseqid']), None) for index, row in df.iterrows()]
    x = np.maximum(df['qlen'], df['slen'])
    df['len_ident'] = (df['length'] / x) * df['pident']
    df_filtered = df[df['len_ident'] >= 70].copy()
    result = df_filtered.groupby('qseqid')['sseqid'].apply(list).reset_index(name='sseqid_list')
    result['sseqid_list'] = result['sseqid_list'].apply(lambda x: sorted(x))
    result_list_of_lists = result['sseqid_list'].to_list()
    unique_results = []
    for r in result_list_of_lists:
        if r not in unique_results:
            unique_results.append(r)
    
    filtered_results = []
    for i in range(len(unique_results)):
        is_subset = False
        for j in range(len(unique_results)):
            if i != j and set(unique_results[i]).issubset(set(unique_results[j])):
                is_subset = True
                break
        if not is_subset:
            filtered_results.append(unique_results[i])

    for i, l in enumerate(unique_results):
        with open(f"{fasta_out}_{i}.fasta", 'w') as file:
            for identifier, sequence in read_fasta(fasta_in, count):
                if identifier[1:] in l:
                    # if all(base in 'ATGC' for base in sequence):
                    # print(identifier.rsplit('_', 1)[0], file=file)
                    print(identifier, file=file)
                    formatted_sequence = '\n'.join([sequence[i:i+FASTA_LINE_LENGTH] for i in range(0, len(sequence), FASTA_LINE_LENGTH)])
                    print(formatted_sequence, file=file)

def split_fasta(fasta_in, num_files, count, fasta_out):
    with ExitStack() as stack:
        files = [stack.enter_context(open(f"{fasta_out}_{i+1}.fasta", 'w')) for i in range(num_files)]
        
        for seq_num, (identifier, sequence) in enumerate(read_fasta(fasta_in, count)):
            target_file = files[seq_num % num_files]
            print(identifier, file=target_file)
            
            formatted_sequence = '\n'.join([sequence[j:j+FASTA_LINE_LENGTH] for j in range(0, len(sequence), FASTA_LINE_LENGTH)])
            print(formatted_sequence, file=target_file)

def fasta_stats(fasta_in, count):
    for identifier, sequence in read_fasta(fasta_in, count):
        print(identifier, sum(1 for char in sequence if char not in ['A', 'T', 'G', 'C']))

def restore_duplicate_sequences(fasta_in, count, json_in, fasta_out):
    with open(json_in) as file:
        duplicates = json.load(file)
    
    with open(fasta_out, 'w') as file:
        for identifier, sequence in read_fasta(fasta_in, count):
            for duplicate in duplicates[identifier]:
                print(duplicate, file=file)
                formatted_sequence = '\n'.join([sequence[j:j+FASTA_LINE_LENGTH] for j in range(0, len(sequence), FASTA_LINE_LENGTH)])
                print(formatted_sequence, file=file)

def remove_gaps(fasta_in, count, fasta_out):
    with open(fasta_out, 'w') as file:
        for identifier, sequence in read_fasta(fasta_in, count):
            sequence = sequence.replace('-', '').upper()
            print(identifier, file=file)
            formatted_sequence = '\n'.join([sequence[j:j+FASTA_LINE_LENGTH] for j in range(0, len(sequence), FASTA_LINE_LENGTH)])
            print(formatted_sequence, file=file)

def parse_args():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='command')
    
    remove_duplicate_sequences_parser = subparsers.add_parser('remove_duplicate_sequences')
    remove_duplicate_sequences_parser.add_argument('-i', '--input')
    remove_duplicate_sequences_parser.add_argument('-c', '--count', type=int, default=None)
    remove_duplicate_sequences_parser.add_argument('-f', '--fasta')
    remove_duplicate_sequences_parser.add_argument('-j', '--json')

    extract_orfs_parser = subparsers.add_parser('extract_orfs')
    extract_orfs_parser.add_argument('-i', '--input')
    extract_orfs_parser.add_argument('-c', '--count', type=int, default=None)
    extract_orfs_parser.add_argument('-o', '--output')
    extract_orfs_parser.add_argument('-p', '--prefix', default=None)

    collect_homologs_parser = subparsers.add_parser('collect_homologs')
    collect_homologs_parser.add_argument('-f', '--fasta')
    collect_homologs_parser.add_argument('-b', '--blast')
    collect_homologs_parser.add_argument('-c', '--count', type=int, default=None)
    collect_homologs_parser.add_argument('-o', '--output')

    ava_homologs_parser = subparsers.add_parser('ava_homologs')
    ava_homologs_parser.add_argument('-f', '--fasta')
    ava_homologs_parser.add_argument('-b', '--blast')
    ava_homologs_parser.add_argument('-c', '--count', type=int, default=None)
    ava_homologs_parser.add_argument('-o', '--output')

    split_fasta_parser = subparsers.add_parser('split_fasta')
    split_fasta_parser.add_argument('-i', '--input')
    split_fasta_parser.add_argument('-t', '--threads', type=int)
    split_fasta_parser.add_argument('-c', '--count', type=int, default=None)
    split_fasta_parser.add_argument('-o', '--output')

    fasta_stats_parser = subparsers.add_parser('fasta_stats')
    fasta_stats_parser.add_argument('-i', '--input')
    fasta_stats_parser.add_argument('-c', '--count', type=int, default=None)

    restore_duplicate_sequences_parser = subparsers.add_parser('restore_duplicate_sequences')
    restore_duplicate_sequences_parser.add_argument('-i', '--input')
    restore_duplicate_sequences_parser.add_argument('-c', '--count', type=int, default=None)
    restore_duplicate_sequences_parser.add_argument('-j', '--json')
    restore_duplicate_sequences_parser.add_argument('-o', '--output')

    remove_gaps_parser = subparsers.add_parser('remove_gaps')
    remove_gaps_parser.add_argument('-i', '--input')
    remove_gaps_parser.add_argument('-c', '--count', type=int, default=None)
    remove_gaps_parser.add_argument('-o', '--output')

    return parser.parse_args()

def main():
    args = parse_args()
    if args.command == 'remove_duplicate_sequences':
        remove_duplicate_sequences(args.input, args.count, args.fasta, args.json)
    elif args.command == 'extract_orfs':
         extract_orfs(args.input, args.count, args.output, args.prefix)
    elif args.command == 'collect_homologs':
         collect_homologs(args.fasta, args.blast, args.count, args.output)
    elif args.command == 'ava_homologs':
         ava_homologs(args.fasta, args.blast, args.count, args.output)
    elif args.command == 'split_fasta':
         split_fasta(args.input, args.threads, args.count, args.output)
    elif args.command == 'fasta_stats':
         fasta_stats(args.input, args.count)
    elif args.command == 'restore_duplicate_sequences':
         restore_duplicate_sequences(args.input, args.count, args.json, args.output)
    elif args.command == 'remove_gaps':
         remove_gaps(args.input, args.count, args.output)

if __name__ == '__main__':
    main()
