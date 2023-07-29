import argparse
import sqlite3
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_args():
    # Create the main parser
    parser = argparse.ArgumentParser(description='Database operations.')

    # Add the subparsers
    subparsers = parser.add_subparsers(dest='command')

    # Create the parser for the "build_metadata_database" command
    parser_build = subparsers.add_parser('build_metadata_database', help='Build metadata database.')
    parser_build.add_argument('-i', '--input', required=True, help='Input file.')
    parser_build.add_argument('-o', '--output', required=True, help='Output SQLite database file.')

    # Create the parser for the "add_alignment_table" command
    parser_add = subparsers.add_parser('add_alignment_table', help='Add alignment table.')
    parser_add.add_argument('-i', '--input', required=True, help='Input SQLite database file.')
    parser_add.add_argument('-a', '--alignment', required=True, help='Alignment file in Multi-FASTA format.')
    parser_add.add_argument('-t', '--table_name', required=True)

    # Create the parser for the "filter_alignment_table" command
    parser_filter = subparsers.add_parser('filter_alignment_table', help='Filter alignment table.')
    parser_filter.add_argument('-i', '--input', required=True, help='Input SQLite database file.')
    parser_filter.add_argument('-c', '--condition', required=True, help='SQL condition to filter the alignment table.')
    parser_filter.add_argument('-t', '--table_name', required=True)
    parser_filter.add_argument('-o', '--output', required=True, help='Output file in Multi-FASTA format.')

    # Parse the command-line arguments
    return parser.parse_args()

def build_metadata_database(input_file, output_file):
    # Read the input file into a pandas DataFrame
    df = pd.read_csv(input_file, sep='\t')

    # Connect to the SQLite database
    conn = sqlite3.connect(output_file)

    # Write the data to the SQLite database
    df.to_sql('metadata', conn, if_exists='replace', index=False)

    # Close the connection
    conn.close()

def add_alignment_table(database_file, alignment_file, table_name):
    # Read the alignment file into a list of SeqRecord objects
    records = list(SeqIO.parse(alignment_file, 'fasta'))

    # Create a pandas DataFrame from the SeqRecord objects
    df = pd.DataFrame([(record.id, str(record.seq)) for record in records], columns=['id', 'sequence'])

    # Connect to the SQLite database
    conn = sqlite3.connect(database_file)

    # Write the data to the SQLite database
    df.to_sql(table_name, conn, if_exists='replace', index=False)

    # Close the connection
    conn.close()

def filter_alignment_table(database_file, condition, table_name, output_file):
    # Connect to the SQLite database
    conn = sqlite3.connect(database_file)

    # Query the alignment table and the metadata table with the provided condition
    df = pd.read_sql(f'SELECT {table_name}.* FROM {table_name} JOIN metadata ON {table_name}.id = metadata.covv_accession_id WHERE {condition}', conn)

    # Close the connection
    conn.close()

    with open(output_file, 'w') as file:
        for index, row in df.iterrows():
            print(f'>{row["id"]}', file=file)
            formatted_sequence = '\n'.join([row['sequence'][i:i+60] for i in range(0, len(row['sequence']), 60)])
            print(formatted_sequence, file=file)        


if __name__ == "__main__":
    # Parse the command-line arguments
    args = parse_args()

    # Call the appropriate function based on the provided command
    if args.command == 'build_metadata_database':
        build_metadata_database(args.input, args.output)
    elif args.command == 'add_alignment_table':
        add_alignment_table(args.input, args.alignment, args.table_name)
    elif args.command == 'filter_alignment_table':
        filter_alignment_table(args.input, args.condition, args.table_name, args.output)
