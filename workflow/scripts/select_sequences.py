import sys
import random
import re
from Bio import SeqIO

def random_fasta_selection(input_fasta, output_fasta, output_fasta_renamed, num_sequences):
    # Read all sequences from the input fasta file
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    
    # Check if there are enough sequences to select from
    if len(sequences) < num_sequences:
        print(f"Error: The input file contains less than {num_sequences} sequences.")
        sys.exit(1)
    
    # Randomly select 'num_sequences' sequences
    selected_sequences = random.sample(sequences, num_sequences)

    # Write the selected sequences to the output fasta file
    with open(output_fasta, "w") as output_handle:
        SeqIO.write(selected_sequences, output_handle, "fasta")
    
    # Replace all non-alphanumeric characters in the sequence IDs and descriptions
    for seq_record in selected_sequences:
        seq_record.id = re.sub(r'\W', '_', seq_record.id)
        seq_record.description = re.sub(r'\W', '_', seq_record.description)

    # Write the selected sequences to the output fasta file
    with open(output_fasta_renamed, "w") as output_handle:
        SeqIO.write(selected_sequences, output_handle, "fasta")
    
    print(f"Successfully wrote {num_sequences} random sequences to {output_fasta}.")

# Main script logic
if __name__ == "__main__":
    # Access input, output, and wildcards from snakemake object
    input_fasta = snakemake.input[0]
    output_fasta = snakemake.output[0]
    output_fasta_renamed = snakemake.output[1]
    num_sequences = snakemake.params.hap_count
    random_fasta_selection(input_fasta, output_fasta, output_fasta_renamed, num_sequences)
