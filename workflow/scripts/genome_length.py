from Bio import SeqIO

def calculate_average_sequence_length(fasta_file, output_file):
    total_length = 0
    sequence_count = 0

    # Parse the fasta file and accumulate sequence lengths
    for record in SeqIO.parse(fasta_file, "fasta"):
        total_length += len(record.seq)
        sequence_count += 1

    # Calculate average length if there are sequences
    if sequence_count == 0:
        print("No sequences found in the file.")
        return None

    average_length = int(round(total_length / sequence_count))

    # Write the average length to the output file
    with open(output_file, "w") as f:
        f.write(f"{average_length}")

    return average_length

# Example usage
fasta_file = snakemake.input[0]
output_file = snakemake.output[0]
calculate_average_sequence_length(fasta_file, output_file)
