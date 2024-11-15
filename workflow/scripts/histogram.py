import pandas as pd
import matplotlib.pyplot as plt

# Save textfile (so Snakemake initiates rule)
text_to_save = "Creating histograms."

# Specify the filename and path
file_path = snakemake.output[0]

# Open the file in write mode and save the text
with open(file_path, 'w') as file:
    file.write(text_to_save)

# Load the data
data = pd.read_csv(snakemake.input[0], sep='\t')

# Get the maximum number of replicates for normalization
max_replicates = data['Replicate'].max()

# Define the parameter columns for grouping
parameter_columns = ['Tree', 'Read_Simulator', 'Reads_or_Depth', 'Panmap_Params', 'Num_Haplotypes_Param']

# Loop through each unique combination of parameters
for parameter_values, group_data in data.groupby(parameter_columns):
    # Calculate normalized sum of Abundance_est for each Newick_Dist
    abundance_by_dist = group_data.groupby('Newick_Dist')['Abundance_est'].sum() / max_replicates

    # Create the histogram plot
    plt.figure(figsize=(10, 6))
    plt.bar(abundance_by_dist.index, abundance_by_dist.values, color='skyblue', edgecolor='black')
    
    # Set x and y labels and title
    plt.xlabel('Newick Distance')
    plt.ylabel('Normalized Abundance Estimate')
    title = " - ".join([f"{col}: {val}" for col, val in zip(parameter_columns, parameter_values)])
    plt.title(f"Histogram of Normalized Abundance by Newick Distance\n{title}")

    # Save the plot with a unique filename based on parameter values
    filename = f"histogram_{'_'.join(map(str, parameter_values))}.png"
    output_path = f"plots/histograms/{filename}"
    plt.tight_layout()
    plt.savefig(output_path, format='png', dpi=300)
    plt.close()  # Close the plot to avoid display overlap
