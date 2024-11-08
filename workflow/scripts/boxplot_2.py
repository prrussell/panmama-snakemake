import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import os

# Set up argparse to take input and output file paths as command-line arguments
parser = argparse.ArgumentParser(description="Plot estimated vs actual abundance by haplotype.")
parser.add_argument('input_file_path', type=str, help='Path to the input TSV file')
parser.add_argument('output_directory', type=str, help='Directory to save the output plots')
args = parser.parse_args()

# Load your data
data = pd.read_csv(args.input_file_path, sep='\t')

# Define the column names for flexibility
haplotype_column_name = 'Haplotype_Index'
actual_abundance_column_name = 'Abundance_actual'
estimated_abundance_column_name = 'Abundance_est'
grouping_column_name = 'Reads_or_Depth'  # Flexible column, change grouping as needed
loop_column_name = 'Num_Haplotypes_Param'  # Hardcoded column to loop through for separate plots

# Check that these columns exist in your DataFrame
required_columns = [haplotype_column_name, actual_abundance_column_name, estimated_abundance_column_name, grouping_column_name, loop_column_name]
for col in required_columns:
    assert col in data.columns, f"{col} not found in data columns."

# Ensure output directory exists
os.makedirs(args.output_directory, exist_ok=True)

# Loop through each unique value in the specified loop column
for group_value in data[loop_column_name].unique():
    # Filter the data for the current group
    subset_data = data[data[loop_column_name] == group_value]

    # Set up the plot
    plt.figure(figsize=(12, 8))

    # Overlay with box plots for estimated abundances
    ax = sns.boxplot(
        x=subset_data[haplotype_column_name],
        y=subset_data[estimated_abundance_column_name],
        hue=subset_data[grouping_column_name],
        palette='Pastel1',
        dodge=True
    )
    for artist in ax.lines:
        if artist.get_linestyle() == "None":
            x_pos = artist.get_xdata()
            y_pos = artist.get_ydata()

            # Add jitter to both x and y positions
            artist.set_xdata(x_pos + np.random.uniform(-0.01, 0.01, len(x_pos)))
            artist.set_ydata(y_pos + np.random.uniform(-0.01, 0.01, len(y_pos)))

    # Get the x-tick positions for each haplotype to dynamically place lines
    xticks = ax.get_xticks()
    unique_haplotypes = subset_data[haplotype_column_name].unique()

    # Map each x-tick to the haplotype so we can align the lines
    for i, haplotype in enumerate(unique_haplotypes):
        # Get the actual abundance value for this haplotype
        actual_abundance = subset_data.loc[subset_data[haplotype_column_name] == haplotype, actual_abundance_column_name].iloc[0]

        # Get the x-position for this haplotype
        x_position = xticks[i]

        # Draw a dotted line across the length of the x-tick for each haplotype
        ax.hlines(
            y=actual_abundance,
            xmin=x_position - 0.4,  # Adjust this value to control span within box plot dodge
            xmax=x_position + 0.4,  # Adjust this value to control span within box plot dodge
            colors='red',
            linestyles='dotted',
            linewidth=3
        )

    # Add titles and labels
    plt.xlabel('Haplotypes')
    plt.ylabel('Abundance')
    plt.title(f'Estimated vs Actual Abundance by Haplotype (Group: {group_value})')
    plt.legend(title=grouping_column_name)

    # Save each plot with a unique filename in the specified output directory
    output_path = os.path.join(args.output_directory, f"boxplot_{group_value}.png")
    plt.tight_layout()
    plt.savefig(output_path, format='png', dpi=300)
    plt.close()  # Close the plot to avoid display overlap
