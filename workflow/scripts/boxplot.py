import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import argparse
import os

# Set up argparse to take input and output file paths as command-line arguments
parser = argparse.ArgumentParser(description="Plot estimated vs actual abundance by grouped haplotypes.")
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

# Helper function to sort haplotype indices with "other" last
def custom_sort_key(value):
    try:
        # Try to convert to integer if possible
        return (0, int(value))
    except ValueError:
        # Place "other" after numeric values
        return (1, value)

# Loop through each unique value in the specified loop column
for group_value in data[loop_column_name].unique():
    # Filter the data for the current group
    subset_data = data[data[loop_column_name] == group_value]

    # Group data by actual abundance, sorting by abundance in descending order
    grouped_data = subset_data.groupby(actual_abundance_column_name, sort=False)
    sorted_abundances = sorted(grouped_data, key=lambda x: -x[0])  # Sort by actual abundance descending

    # Prepare data for plotting by assigning each actual abundance to a unique x position
    plot_data = []
    x_labels = []  # Store custom labels for each x-tick
    x_positions = {}  # Track x positions for each unique actual abundance
    for i, (abundance, group) in enumerate(sorted_abundances):
        x_positions[abundance] = i
        group["x_position"] = i
        plot_data.append(group)

        # Sort haplotype indices with custom sort, placing "other" last
        haplotypes = sorted(group[haplotype_column_name].unique(), key=custom_sort_key)
        
        # Create x-tick label based on haplotype indices, now sorted
        if len(haplotypes) > 1:
            x_labels.append(f"{haplotypes[0]}-{haplotypes[-1]}")
        else:
            x_labels.append(str(haplotypes[0]))

    plot_data = pd.concat(plot_data)

    # Set up the plot
    plt.figure(figsize=(12, 8))
    ax = sns.boxplot(
        x=plot_data["x_position"],
        y=plot_data[estimated_abundance_column_name],
        hue=plot_data[grouping_column_name],
        palette='Pastel1',
        dodge=True
    )

    # Apply jitter to both x and y positions for box plots
    for artist in ax.lines:
        if artist.get_linestyle() == "None":
            x_pos = artist.get_xdata()
            y_pos = artist.get_ydata()
            # Add jitter
            artist.set_xdata(x_pos + np.random.uniform(-0.01, 0.01, len(x_pos)))
            artist.set_ydata(y_pos + np.random.uniform(-0.01, 0.01, len(y_pos)))

    # Draw a single red dotted line per unique actual abundance
    for abundance, x_position in x_positions.items():
        ax.hlines(
            y=abundance,
            xmin=x_position - 0.4,
            xmax=x_position + 0.4,
            colors='red',
            linestyles='dotted',
            linewidth=3
        )

    # Set x-tick labels to the list of haplotype ranges or indices
    ax.set_xticks(list(x_positions.values()))
    ax.set_xticklabels(x_labels)

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