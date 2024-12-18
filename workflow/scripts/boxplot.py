import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load your data
data = pd.read_csv(snakemake.input[0], sep='\t')

# Save textfile (so Snakemake initiates rule)
text_to_save = "Creating boxplots."

# Specify the filename and path
file_path = snakemake.output[0]

# Open the file in write mode and save the text
with open(file_path, 'w') as file:
    file.write(text_to_save)

# Define the column names for flexibility (can be single string or list of strings)
parameter_columns = ['Tree', 'Read_Simulator', 'Reads_or_Depth', 'Panmap_Params', 'Num_Haplotypes_Param']
haplotype_column_name = 'Haplotype_Index'
actual_abundance_column_name = 'Abundance_actual'
estimated_abundance_column_name = 'Abundance_est'
grouping_column_name = ['Reads_or_Depth']  # List of columns to group for boxplots. Determines number of different colored boxplots and legend values.
# Create separate plots for any parameter variables not captured in grouping column
loop_column_name = sorted(list(set(parameter_columns) - set(grouping_column_name)))

# Helper function to handle single vs multiple columns
def get_combined_column(df, columns):
    if isinstance(columns, list):
        combined = df[columns].astype(str).agg('_'.join, axis=1)
    else:
        combined = df[columns]
    return combined

# Create combined columns if needed
data['combined_group'] = get_combined_column(data, grouping_column_name)
data['combined_loop'] = get_combined_column(data, loop_column_name)

# Helper function to sort haplotype indices with "other" last
def custom_sort_key(value):
    try:
        return (0, int(value))
    except ValueError:
        return (1, value)

# Loop through each unique value in the combined loop column
for group_value in data['combined_loop'].unique():
    subset_data = data[data['combined_loop'] == group_value]

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
        hue=plot_data["combined_group"],
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
    plt.xlabel('Haplotypes', fontsize=14)
    plt.ylabel('Abundance', fontsize = 14)
    plt.xticks(fontsize=12) 
    plt.yticks(fontsize=12) 
    # Get the parameter values for the current group
    loop_values = [subset_data[col].iloc[0] for col in loop_column_name]
    title_details = " - ".join([f"{col}: {val}" for col, val in zip(loop_column_name, loop_values)])
    plt.title(f"Estimated vs Actual Abundance by Haplotype\n{title_details}", fontsize = 14)
    plt.legend(title=" & ".join(grouping_column_name) if isinstance(grouping_column_name, list) else grouping_column_name,
               fontsize = 12,
               title_fontsize = 14)

    # Save each plot with a unique filename in the specified output directory
    path_details = "_".join(str(val) for val in loop_values)
    output_path = f"plots/boxplots/boxplot_{path_details}.png"
    plt.tight_layout()
    plt.savefig(output_path, format='png', dpi=300)
    plt.close()  # Close the plot to avoid display overlap