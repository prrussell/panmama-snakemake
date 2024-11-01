import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Load your data
input_file_path = 'results/initial_runs/run18/results/results_boxplot.tsv'
data = pd.read_csv(input_file_path, sep='\t')

# Define the column names for flexibility
haplotype_column_name = 'Haplotype_Index'
actual_abundance_column_name = 'Abundance_actual'
estimated_abundance_column_name = 'Abundance_est'
grouping_column_name = 'Reads_or_Depth'  # Flexible column, change grouping as needed

# Check that these columns exist in your DataFrame
assert haplotype_column_name in data.columns, f"{haplotype_column_name} not found in data columns."
assert actual_abundance_column_name in data.columns, f"{actual_abundance_column_name} not found in data columns."
assert estimated_abundance_column_name in data.columns, f"{estimated_abundance_column_name} not found in data columns."
assert grouping_column_name in data.columns, f"{grouping_column_name} not found in data columns."

# Check if there is more than one unique level in the grouping column to avoid ZeroDivisionError
n_hue_levels = data[grouping_column_name].nunique()
if n_hue_levels > 1:
    # Set up the plot
    plt.figure(figsize=(12, 8))

    # Overlay with box plots for estimated abundances
    ax = sns.boxplot(
        x=data[haplotype_column_name],
        y=data[estimated_abundance_column_name],
        hue=data[grouping_column_name],
        palette='Pastel1',
        dodge=True
    )
    for artist in ax.lines:
        if artist.get_linestyle() == "None":
            pos = artist.get_xdata()
            artist.set_xdata(pos + np.random.uniform(-.05, .05, len(pos)))

    # Get the x-tick positions for each haplotype to dynamically place lines
    xticks = ax.get_xticks()
    unique_haplotypes = data[haplotype_column_name].unique()
    
    # Map each x-tick to the haplotype so we can align the lines
    for i, haplotype in enumerate(unique_haplotypes):
        # Get the actual abundance value for this haplotype
        actual_abundance = data.loc[data[haplotype_column_name] == haplotype, actual_abundance_column_name].iloc[0]

        # Get the x-position for this haplotype
        x_position = xticks[i]

        # Draw a dotted line across the length of the x-tick for each haplotype
        ax.hlines(
            y=actual_abundance,
            xmin=x_position - 0.4,  # Adjust this value to control span within box plot dodge
            xmax=x_position + 0.4,  # Adjust this value to control span within box plot dodge
            colors='red',
            linestyles='dotted',
            linewidth=2
        )

    # Add titles and labels
    plt.xlabel('Haplotypes')
    plt.ylabel('Abundance')
    plt.title('Estimated vs Actual Abundance by Haplotype')
    plt.legend(title=grouping_column_name)

    plt.tight_layout()
    plt.savefig('test_plot.png', format='png', dpi=300)  # Higher dpi for better resolution
else:
    print("Insufficient unique levels in the grouping column to create a plot with hue.")
