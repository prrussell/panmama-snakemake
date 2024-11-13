import pandas as pd
import re
from compute_distance import newick_distances

# List of input/output files passed by Snakemake
input_files = snakemake.input
input_newick_files = {
    re.match(r"([^/]+)\.tree$", file.split("/")[-1]).group(1): file
    for file in input_files if file.endswith(".tree")
}
actual_abundance_file = snakemake.output[0]
estimated_abundance_file = snakemake.output[1]
joined_output_file = snakemake.output[2]
plot_output_file = snakemake.output[3]

# Input file column headers
column_headers = ["Haplotype", "Abundance"]

# Initialize lists to hold DataFrames for each category
panmap_data_frames = []
panman_data_frames = []

# Updated regular expressions for parsing filenames
panmap_pattern = re.compile(
    r"panmap_outputs/(?P<tree>[^_]+)_(?P<num_hap>\d+hap-[a-z])_(?P<simulator>[^_]+)_(?P<depth>\d+)_panmap(?P<param>\d+)_abundance_rep(?P<replicate>\d+)\.tsv"
)
panman_pattern = re.compile(
    r"panman_outputs/(?P<tree>[^_]+)_(?P<num_hap>\d+hap-[a-z])_abundance_rep(?P<replicate>\d+)\.tsv"
)

# Separate and load files based on their prefixes
for file in input_files:
    if "panmap" in file:
        # Extract information using the panmap pattern
        match = panmap_pattern.search(file)
        if match:
            # Load the DataFrame and add additional columns
            df = pd.read_csv(file, sep='\t', names=column_headers)
            df["Tree"] = match.group("tree")
            df["Num_Haplotypes_Param"] = match.group("num_hap")  # Capture entire hap segment (e.g., 2hap-a)
            df["Num_Haplotypes"] = int(re.match(r"\d+", match.group("num_hap")).group())  # Capture only the numeric part
            df["Read_Simulator"] = match.group("simulator")
            df["Reads_or_Depth"] = int(match.group("depth"))
            df["Panmap_Params"] = int(match.group("param"))
            df["Replicate"] = int(match.group("replicate"))
            df["Num_Haplotypes_est"] = len(df)  # Count rows for estimated haplotype count
            panmap_data_frames.append(df)
    elif "panman" in file:
        # Extract information using the panman pattern
        match = panman_pattern.search(file)
        if match:
            # Load the DataFrame and add additional columns
            df = pd.read_csv(file, sep='\t', names=column_headers)
            df["Tree"] = match.group("tree")
            df["Num_Haplotypes_Param"] = match.group("num_hap")  # Capture entire hap segment (e.g., 2hap-a)
            df["Num_Haplotypes"] = int(re.match(r"\d+", match.group("num_hap")).group())  # Capture only the numeric part
            df["Replicate"] = int(match.group("replicate"))
            df["Haplotype_Index"] = df.index + 1
            panman_data_frames.append(df)

# Concatenate all data into two separate DataFrames
panmap_summary_df = pd.concat(panmap_data_frames, ignore_index=True)
panman_summary_df = pd.concat(panman_data_frames, ignore_index=True)

panman_summary_df.to_csv(actual_abundance_file, sep='\t', index=False)
panmap_summary_df.to_csv(estimated_abundance_file, sep='\t', index=False)

# Extract unique combinations of the specified columns
params_df = panmap_summary_df[
    ["Tree", "Num_Haplotypes_Param", "Num_Haplotypes", "Read_Simulator", "Reads_or_Depth", "Panmap_Params", "Replicate"]
].drop_duplicates().reset_index(drop=True)

# Expand true haplotype df to include all pipeline parameter conditions
panman_expanded_df = pd.merge(
    params_df,
    panman_summary_df,
    how="outer"
)

# Step to expand Haplotype column on commas and preserve original list
panmap_summary_df["Haplotypes_All"] = panmap_summary_df["Haplotype"]  # Preserve original list

# Expand rows by splitting "Haplotype" on commas
panmap_summary_df = panmap_summary_df.assign(Haplotype=panmap_summary_df["Haplotype"].str.split(',')).explode("Haplotype").reset_index(drop=True)

# Combine panmap and panman DataFrames by row for a summary
merged_df = panmap_summary_df.merge(
    panman_expanded_df,
    on=["Haplotype", "Tree", "Num_Haplotypes", "Num_Haplotypes_Param", "Replicate", "Read_Simulator", "Reads_or_Depth", "Panmap_Params"],
    how="outer",
    suffixes=('_est', '_actual')
)

# Fill nulls in Haplotypes_All column
merged_df["Haplotypes_All"] = merged_df["Haplotypes_All"].fillna(merged_df["Haplotype"])

# Sort merged_df by 'Abundance_actual' in descending order
merged_df = merged_df.sort_values(by="Abundance_actual", ascending=False)

# Drop duplicates based on the specified columns, keeping only the first occurrence
merged_df = merged_df.drop_duplicates(
    subset=[
        "Haplotypes_All", "Abundance_est", "Tree", "Num_Haplotypes", "Num_Haplotypes_Param", 
        "Read_Simulator", "Reads_or_Depth", "Panmap_Params", "Replicate"],
    keep="first"
)

# Reorder columns in the specified order
merged_df = merged_df[[
    "Tree", "Read_Simulator", "Reads_or_Depth", "Panmap_Params", 
    "Replicate", "Num_Haplotypes_Param", "Num_Haplotypes", "Num_Haplotypes_est",
    "Haplotype_Index", "Haplotype", "Abundance_actual", "Abundance_est", "Haplotypes_All"
]]

# Sort the DataFrame by the specified columns and order
merged_df = merged_df.sort_values(
    by=["Tree", "Num_Haplotypes_Param", "Read_Simulator", "Reads_or_Depth", "Panmap_Params", "Replicate", "Haplotype_Index"]
)

# Fill in missing values for Num_Haplotypes_est
merged_df["Num_Haplotypes_est"] = merged_df.groupby(
    ["Tree", "Read_Simulator", "Reads_or_Depth", "Panmap_Params", "Replicate", "Num_Haplotypes_Param"]
)["Num_Haplotypes_est"].transform("first")

# Fill nulls in Abundance_est column
merged_df["Abundance_est"] = merged_df["Abundance_est"].fillna(0)
merged_df['Haplotype_Index'] = merged_df['Haplotype_Index'].astype('Int64')

# Use newick files to calculate distance 
# Initialize a list to store processed DataFrames for each tree
tree_dfs = []

# Loop over each unique tree in merged_df
for tree in merged_df["Tree"].unique():
    # Extract the subset of merged_df for the current tree
    tree_df = merged_df[merged_df["Tree"] == tree].copy()
    
    # Get the corresponding Newick file for the current tree
    tree_newick_file = input_newick_files.get(tree)
    if tree_newick_file:
        # Apply the newick_distances function to this subset
        tree_df = newick_distances(tree_df, tree_newick_file)
    
    # Append the processed DataFrame to the list
    tree_dfs.append(tree_df)

# Concatenate all the processed DataFrames back into one DataFrame
merged_df = pd.concat(tree_dfs, ignore_index=True)

# Save the combined summary DataFrame to the output file
merged_df.to_csv(joined_output_file, sep='\t', index=False)


## Create dataframe for box plots
plot_df = merged_df[["Tree", "Read_Simulator", "Reads_or_Depth", "Panmap_Params", "Replicate", "Num_Haplotypes_Param",
                     "Num_Haplotypes", "Haplotype_Index", "Abundance_actual", "Abundance_est"]]

# Create "other haplotype" rows for grouping
params_df["Haplotype_Index"] = "others"
params_df["Abundance_est"] = 0
plot_df = pd.concat([plot_df, params_df], ignore_index=True)

# Fill nulls in Abundance_est column
plot_df["Abundance_actual"] = plot_df["Abundance_actual"].fillna(0)

# Fill nulls in Haplotype_Index column
plot_df["Haplotype_Index"] = plot_df["Haplotype_Index"].fillna('others')

# Group data by replicates
plot_df = plot_df.groupby(["Tree", "Read_Simulator", "Reads_or_Depth", "Panmap_Params", "Replicate", "Num_Haplotypes_Param",
                     "Num_Haplotypes", "Haplotype_Index"]).agg({
    'Abundance_actual': 'sum',
    'Abundance_est': 'sum'
}).reset_index()

plot_df.to_csv(plot_output_file, sep='\t', index=False)
