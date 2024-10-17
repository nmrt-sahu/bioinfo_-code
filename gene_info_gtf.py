import pandas as pd 

def id_version_column(base_id, version):
    if version == 'NA':
        return base_id
    else:
        return f"{base_id}.{version}"
    
def transform_gtf_to_df(input_df):
    # Define the keys that will be used as headers in the CSV
    keys = ["gene_id", "gene_version", "transcript_id", "transcript_version", "gene_name", "gene_source",
            "gene_biotype", "transcript_name", "transcript_source",
            "transcript_biotype", "exon_number", "exon_version", "exon_id", "tag"]

    # Initialize an empty list to store the rows of the output DataFrame
    output_rows = []

    for _, row in input_df.iterrows():
        # Split the 9th column by semicolons to get the key-value pairs
        kv_pairs = row[8].split('; ')

        # Create a dictionary to store the key-value pairs, initialize with 'NA'
        kv_dict = {key: 'NA' for key in keys}

        for pair in kv_pairs:
            # Split each pair by space to separate key and value
            key_value = pair.split(' ', 1)
            # Only process if we have both key and value (some pairs may not be well-formed)
            if len(key_value) == 2:
                key, value = key_value
                value = value.strip('"').strip('";')
                # Add to dictionary if key is in the predefined list
                if key in keys:
                    kv_dict[key] = value

        # Combine the values from columns 0 to 7 with the kv_dict
        combined_row = {f"col_{i}": row[i] for i in range(8)}
        combined_row.update(kv_dict)
        
        # Create new columns using the create_combined_column function
        combined_row['gene_id_new'] = id_version_column(combined_row.get('gene_id', 'NA'), combined_row.get('gene_version', 'NA'))
        combined_row['transcript_id_new'] = id_version_column(combined_row.get('transcript_id', 'NA'), combined_row.get('transcript_version', 'NA'))
        combined_row['exon_id_new'] = id_version_column(combined_row.get('exon_id', 'NA'), combined_row.get('exon_version', 'NA'))
        
        # Append the dictionary values as a new row to the output list
        output_rows.append(combined_row)

    # Create a new DataFrame from the output rows
    output_df = pd.DataFrame(output_rows)

    # Drop unnecessary columns
    columns_to_drop = ['gene_id', 'gene_version', 'transcript_id', 'transcript_version', 'exon_id', 'exon_version', 'exon_number']
    output_df = output_df.drop(columns=columns_to_drop, errors='ignore')
    new_column_names = {
        "col_0": "chr",
        "col_1": "source",
        "col_2": "feature",
        "col_3": "start",
        "col_4": "end",
        "col_5": "score",
        "col_6": "strand",
        "col_7": "frame",
        "transcript_id_new": "transcript_id",
        "gene_id_new": "gene_id"
    }
    output_df = output_df.rename(columns=new_column_names)

    return output_df


gtf_df = pd.read_csv("Arabidopsis_thaliana.TAIR10.56.gtf", sep='\t', header=None, comment='#', low_memory=False)

# Filter rows where the 'feature' column value is "transcript"
transcript_df = gtf_df[gtf_df.iloc[:, 2] == "gene"]
print("Extracting transcripts information from the GTF file for the next step...")

# Transform the DataFrame
new_transcript_df = transform_gtf_to_df(transcript_df)

#print(new_transcript_df)


# Assuming df1 and new_transcript_df are the two dataframes
# df1: DataFrame with 'bpt' and 'chr'
# new_transcript_df: DataFrame with 'start', 'end', 'chr', and 'gene_id'

# Load df1 and df2
df1 = pd.read_csv("sakshi_fusion.csv", sep=',')

# Convert the chr1 and chr2 columns in df1 to strings
df1['chr1'] = df1['chr1'].astype(str)
df1['chr2'] = df1['chr2'].astype(str)

# Assuming new_transcript_df is already defined and loaded
df2 = new_transcript_df[['chr', 'gene_id', 'start', 'end']]

# Convert the chr column in df2 to string
df2['chr'] = df2['chr'].astype(str)

# Print df1 to verify its structure
#print("df1:")
#print(df1.head())

#print(df1['chr1'])
# Create empty columns for gene_id, start, and end in df1
df1['gene_1'] = None
df1['start_1'] = None
df1['end_1'] = None
df1['gene_2'] = None
df1['start_2'] = None
df1['end_2'] = None

# Define a function to perform the matching and assign gene_id, start, end to the dataframe
def match_bpt_and_add_info(df1, df2, chr_col, bpt_col, gene_id_col, start_col, end_col):
    for idx, row1 in df1.iterrows():
        # Filter rows in df2 where 'chr' matches
        df2_filtered = df2[df2['chr'] == row1[chr_col]]
        
        # Check if any rows in df2 match the bpt range
        for _, row2 in df2_filtered.iterrows():
            if row1[bpt_col] >= row2['start'] and row1[bpt_col] <= row2['end']:
                # Assign gene_id, start, end from df2 to the corresponding columns in df1
                df1.at[idx, gene_id_col] = row2['gene_id']
                df1.at[idx, start_col] = row2['start']
                df1.at[idx, end_col] = row2['end']
                if row1['chr1']==row1['chr2']:
                    df1.at[idx, 'feature'] = 'Intra'
                else:
                    df1.at[idx, 'feature'] = 'Inter'

                break  # Exit loop after the first match

# Match chr1 and bpt1 from df1 with df2
match_bpt_and_add_info(df1, df2, 'chr1', 'bpt1', 'gene_1', 'start_1', 'end_1')

# Match chr2 and bpt2 from df1 with df2
match_bpt_and_add_info(df1, df2, 'chr2', 'bpt2', 'gene_2', 'start_2', 'end_2')

# Print the updated df1
print("Updated df1:")
print(df1)
df1.to_csv('sakshi_fusion_updated.csv',sep='\t')