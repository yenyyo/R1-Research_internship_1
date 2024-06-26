import matplotlib as plt
import time
import pandas as pd
import os
from concurrent.futures import ThreadPoolExecutor
import tqdm

# Function to create a pie plot
def pieplot(field):
    # Plot the pie chart
    plt.figure(figsize=(6, 6))
    plt.pie(field, labels=field.index, autopct='%1.1f%%', colors=['skyblue', 'lightcoral'])
    plt.title('Distribution of Sex')
    return plt

# Function to create a bar plot
def barplot(series):
    # Create a bar chart
    plt.bar(series.index, series.values)

    # Set axis labels and title
    plt.xlabel("Racial Group")
    plt.ylabel("Count")
    plt.title("Racial Group Distribution")

    # Rotate x-axis labels for better readability with long labels
    plt.xticks(rotation=45)

    # Display the plot
    plt.show()

# Function to count occurrences of TCRs (T-cell receptors) from a list of files
def count_list(list_files):
    # Initialize an empty dictionary to store counts
    counts = {}

    # Loop through the sample names
    for idx, sample_id in enumerate(list_files, start=1):
        # Start a timer, used later
        start_time = time.time()

        # First, fetch the file corresponding to the current sample code
        try:
            # Read the data in chunks
            chunks = pd.read_csv(f'data/processed_files/{sample_id}_clonotypes.tsv', sep='\t', chunksize=1000000) # chunksize = 1,000,000
            # Process each chunk individually
            for chunk in chunks:
                # Clean up some sequences in cdr3 column
                # Remove everything after the asterisk in the "v_resolved" column
                chunk["v_resolved"] = chunk["v_resolved"].str.replace(r"\*.*$", "", regex=True)

                unique_tcrs = (chunk["v_resolved"] + chunk["cdr3_amino_acid"]).value_counts()

                # Update counts dictionary
                for tcr, count in unique_tcrs.items():
                    counts[tcr] = counts.get(tcr, 0) + 1

        except FileNotFoundError:
            print("Sample code:", sample_id, "has no file associated")
            continue  # Move to the next iteration if the file is not found

        # Print general information
        elapsed_time = time.time() - start_time
        progress = idx / len(list_files) * 100
        print(f"Processing {sample_id} completed in {elapsed_time:.2f} seconds. Progress: {progress:.2f}%")

    # Convert counts dictionary to a list of lists
    tcr_list = [[tcr, count] for tcr, count in counts.items()]

    return tcr_list

# Function to match and extract sample names based on TCR sequences
def obtain_sample_name(tcr_df,list_files):
    # Initialize an empty dataframe
    tcr_sample = []
    matched_cdr3_sequences = []

    # Loop through the sample names
    for idx, sample_id in enumerate(list_files, start=1):

        # Start a timer, for logging purpouses
        start_time = time.time()

        # First, fetch the file corresponding to the current sample_id (sample_name)
        try:
            # Read the data in chunks
            file = pd.read_csv(f'data/processed_files/{sample_id}_clonotypes.tsv', sep='\t') # chunksize = 100.000 gives me about 3 "dfs" per file

            # Check if cdr3_amino_acid is in tcr_df
            matching_cdr3 = file[file["cdr3_amino_acid"].isin(tcr_df["cdr3"])]

            # If there are matching cdr3 sequences, append them to the result list
            if not matching_cdr3.empty:
                matched_cdr3_sequences.extend(zip(matching_cdr3["cdr3_amino_acid"], [sample_id] * len(matching_cdr3)))
                df = pd.DataFrame(matched_cdr3_sequences, columns=['cdr3', 'file_id'])
                df['file_id'] = sample_id
                tcr_sample.append(df)


        except FileNotFoundError:
            print("Sample code:", sample_id, "has no file associated")
            continue  # Move to the next iteration if the file is not found

        # Print general information
        elapsed_time = time.time() - start_time
        progress = idx / len(list_files) * 100
        print(f"Processing {sample_id} completed in {elapsed_time:.2f} seconds. Progress: {progress:.2f}%")

def obtain_sample_name_opt(tcr_df, list_files):
    tcr_sample = []  # List to store DataFrames for each sample

    for idx, sample_id in enumerate(list_files, start=1):
        start_time = time.time()

        try:
            # Read the data in chunks
            chunk_iter = pd.read_csv(f'data/processed_files/{sample_id}_clonotypes.tsv', sep='\t', chunksize=1000000)

            # Iterate over chunks
            for chunk in chunk_iter:

                chunk["v_resolved"] = chunk["v_resolved"].str.replace(r"\*.*$", "", regex=True)

                matching_cdr3 = chunk[(chunk["v_resolved"] + chunk["cdr3_amino_acid"]).isin(tcr_df["cdr3"])].copy()
                if not matching_cdr3.empty:
                    # Add a new column with the combined value
                    matching_cdr3["combined"] = matching_cdr3["v_resolved"] + matching_cdr3["cdr3_amino_acid"]
                    matching_cdr3.loc[:, 'file_id'] = sample_id
                    # Append the new column 'combined' along with 'file_id'
                    tcr_sample.append(matching_cdr3[['combined', 'file_id']])

        except FileNotFoundError:
            print("Sample code:", sample_id, "has no file associated")
            continue

        elapsed_time = time.time() - start_time
        progress = idx / len(list_files) * 100
        print(f"Processing {sample_id} completed in {elapsed_time:.2f} seconds. Progress: {progress:.2f}%")

    return tcr_sample

def print_example(sample_id,tcr_df_threshold):

    file = pd.read_csv(f'data/processed_files/{sample_id}_clonotypes.tsv', sep='\t') # chunksize = 100.000 gives me about 3 "dfs" per file
    print(f"original file {sample_id}_clonotypes.tsv")
    display(file[:5])

    print("tcrs that repeat more than threshold value")
    display(tcr_df_threshold["cdr3"][:5])

    matching_cdr3 = file[file["cdr3_amino_acid"].isin(tcr_df_threshold["cdr3"])]
    print("matiching values inbetween the list and original data")
    display(matching_cdr3[:5])

    matched_cdr3_sequences = []
    matched_cdr3_sequences.extend(zip(matching_cdr3["cdr3_amino_acid"], [sample_id] * len(matching_cdr3)))
    matched_cdr3_sequences

    df = pd.DataFrame(matched_cdr3_sequences, columns=['cdr3', sample_id])
    df[sample_id] = True
    print("join operation on both")
    display(df[:5])

def get_count_threshold(df, unique_value_count):
    """
    This function calculates the count threshold to keep the top tcrs by count.

    Args:
        df (pandas.DataFrame): The dataframe containing tcr and count columns.
        unique_value_count (int): The desired number of unique tcr values to keep in the pruned dataframe.

    Returns:
        float: The count threshold for pruning.
    """
    # Sort the dataframe by count in descending order
    sorted_df = df.sort_values(by="count", ascending=False)

    # Get the count threshold by taking the value at the desired unique_value_count index
    count_threshold = sorted_df.iloc[unique_value_count - 1]['count']

    return count_threshold
