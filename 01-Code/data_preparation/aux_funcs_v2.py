## I have commented out the methods which are not used in the moment, left them here in case i need them later (or to copy pieces of it)

import matplotlib as plt
import time
import pandas as pd
import os
from concurrent.futures import ThreadPoolExecutor
import tqdm

def pieplot(field):
    # Plot the pie chart
    plt.figure(figsize=(6, 6))
    plt.pie(field, labels=field.index, autopct='%1.1f%%', colors=['skyblue', 'lightcoral'])
    plt.title('Distribution of Sex')
    return plt

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
            chunks = pd.read_csv(f'data/processed_files/{sample_id}_clonotypes.tsv', sep='\t', chunksize=1000000) # chunksize = 100.000 gives me about 3 "dfs" per file
            # Process each chunk individually
            for chunk in chunks:

                # Clean up some sequences in cdr3 column
                # chunk["cdr3_amino_acid"] = chunk[chunk['cdr3_amino_acid'].str.contains('^C[^*_]*[FW]$', na=False)]                # Remove everything after the asterisk in the "v_resolved" column
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

    return tcr_sample

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

def obtain_sample_name_opt2(tcr_df, list_files):
    tcr_sample = []  # List to store DataFrames for each sample

    for idx, sample_id in enumerate(list_files, start=1):
        start_time = time.time()

        try:
            # Read the data in chunks
            chunk_iter = pd.read_csv(f'data/processed_files/{sample_id}_clonotypes.tsv', sep='\t', chunksize=100000)

            for chunk in chunk_iter:
                matching_cdr3_index = chunk["cdr3_amino_acid"].isin(tcr_df["cdr3"])
                matching_cdr3 = chunk.loc[matching_cdr3_index, :]
                if not matching_cdr3.empty:
                    matching_cdr3.loc[:, 'file_id'] = sample_id
                    tcr_sample.append(matching_cdr3[['cdr3_amino_acid', 'file_id']])

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


# class ChunkProcessing:
#     def __init__(self):
#         pass

#     def process_chunks(self, list_files, process_chunk_func):
#             # Initialize an empty dictionary to store counts
#             counts = {}

#             # Loop through the sample names
#             for idx, sample_id in enumerate(list_files, start=1):
#                 # Start a timer, used later
#                 start_time = time.time()

#                 # First, fetch the file corresponding to the current sample code
#                 try:
#                     # Read the data in chunks
#                     chunks = pd.read_csv(f'data/processed_files/{sample_id}_clonotypes.tsv', sep='\t', chunksize=1000000) # chunksize = 100.000 gives me about 3 "dfs" per file
#                     # Process each chunk individually
#                     for chunk in chunks:
#                         process_chunk_func(chunk, counts)
#                 except FileNotFoundError:
#                     print("Sample code:", sample_id, "has no file associated")
#                     continue  # Move to the next iteration if the file is not found

#                 # Print general information
#                 elapsed_time = time.time() - start_time
#                 progress = idx / len(list_files) * 100
#                 print(f"Processing {sample_id} completed in {elapsed_time:.2f} seconds. Progress: {progress:.2f}%")


# def count_list(chunk, counts):
   
#     unique_tcrs = chunk["cdr3_amino_acid"].value_counts()
#     # Update counts dictionary
#     for tcr, count in unique_tcrs.items():
#         counts[tcr] = counts.get(tcr, 0) + 1


# def concat_df(tcr_list, list_files):
    
#     # Concatenate all DataFrames in tcr_list into a single DataFrame
#     tcr_list = pd.concat(tcr_list, ignore_index=True)
#     # Group by TCRs and sum counts across all files
#     tcr_count = tcr_list.groupby("cdr3_amino_acid").sum().reset_index()
#     # How many times does it appear across files?
#     tcr_count['count'] = tcr_count[list_files].sum(axis=1)

#     # Print the resulting DataFrame tcr_counts
#     return tcr_count



# # Define a function to process dataframes in batches
# def process_batches(dataframes, batch_size):
#     for i in range(0, len(dataframes), batch_size):
#         yield pd.concat(dataframes[i:i+batch_size], ignore_index=True)


# def concat_batch_df(tcr_list,list_files,batch_size):
#     # Initialize an empty list to store the results
#     tcr_counts = []
#     file_counter = 0

#     for batch_df in process_batches(tcr_list, batch_size):
#         # Filter out columns that are not present in the current dataframe
#         relevant_columns = [col for col in list_files if col in batch_df.columns]
#         # Group by TCRs and sum counts across relevant columns
#         tcr_count = batch_df.groupby("cdr3_amino_acid")[relevant_columns].sum().reset_index()
#         # How many times does it appear across files?
#         tcr_count['count'] = tcr_count[relevant_columns].sum(axis=1)
#         file_name = f'final_tcr_counts_{file_counter}.pkl'
#         # Save the batch DataFrame to disk with the constructed file name
#         tcr_count.to_pickle(os.path.join('data/batches', file_name))
#         # Increment the file counter
#         file_counter += 1

#         # Append the result to the list
#         tcr_counts.append(tcr_count)

#     # Concatenate the list of results into a single DataFrame
#     final_tcr_counts = pd.concat(tcr_counts, ignore_index=True)

#     # Print the resulting DataFrame
#     return final_tcr_counts
    

# def dataframe_list(list_files):
#     # Initialize an empty list to store dataframes containing our tcr info that have been found
#     tcr_list = []

#     # Loop through the sample names
#     for idx, sample_id in enumerate(list_files, start=1):
#         # Start a timer, used later
#         start_time = time.time()

#         # First, fetch the file corresponding to the current sample code
#         try:
#             # Read the data in chunks
#             chunks = pd.read_csv(f'data/processed_files/{sample_id}_clonotypes.tsv', sep='\t', chunksize=1000000) # chunksize = 100.000 gives me about 3 "dfs" per file
#             # Process each chunk individually
#             for chunk in chunks:
#                 unique_tcrs = chunk["cdr3_amino_acid"].unique()
#                 # Create a DataFrame with TCRs and their corresponding filename where we have found it
#                 tcr_df = pd.DataFrame({"cdr3_amino_acid": unique_tcrs, sample_id: 1})
#                 # Append unique TCRs with filename to tcr_list 
#                 tcr_list.append(tcr_df)
#         except FileNotFoundError:
#             print("Sample code:", sample_id, "has no file associated")
#             continue  # Move to the next iteration if the file is not found

#         # Print general information
#         elapsed_time = time.time() - start_time
#         progress = idx / len(list_files) * 100
#         print(f"Processing {sample_id} completed in {elapsed_time:.2f} seconds. Progress: {progress:.2f}%")
    
#     #tcr_list = [df_1, df_2, df_3]
#     return tcr_list

# def dataframe_list_update(list_files):

#     tcr_df = pd.DataFrame(columns=["cdr3_amino_acid","sample_file"])

#     # First, fetch the file corresponding to the current sample code
#     for idx, sample_id in enumerate(list_files, start=1):
#         # Start a timer, used later
#         start_time = time.time()

#         # First, fetch the file corresponding to the current sample code
#         try:
#             file_path = f'data/processed_files/{sample_id}_clonotypes.tsv'
#             if not os.path.exists(file_path):
#                 # Handle missing file
#                 continue
            
#             # Read the data in chunks
#             chunks = pd.read_csv(file_path, sep='\t', chunksize=1000000)  # Adjust chunk size as needed
#             # Process each chunk individually
#             for chunk in chunks:
#                 # Get unique TCRs in the chunk
#                 unique_tcrs = chunk["cdr3_amino_acid"].unique()

#                 # Iterate over unique TCRs
#                 for tcr in unique_tcrs:
#                     # Check if TCR already exists in tcr_df
#                     if tcr in tcr_df["cdr3_amino_acid"].values:
#                         # Update sample_file
#                         tcr_df.loc[tcr_df["cdr3_amino_acid"] == tcr, "sample_file"] += f" + {sample_id}"
#                     else:
#                         # Add new TCR to tcr_df
#                         tcr_df = tcr_df.append({"cdr3_amino_acid": tcr, "sample_file": sample_id}, ignore_index=True)

#         except FileNotFoundError:
#             print("Sample code:", sample_id, "has no file associated")
#             continue  # Move to the next iteration if the file is not found

#         # Print general information
#         elapsed_time = time.time() - start_time
#         progress = idx / len(list_files) * 100
#         print(f"Processing {sample_id} completed in {elapsed_time:.2f} seconds. Progress: {progress:.2f}%")
    
#     #tcr_list = [df_1, df_2, df_3]
#     return tcr_df