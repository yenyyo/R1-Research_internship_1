import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from collections import defaultdict
import numpy as np
from scipy.stats import fisher_exact
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture


def graph_clusters_subplt(data, num_clusters=5, plots_per_row=5, y_axis_as_percentage=False):
    """
    Generates subplots for given clusters with scatter plots of file sizes vs. counts or percentages of TCRs.

    Parameters:
    - data: A tuple containing:
        - unique_clusters: List of unique cluster IDs
        - patient_tcrs: DataFrame with columns ['combined', 'file_id']
        - patient_cluster: DataFrame with columns ['cluster_id', 'combined']
        - cluster_counts: DataFrame with columns ['cluster_id', 'total_tcrs']
        - file_sizes: DataFrame with columns ['Unnamed: 0', 'unique_tcrs']
    - num_clusters: Number of clusters to plot (default is 5)
    - plots_per_row: Number of plots per row in the figure (default is 5)
    - y_axis_as_percentage: Boolean flag to determine if y-axis should be percentage (default is False, meaning count)

    Returns:
    - cluster_data: Dictionary with cluster ID as key and a tuple of lists (sizes, y_values) as values
    """

    unique_clusters, patient_tcrs, patient_cluster, cluster_counts, file_sizes = data

    num_rows = num_clusters // plots_per_row + (1 if num_clusters % plots_per_row != 0 else 0)

    fig, axs = plt.subplots(num_rows, plots_per_row, figsize=(plots_per_row * 5, num_rows * 5))
    axs = axs.flatten()

    cluster_data = {}

    for idx, cluster in enumerate(unique_clusters[:num_clusters]):
        file_id_counts = defaultdict(int)
        file_id_tcr_map = defaultdict(set)
        filtered_patients = patient_cluster.loc[patient_cluster["cluster_id"] == cluster]
        
        for tcr in filtered_patients["combined"]:
            files = patient_tcrs.loc[patient_tcrs["combined"] == tcr]["file_id"]
            for file_id_list in files:
                for file_id in file_id_list.split(','):
                    file_id_counts[file_id] += 1
                    file_id_tcr_map[file_id].add(tcr)

        sizes = []
        y_values = []
        colors = []
        file_ids = []

        unique_tcrs = list({tcr for tcrs in file_id_tcr_map.values() for tcr in tcrs})
        tcr_color_map = {tcr: idx for idx, tcr in enumerate(unique_tcrs)}
        tcrs_in_cluster = cluster_counts.loc[cluster_counts["cluster_id"] == cluster]["total_tcrs"].iloc[0]

        for file_id, count in file_id_counts.items():
            size = file_sizes.loc[file_sizes["Unnamed: 0"] == file_id]["unique_tcrs"].values
            if size.size > 0:
                sizes.append(size[0])
                if y_axis_as_percentage:
                    y_values.append((count / tcrs_in_cluster) * 100)
                else:
                    y_values.append(count)
                tcrs_for_file = list(file_id_tcr_map[file_id])
                color_idx = tcr_color_map[tcrs_for_file[0]]
                colors.append(color_idx)
                
                # Append the current file_id to the file_ids list
                file_ids.append(file_id)

        cluster_data[cluster] = (sizes, y_values, file_ids)

        ax = axs[idx]
        scatter = ax.scatter(sizes, y_values, c=colors, cmap='viridis', marker='o', alpha=0.7, edgecolors='black')
        ylabel = 'Percentage of TCRs in Cluster (%)' if y_axis_as_percentage else 'Count of File ID Occurrences'
        ax.set_title(f'Cluster {cluster}\nSize: {tcrs_in_cluster} TCRs in {len(file_id_counts)} files', fontsize=10)
        ax.set_xlabel('File Size (unique_tcrs)', fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.grid(True, linestyle='--', linewidth=0.5, color='gray')

        ax.tick_params(bottom=True, top=True, left=True, right=True, labelsize=8)
        ax.set_xlim(min(sizes) - 0.1 * (max(sizes) - min(sizes)), max(sizes) + 0.1 * (max(sizes) - min(sizes)))
        ax.set_ylim(0, max(y_values) + 0.1 * max(y_values))

    plt.tight_layout()
    plt.colorbar(scatter, ax=axs, orientation='horizontal', fraction=0.02, pad=0.04)
    plt.show()
    
    return cluster_data


def plt_single_graph(data, number_to_plot, rainbow):
    """
    Function to print individual plots (not in subplot).

    Args:
        data: A tuple containing data for analysis.
        number_to_plot: Number of clusters to plot.
        colors: If 1, uses colors associated with IDs, else uses blue.
    """

    unique_clusters, patient_tcrs, patient_cluster, cluster_counts, file_sizes = data

    for cluster in unique_clusters[:number_to_plot]:
        file_id_counts = defaultdict(int)

        filtered_patients = patient_cluster.loc[patient_cluster["cluster_id"] == cluster]

        for tcr in filtered_patients["combined"]:
            files = patient_tcrs.loc[patient_tcrs["combined"] == tcr]["file_id"]
            for file_id_list in files:
                for file_id in file_id_list.split(','):
                    file_id_counts[file_id] += 1

        sizes, counts, file_ids = [], [], []
        for file_id, count in file_id_counts.items():
            size = file_sizes.loc[file_sizes["Unnamed: 0"] == file_id]["unique_tcrs"].values
            if size.size > 0:
                sizes.append(size[0])
                counts.append(count)
                file_ids.append(file_id)

        tcrs_in_cluster = cluster_counts.loc[cluster_counts["cluster_id"] == cluster]["total_tcrs"].iloc[0]

        plt.figure(figsize=(10, 5))

        if rainbow == 1:
            cmap = plt.get_cmap("viridis")
            colors = cmap(np.linspace(0, 1, len(file_ids)))
            for i, file_id in enumerate(file_ids):
                plt.scatter(sizes[i], counts[i], marker='o', color=colors[i], alpha=0.7, edgecolors='black', label=file_id)
        else:
            plt.scatter(sizes, counts, marker='o', c='blue', alpha=0.7, edgecolors='black')

        plt.title(f'Cluster {cluster} - File Size vs Count\nSize: {tcrs_in_cluster} distinct TCRs in {len(file_id_counts)} different files', fontsize=14)
        plt.xlabel('File Size (unique_tcrs)', fontsize=12)
        plt.ylabel('Count of File ID Occurrences', fontsize=12)
        plt.grid(True, linestyle='--', linewidth=0.5, color='gray')
        plt.tick_params(bottom=True, top=True, left=True, right=True, labelsize=10)
        plt.xlim(min(sizes) - 0.1 * (max(sizes) - min(sizes)), max(sizes) + 0.1 * (max(sizes) - min(sizes)))
        plt.ylim(0, max(counts) + 0.1 * max(counts))
        plt.tight_layout()
        plt.show()

def apply_kmeans_clustering(cluster_data, n_clusters=2, plots_per_row=5):
    """
    Applies K-means clustering to each cluster's data and plots the results.

    Parameters:
    - cluster_data: Dictionary with cluster ID as key and a tuple of lists (sizes, percentages) as values
    - n_clusters: Number of clusters for K-means (default is 2)
    - plots_per_row: Number of plots per row in the figure (default is 5)
    """
    num_clusters = len(cluster_data)
    num_rows = num_clusters // plots_per_row + (1 if num_clusters % plots_per_row != 0 else 0)

    fig, axs = plt.subplots(num_rows, plots_per_row, figsize=(plots_per_row * 5, num_rows * 5))
    axs = axs.flatten()

    for idx, (cluster, (sizes, percentages)) in enumerate(cluster_data.items()):
        data = np.column_stack((sizes, percentages))
        
        # Apply K-means clustering with specified number of clusters
        kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(data)
        labels = kmeans.labels_
        centers = kmeans.cluster_centers_

        # Plot the data points and color them by cluster
        ax = axs[idx]
        scatter = ax.scatter(sizes, percentages, c=labels, cmap='coolwarm', marker='o', edgecolors='black')
        ax.scatter(centers[:, 0], centers[:, 1], c='yellow', s=200, alpha=0.75, marker='X')

        # Optionally, draw a line between the centroids (for visual aid)
        for i in range(n_clusters - 1):
            for j in range(i + 1, n_clusters):
                ax.plot([centers[i, 0], centers[j, 0]], [centers[i, 1], centers[j, 1]], '--', color='red')

        ax.set_title(f'Cluster {cluster} - K-means Clustering with k={n_clusters}', fontsize=10)
        ax.set_xlabel('File Size (unique_tcrs)', fontsize=8)
        ax.set_ylabel('Percentage of TCRs in Cluster (%)', fontsize=8)
        ax.grid(True, linestyle='--', linewidth=0.5, color='gray')

        ax.tick_params(bottom=True, top=True, left=True, right=True, labelsize=8)

    # Remove empty subplots
    for j in range(idx + 1, len(axs)):
        fig.delaxes(axs[j])

    plt.tight_layout()
    plt.show()

def apply_dbscan_clustering(cluster_data, eps=0.05, min_samples=5, num_clusters=2):
    """
    function recieves a cluster_data parameter:
    This contains information of sizes and percentages for each of the clusters
    Its the return value of graph_clusters_subplt

    It applies dbscan (density based clustering) clustering to each of the plots
    """
    # Limit the number of clusters to the specified value
    selected_clusters = list(cluster_data.items())[:num_clusters]

    fig, axs = plt.subplots(1, num_clusters, figsize=(num_clusters * 5, 5))
    if num_clusters == 1:
        axs = [axs]  # Ensure axs is a list if there is only one subplot

    for idx, (cluster, (sizes, percentages)) in enumerate(selected_clusters):
        data = np.column_stack((sizes, percentages))
        
        # Normalize data
        data = (data - data.mean(axis=0)) / data.std(axis=0)
        
        # Apply DBSCAN clustering
        dbscan = DBSCAN(eps=eps, min_samples=min_samples).fit(data)
        labels = dbscan.labels_

        ax = axs[idx]
        unique_labels = set(labels)
        colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]

        for k, col in zip(unique_labels, colors):
            class_member_mask = (labels == k)
            xy = data[class_member_mask]
            ax.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=tuple(col), markeredgecolor='k', markersize=6)

        ax.set_title(f'Cluster {cluster} - DBSCAN Clustering')
        ax.set_xlabel('Normalized File Size (unique_tcrs)')
        ax.set_ylabel('Normalized Percentage of TCRs in Cluster (%)')
        ax.grid(True, linestyle='--', linewidth=0.5, color='gray')

    plt.tight_layout()
    plt.show()


def find_interesting_clusters(cluster_data, threshold = 1.5, min_mean_diff=0.5, min_var_ratio=2.0, n_components=2):
    """
    Identifies interesting clusters based on slope differences using GMM.

    Args:
        cluster_data (dict): Dictionary containing data for each cluster.
            - 'counts' (np.array): Array of x-axis values (counts).
            - 'sizes' (np.array): Array of y-axis values (sizes).
        
        min_mean_diff (float, optional): Minimum difference in means for interesting clusters. Defaults to 0.5.
        min_var_ratio (float, optional): Minimum ratio of variances for interesting clusters. Defaults to 2.0.
        n_components (int, optional): Number of components for GMM. Defaults to 2.

    Returns:
        tuple: Tuple containing a list of interesting cluster IDs and a dictionary of GMM labels and means.
    """
    interesting_clusters = []
    cluster_predictions = []

    # Loop through each cluster
    for cluster_id, data in cluster_data.items():
        counts = data[0]
        sizes = data[1]
        file_id = data[2]

        # Calculate slopes
        slopes = np.array(sizes) / np.array(counts)
        slopes = slopes * 1000 # for some reason, as slopes are too small it dosent compute them well

        # Fit GMM to slopes
        gmm = GaussianMixture(n_components=n_components, random_state=42)
        predictions = gmm.fit_predict(slopes.reshape(-1, 1))

        # Extract means and variances
        means = gmm.means_.squeeze()
        variances = gmm.covariances_.squeeze()

        # Check for significant difference
        #rank them, and get top-k. determine threshold based = higher slope = lower slope x 1.5
        if (means[1] >= threshold * means[0]): #and abs(variances[0] / variances[1]) >= min_var_ratio :
            interesting_clusters.append(cluster_id)
            cluster_prediction = pd.DataFrame({'slopes':slopes, 'predictions' : predictions, 'percentages': sizes, 'counts': counts, 'file_id': file_id})
            cluster_predictions.append(cluster_prediction)


        # print(means)
    return interesting_clusters, cluster_predictions


def plot_gaussian_clusters(cluster_data, gmm_results, plots_per_row=5):
    """
    Plots interesting clusters and colors points based on GMM component labels.

    Args:
        cluster_data (dict): Dictionary containing data for each cluster.
        gmm_results (dict): Dictionary containing GMM results for interesting clusters.
        plots_per_row (int, optional): Number of plots per row. Defaults to 5.
    """
    interesting_clusters = list(gmm_results.keys())
    num_clusters = len(interesting_clusters)
    num_rows = num_clusters // plots_per_row + (1 if num_clusters % plots_per_row != 0 else 0)

    fig, axs = plt.subplots(num_rows, plots_per_row, figsize=(plots_per_row * 5, num_rows * 5))
    
    # Ensure axs is always a flattened array for consistent handling
    if num_clusters == 1:
        axs = np.array([axs])
    else:
        axs = axs.flatten()

    for idx, cluster_id in enumerate(interesting_clusters):
        counts = cluster_data[cluster_id][0]
        sizes = cluster_data[cluster_id][1]
        slopes = gmm_results[cluster_id]['slopes']
        labels = gmm_results[cluster_id]['labels']
        means = gmm_results[cluster_id]['means']

        ax = axs[idx]
        scatter = ax.scatter(np.arange(len(slopes)), slopes, c=labels, cmap='coolwarm', marker='o', edgecolors='black')
        
        # Plot GMM means as horizontal lines
        for mean in means:
            ax.axhline(y=mean, color='red', linestyle='--')

        ax.set_title(f'Cluster {cluster_id} - GMM Components', fontsize=10)
        ax.set_xlabel('Slope Index', fontsize=8)
        ax.set_ylabel('Slope', fontsize=8)
        ax.grid(True, linestyle='--', linewidth=0.5, color='gray')
        ax.tick_params(bottom=True, top=True, left=True, right=True, labelsize=8)

    # Remove empty subplots
    for j in range(idx + 1, len(axs)):
        fig.delaxes(axs[j])

    plt.tight_layout()
    plt.show()

    
def plot_clusters(cluster_data, interesting_clusters, plots_per_row=5, y_axis_label='Count of File ID Occurrences'):
    """
    Plots only the specified interesting clusters with sizes and counts/percentages.

    Parameters:
    - cluster_data: Dictionary with cluster ID as key and a tuple of lists (sizes, y_values) as values
    - interesting_clusters: List of cluster IDs to plot
    - plots_per_row: Number of plots per row in the figure (default is 5)
    - y_axis_label: Label for the y-axis (default is 'Count of File ID Occurrences')
    """
    num_clusters = len(interesting_clusters)
    num_rows = num_clusters // plots_per_row + (1 if num_clusters % plots_per_row != 0 else 0)

    fig, axs = plt.subplots(num_rows, plots_per_row, figsize=(plots_per_row * 5, num_rows * 5))
    axs = axs.flatten()

    for idx, cluster in enumerate(interesting_clusters):
        if cluster not in cluster_data:
            print(f"Cluster {cluster} not found in cluster_data.")
            continue
        
        sizes, y_values = cluster_data[cluster]
        ax = axs[idx]
        scatter = ax.scatter(sizes, y_values, c=y_values, cmap='viridis', marker='o', alpha=0.7, edgecolors='black')
        ax.set_title(f'Cluster {cluster}', fontsize=10)
        ax.set_xlabel('File Size (unique_tcrs)', fontsize=8)
        ax.set_ylabel(y_axis_label, fontsize=8)
        ax.grid(True, linestyle='--', linewidth=0.5, color='gray')

        ax.tick_params(bottom=True, top=True, left=True, right=True, labelsize=8)
        ax.set_xlim(min(sizes) - 0.1 * (max(sizes) - min(sizes)), max(sizes) + 0.1 * (max(sizes) - min(sizes)))
        ax.set_ylim(0, max(y_values) + 0.1 * max(y_values))

    # Hide any unused subplots
    for j in range(idx + 1, len(axs)):
        fig.delaxes(axs[j])

    plt.tight_layout()
    plt.colorbar(scatter, ax=axs, orientation='horizontal', fraction=0.02, pad=0.04)
    plt.show()

def plot_predictions(gmm_results, n_columns=5):
    
    num_plots = len(gmm_results)  # Number of plots
    cols = n_columns  # Number of columns in the subplot grid
    rows = (num_plots + 1) // cols  # Calculate the number of rows needed

    # Create a figure and a grid of subplots
    fig, axes = plt.subplots(rows, cols, figsize=(25, rows * 5))  # Adjust figsize as needed


    # Flatten axes array for easy indexing
    axes = axes.flatten()

    # Loop through the gmm_results and create scatter plots
    for idx, data in enumerate(gmm_results):
        scatter = sns.scatterplot(data=data, x='counts', y='percentages', hue='predictions', 
                                marker='o', alpha=0.7, edgecolor='black', ax=axes[idx])
        axes[idx].set_title(f'Plot {idx + 1}')
        
    # Remove any unused subplots if gmm_results is not a perfect multiple of cols
    for ax in axes[num_plots:]:
        ax.remove()

    plt.tight_layout()  # Adjust layout to prevent overlap
    plt.show()