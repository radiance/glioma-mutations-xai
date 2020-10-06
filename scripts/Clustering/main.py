import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans


# --------------------------------------------------------------------------------
def init_data(n):
    data = pd.read_csv("..\data\mutations_merged_filtered_and_processed.csv", delimiter=';')

    # define dimension list
    dim_list = ['Age', 'TP53', 'IDH1', 'ATRX', 'NF1', 'PIK3R2', 'TERT', 'KMT2A', 'ERBB2', 'ATR', 'BCORL1', 'FLG', 'PIK3CG', 'KDM6A',
                'IL7R', 'RAMP2', 'AXL', 'BRCA1', 'BARD1', 'PBRM1', 'EP300', 'KMT2D', 'RPTOR', 'U2AF1', 'EZH2', 'RELN',
                'PHLPP1', 'SMO', 'CREBBP', 'KIT', 'H3F3B', 'USH2A', 'GNAS', 'SUZ12', 'MAP2K2', 'EPHA3', 'IRF4',
                'SMARCA4', 'PTCH1', 'RYR2', 'EPHA5', 'PDGFRB', 'RAD21', 'CARD11', 'PPM1D', 'IRS2', 'SF3A1', 'SOX1',
                'BLM', 'CDKN2A', 'FLT3', 'TCHH', 'SMC3', 'NTRK1', 'SUSD2', 'FUBP1', 'MED12', 'ASXL1', 'MUC4', 'PRDM1',
                'FGFR3', 'PALB2', 'TEX13D', 'BCL6', 'APC', 'SVIL', 'ASXL2', 'ERBB3', 'MUTYH', 'TNFRSF14', 'OBSCN',
                'RET', 'TNFAIP3', 'HMCN1', 'RICTOR', 'PCLO', 'FOXL2', 'ISM2', 'ABL1', 'CTNNB1', 'FGFR1', 'MPL', 'NF2',
                'TSHR', 'ACVR1', 'MKI67', 'GRIN2A', 'JAK3']

    # create 2D array for clustering (Age and DifferentMutatedGenesCount)
    age = data['Age']
    different_genes_count = data['DifferentMutatedGenesCount']
    size = age.size
    data_2dim = []
    for i in range(0, size):
        data_2dim.append((age.values[i], different_genes_count.values[i]))
    data_2dim = np.asarray(data_2dim)

    # find top genes and sort them in a list
    sorted_genes_list = find_top_n(n, data, dim_list)
    data_top_genes = []
    data_top_genes.append(age)
    for i in range(0, 10):
        data_top_genes.append(data[sorted_genes_list[i][1]])

    # count genes per age and create array
    data_genes_count = np.zeros((100, n+1))
    for row in range(0, 100):
        data_genes_count[row][0] = row
    for col in range(1, n+1):
        for row in range(0, size):
            if data_top_genes[col][row] != '0':
                age = int(data_top_genes[0][row])
                data_genes_count[age][col] += 1

    return data_2dim, sorted_genes_list, data_genes_count


# --------------------------------------------------------------------------------
def find_top_n(n, data, dim_list):
    """
    ----------------------------------------------------------------------------------
    - finds top genes and sort them in a list
    Input:
        n... number of genes, scalar
        data... data array from pandas
        dim_list... a list which contains the names of the dimensions
    Returns:
        sum_list... a list of all genes sorted from the top
    ----------------------------------------------------------------------------------
    """

    list_size = len(dim_list)
    sum_list = []
    for gene in range(1, list_size):
        dim_name = str(dim_list[gene])
        dim_data = data[dim_name]
        sum = 0
        for i in range(0, len(dim_data)):
            if dim_data[i] != '0':
                sum = sum + 1
        entry = [sum, dim_name]
        sum_list.append(entry)

    sum_list.sort(reverse=True)

    for i in range(0, n):
        print("Top-" + str(i+1) + ": " + str(sum_list[i]))

    return sum_list


# --------------------------------------------------------------------------------
def kmeans_library(nr_clusters, max_iter, tol, data_2dim):
    kmeans = KMeans(n_clusters=nr_clusters, max_iter=max_iter, tol=tol, random_state=0).fit(data_2dim)

    # plot
    for k in range(0, nr_clusters):
        #plt.figure(figsize=(12.8, 7.2))
        label = kmeans.labels_[k]
        all_x = data_2dim[:, 0]
        all_y = data_2dim[:, 1]
        plt.plot(kmeans.cluster_centers_[k][0], kmeans.cluster_centers_[k][1], 'r+')
        plt.scatter(all_x, all_y, s=5)

    plt.title('k-means algorithm (library)')
    plt.xlabel('Age')
    plt.ylabel('Different Types of Mutated Genes')
    plt.show()


# --------------------------------------------------------------------------------
def init_k_means(dimension=None, nr_clusters=None, X=None):
    """
    ----------------------------------------------------------------------------------
    - initializes the k_means algorithm
    - chooses suitable initial values
    Input:
        dimension... dimension D of the dataset, scalar
        nr_clusters... scalar
        X... samples, nr_samples x dimension (D)
    Returns:
        initial_centers... initial cluster centers,  D x nr_clusters
    ----------------------------------------------------------------------------------
    """

    # initialize
    initial_centers = [np.ones(dimension)] * nr_clusters

    # calculate initial values
    for k in range(0, nr_clusters):
        initial_centers[k] = X[k]

    initial_centers = np.asarray(initial_centers)

    return initial_centers


# --------------------------------------------------------------------------------
def k_means(X, K, centers_0, max_iter, tol):
    """
    ----------------------------------------------------------------------------------
    - performs the KMeans-algorithm in order to cluster the data into K clusters
    - iteratively updates the cluster centers and classifies all samples after convergence
    Input:
        X... samples, nr_samples x dimension (D)
        K... nr of clusters, scalar
        centers_0... initial cluster centers,  D x nr_clusters
    Returns:
        centers... final centers, D x nr_clusters
        cumulative_distance... cumulative distance over all iterations, nr_iterations x 1
        labels... labels after performing clustering, nr_samples x 1
    ----------------------------------------------------------------------------------
    """

    # compute the dimension
    # D = X.shape[1]
    # assert D == centers_0.shape[0]

    # get N and dimension
    rows, cols = np.shape(X)
    N = rows
    dimension = cols

    # initialize variables
    centers_plus_one = centers_0
    distance_plus_one = 0
    cumulative_distance = []
    classify_distance = np.zeros(K)

    for i in range(0, max_iter):
        centers = centers_plus_one
        distance = distance_plus_one

        # initialize labels
        labels = [[]] * K
        for k in range(0, K):
            labels[k] = []

        # classify data
        for n in range(0, N):
            for k in range(0, K):
                classify_distance[k] = np.linalg.norm((X[n] - centers[k]))
            min_index = np.argmin(classify_distance)
            labels[min_index].append(X[n])
        labels = np.asarray(labels)

        # calculate and update new cluster centers
        for k in range(0, K):
            sum = 0
            for n in range(0, len(labels[k])):
                sum = sum + labels[k][n]
            centers_plus_one[k] = (1 / len(labels[k])) * sum #todo: ZeroDivisionError: division by zero

        # calculate cumulative distance
        distance_plus_one = 0
        for k in range(0, K):
            distance_plus_one = distance_plus_one + np.linalg.norm((labels[k] - centers_plus_one[k]))
        cumulative_distance.append(distance_plus_one)

        if np.linalg.norm(distance - distance_plus_one) < tol:
            # plot final iteration
            plt.figure(figsize=(12.8, 7.2))
            for k in range(0, K):
                label = np.asarray(labels[k])
                all_x = label[:, 0]
                all_y = label[:, 1]
                plt.plot(centers[k][0], centers[k][1], 'r+')
                plt.scatter(all_x, all_y, s=10)

            plt.title('k-means algorithm - iteration: ' + str(i+1) + ' (final)')
            plt.xlabel('Age')
            plt.ylabel('Different Types of Mutated Genes')
            plt.show()

            break

        # plot each iteration
        plt.figure(figsize=(12.8, 7.2))
        for k in range(0, K):
            label = np.asarray(labels[k])
            all_x = label[:, 0]
            all_y = label[:, 1]
            plt.plot(centers[k][0], centers[k][1], 'r+')
            plt.scatter(all_x, all_y, s=10)

        plt.title('k-means algorithm - iteration: ' + str(i+1))
        plt.xlabel('Age')
        plt.ylabel('Different Types of Mutated Genes')
        plt.show()

    return centers, cumulative_distance, labels


# --------------------------------------------------------------------------------
def make_k_means_plots(X, K, centers, cumulative_distance, labels):
    """
    ----------------------------------------------------------------------------------
    - creates the plots for the k-means algorithm
    Input:
        X... samples, nr_samples x dimension (D)
        K... nr of components, scalar
        centers... final centers, D x nr_clusters
        cumulative_distance... cumulative distance over all iterations, nr_iterations x 1
        labels... class labels after performing clustering, nr_samples x 1
    ----------------------------------------------------------------------------------
    """

    # plot cumulative-distance function
    # plt.title("cumulative-distance function")
    # plt.xlabel("iterations")
    # plt.ylabel("cumulative-distance")
    # plt.plot(cumulative_distance)
    # plt.show()

    # plot clustering
    # plt.figure(figsize=(12.8, 7.2))
    for k in range(0, K):
        label = np.asarray(labels[k])
        all_x = label[:, 0]
        all_y = label[:, 1]
        plt.plot(centers[k][0], centers[k][1], 'r+')
        plt.scatter(all_x, all_y, s=5)

    plt.title('k-means algorithm')
    plt.xlabel('Age')
    plt.ylabel('Different Types of Mutated Genes')
    plt.show()


# --------------------------------------------------------------------------------
def gene_distribution_plot(n, data_genes_count, sorted_genes):
    """
    ----------------------------------------------------------------------------------
    - plots the distribution of genes
    Input:
        n... number of genes, scalar
        data_genes_count... data array, 100 x (n+1)
        sorted_genes... a list of all genes sorted from the top, total number of genes x 2
    ----------------------------------------------------------------------------------
    """

    data = data_genes_count.transpose()

    # plot top n genes distribution together
    for col in range(1, n + 1):
        x = data[0]
        y = data[col]
        label = sorted_genes[col - 1][1]
        plt.scatter(x, y, label=label, s=10)
        plt.legend()

    plt.title('Top ' + str(n) + ' Genes Distribution')
    plt.xlabel('Age')
    plt.ylabel('Mutation Count')
    plt.show()

    # plot gene distribution separately
    for col in range(1, n+1):
        x = data[0]
        y = data[col]
        label = sorted_genes[col-1][1]
        plt.scatter(x, y, label=label, s=10)
        plt.legend()

        plt.title(label + ' Gene Distribution')
        plt.xlabel('Age')
        plt.ylabel('Mutation Count')
        plt.show()


# --------------------------------------------------------------------------------
def scatterplot_matrix_for_top_n_genes(n, data_genes_count, sorted_genes):
    """
    ----------------------------------------------------------------------------------
    - creates labeled scatterplot matrix
    Input:
        n... number of genes, scalar
        data_genes_count... data array, 100 x (n+1)
        sorted_genes... a list of all genes sorted from the top, total number of genes x 2
    ----------------------------------------------------------------------------------
    """

    # label data based on previous found age clusters
    labeled_data = np.c_[data_genes_count, np.zeros(100)]
    for i in range(0, 23):
        labeled_data[i][N + 1] = 1
    for i in range(23, 51):
        labeled_data[i][N + 1] = 2
    for i in range(51, len(labeled_data)):
        labeled_data[i][N + 1] = 3

    # convert data to a the type dataframe
    data = pd.DataFrame(labeled_data)

    # create dict and rename columns
    dict = {0: "Age"}
    for i in range(1, n+1):
        dict[i] = sorted_genes[i-1][1]
    dict[n+1] = "labels"
    df = data.rename(columns=dict)

    sns.pairplot(df, diag_kind="kde", hue="labels")  # use "kde" or "hist"
    plt.show()



# --------------------------------------------------------------------------------
# compute best nr of clusters
# --------------------------------------------------------------------------------
def cluster_number_comparison():
    # determining the elbow point in the SSE curve to find best cluster number
    sse = []
    kmeans_kwargs = {
        "init": "random",
        "n_init": 10,
        "max_iter": 300,
        "random_state": 42,
    }
    scaler = StandardScaler()
    scaled_features = scaler.fit_transform(data_genes_count, sorted_genes_list)
    for k in range(2, 20):
        kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
        kmeans.fit(scaled_features)
        sse.append(kmeans.inertia_)

    plt.plot(range(2, 20), sse)
    plt.xticks(range(2, 20))
    plt.xlabel("Number of Clusters")
    plt.ylabel("SSE")
    plt.show()

    # A list holds the silhouette coefficients for each k
    silhouette_coefficients = []

    # Notice you start at 2 clusters for silhouette coefficient
    for k in range(2, 20):
         kmeans = KMeans(n_clusters=k, **kmeans_kwargs)
         kmeans.fit(scaled_features)
         score = silhouette_score(scaled_features, kmeans.labels_)
         silhouette_coefficients.append(score)

    plt.plot(range(2, 20), silhouette_coefficients)
    plt.xticks(range(2, 20))
    plt.xlabel("Number of Clusters")
    plt.ylabel("Silhouette Coefficient")
    plt.show()



# --------------------------------------------------------------------------------
def kmeans_own(nr_clusters, max_iter, tol, dimension, data):
    initial_centers = init_k_means(dimension=dimension, nr_clusters=nr_clusters, X=data)
    centers, cumulative_distance, labels = k_means(data, nr_clusters, initial_centers, max_iter, tol)
    make_k_means_plots(data, nr_clusters, centers, cumulative_distance, labels)


# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# load data
N = 10  # parameter for top N genes
data_2dim, sorted_genes_list, data_genes_count = init_data(N)


# parameters for 2D clustering of Age and DifferentMutatedGenesCount
dimension = 2
nr_clusters = 3
tol = 0.0001  # tolerance
max_iter = 200  # maximum iterations
#cluster_number_comparison()
kmeans_own(nr_clusters, max_iter, tol, dimension, data_2dim) # init clustering
#kmeans_library(nr_clusters, max_iter, tol, data_2dim)


# gene distribution plot
#print(sorted_genes_list)
gene_distribution_plot(N, data_genes_count, sorted_genes_list)

# prepare data for scatterplot matrix by adding label based on previous found age clusters
#scatterplot_matrix_for_top_n_genes(N, data_genes_count, sorted_genes_list)



print("clustering finished.")
