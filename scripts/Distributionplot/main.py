# -----------------------------------------------------------
# plot distributions about data
#
# created by Jeanquartier 2020
# Released under GNU Public License (GPL)
# -----------------------------------------------------------

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# =============== distribution plot functions =============== #

def print_age_distribution(data_frame):
    print('\n-----------------------')
    print('showing age distribution')
    ax = data_frame[['Age']].plot.hist(bins=100,legend=False) #plot(kind='bar', title ="samples per age", legend=False)
    ax.set_xlabel("age in years")
    ax.set_ylabel("# samples")
    plt.show()

def print_mutation_count(data_frame):
    print('\n-----------------------')
    print('showing mutation count distribution')
    data_frame.sort_values(by='Age', ascending=True)
    ax = data_frame.plot('Age','Mutation_Count') #Mutation_Count,DifferentMutatedGenesCount
    plt.tight_layout()
    plt.show()

def print_samples_per_study(data_frame):
    print('\n-----------------------')
    print('showing samples per study distribution')
    ax = data_frame['Cancer Study'].value_counts().plot.bar(legend=True) #Mutation_Count,DifferentMutatedGenesCount
    ax.set_ylabel("# samples")
    plt.tight_layout()
    plt.show()

def print_top_n_genes(data_frame, topnr):
    print('\n-----------------------')
    print('top ' + str(topnr) + ' genes in filtered data')
    dim_list = data_frame.columns
    toplist = _find_top_n(topnr, data_frame, dim_list)
    for i in range(0, topnr):
        print("Top-" + str(i+1) + ": " + str(toplist[i]))
    #print(toplist[:topnr])
    plt.scatter(*zip(*toplist[:topnr]))
    plt.tight_layout()
    plt.show()

def _find_top_n(n, data, dim_list):
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
    return sum_list


# =============== render plots =============== #

#importantgenes = ["TP53", "IDH1", "ATRX", "H3F3A", "AHNAK2", "SOX1", "SUSD2", "PIK3CA", "TERT", "RYR2", "KMT2A", "KMT2D"] # important genes
data_frame = pd.read_csv("C:\\Dev\\XAI\\xai-glioma-mutations\\data\\mutations_merged_filtered_and_processed-cut.csv", sep=';')
data_frame_only_genes = pd.read_csv("C:\\Dev\\XAI\\xai-glioma-mutations\\data\\mutations_merged_filtered_and_processed-cut.csv", sep=';', usecols=range(5,len(data_frame.columns))) # genes begin in column 6, first 5 columns are sampleids, age, mutation_count etc.

print_age_distribution(data_frame)
#print_mutation_count(data_frame)
#print_samples_per_study(data_frame)
print_top_n_genes(data_frame_only_genes, 87) # min 1, max 87 (only 87 out of 140 selected via pedcbioportal left in csv)
