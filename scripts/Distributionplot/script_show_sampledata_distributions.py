import pandas as pd
import matplotlib.pyplot as plt

dim_list = ['Cancer Study', 'Age', 'Mutation_Count']

data_frame = pd.read_csv("..\data\mutations_merged_filtered_and_processed.csv", sep=';') #usecols=dim_list,

print('\n-----------------------')
print('showing age distribution')

ax = data_frame[['Age']].plot.hist(bins=100,legend=False) #plot(kind='bar', title ="samples per age", legend=False)
ax.set_xlabel("age in years")
ax.set_ylabel("# samples")
plt.show()

print('\n-----------------------')
print('showing mutation count distribution')
#data_frame.sort_values(by='Age', ascending=True)
ax = data_frame.plot('Age','Mutation_Count') #Mutation_Count,DifferentMutatedGenesCount
plt.tight_layout()
plt.show()

print('\n-----------------------')
print('showing samples per study distribution')
ax = data_frame['Cancer Study'].value_counts().plot.bar(legend=True) #Mutation_Count,DifferentMutatedGenesCount
ax.set_ylabel("# samples")
plt.tight_layout()
plt.show()
