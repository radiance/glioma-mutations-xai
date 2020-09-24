import pandas as pd
import matplotlib.pyplot as plt

data_frame = pd.read_csv("C:\\Dev\\XAI\\HCI-Project-19\\Project\\data\\mutations_stripped.csv", sep=';')

print('\n-----------------------')
print('showing age distribution')

ax = data_frame[['Age']].plot.hist(bins=100,legend=False) #plot(kind='bar', title ="samples per age", legend=False)
ax.set_xlabel("age in years")
ax.set_ylabel("# samples")
plt.show()

print('\n-----------------------')
print('showing mutation count distribution')
ax = data_frame[['DifferentMutatedGenesCount']].plot.hist(bins=100,legend=False) #Mutation_Count,DifferentMutatedGenesCount
ax.set_xlabel("mutationcount per age")
ax.set_ylabel("# samples")
plt.show()
