import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

# load data
t15df = pd.read_csv(os.path.join(os.path.dirname(__file__), "..\..\data\mutations_merged_grouped_ranked_filtered_top15_transposed.csv"), delimiter=';', header=0,)
t15df.columns = t15df.columns.str.strip()
df = t15df.T

agegroups = t15df.iloc[0:0,1:10].columns.tolist()
print("agegroups:")
print(agegroups)

# important genes in paper: 12 selected
importantgeneslist = ["TP53", "TERT", "IDH1", "AHNAK2", "ATRX", "KMT2D", "RYR2", "SOX1", "SUSD2", "H3F3A", "PIK3CA", "KMT2C"]

# ---------------------------
# # get top selected genes per age group
#print(importantgeneslist) #print(genelist)
#print(t15df)
#
#each row is a gene with listing values for age group
#removing first item in list because label
AHNAK2 = t15df.iloc[[1]].values.tolist()[0]
AHNAK2.pop(0)
TP53 = t15df.iloc[[2]].values.tolist()[0]
TP53.pop(0)
H3F3A = t15df.iloc[[4]].values.tolist()[0]
H3F3A.pop(0)
KMT2D = t15df.iloc[[10]].values.tolist()[0]
KMT2D.pop(0)
SUSD2 = t15df.iloc[[11]].values.tolist()[0]
SUSD2.pop(0)
SOX1 = t15df.iloc[[12]].values.tolist()[0]
SOX1.pop(0)
ATRX = t15df.iloc[[14]].values.tolist()[0]
ATRX.pop(0)
IDH1 = t15df.iloc[[22]].values.tolist()[0]
IDH1.pop(0)
PIK3CA = t15df.iloc[[26]].values.tolist()[0]
PIK3CA.pop(0)
TERT = t15df.iloc[[28]].values.tolist()[0]
TERT.pop(0)
RYR2 = t15df.iloc[[35]].values.tolist()[0]
RYR2.pop(0)
KMT2C = t15df.iloc[[15]].values.tolist()[0]
KMT2C.pop(0)
print(KMT2D)
print(KMT2C)

rows = 3
cols = 4
i = 0
xlabels = agegroups
xlabels.insert(0, '')
fig, ax = plt.subplots(rows, cols, sharex='col', sharey='row')
for row in range(rows):
    for col in range(cols):

        ax[row, col].invert_yaxis()
        ax[row, col].grid(True)
        ax[row, col].set_xticklabels(xlabels)
        ax[row, col].locator_params(axis = 'x', nbins = 9)
        ax[row, col].locator_params(axis = 'y', nbins = 5)
        ax[row, col].tick_params(axis='x', labelsize='5')
        ax[row, col].set_xlabel('age groups', fontsize='8')
        ax[row, col].set_ylabel('rankings', fontsize='8')

        if i < 12:
            ax[row,col].set_title(importantgeneslist[i])
            ax[row,col].set_ylim(bottom=20, top=0)
            i = i+1

#add the linecharts to the grid plot
ax[0, 0].plot(TP53, marker='o', linestyle='solid', linewidth=1, markersize=3, color='green')
ax[0, 0].fill_between(range(0,9), TP53, 20, alpha=0.2)

ax[0, 1].plot(TERT, marker='o', linestyle='solid', linewidth=1, markersize=3, color='green')
ax[0, 1].fill_between(range(0,9), TERT, 20, alpha=0.2)

ax[0, 2].plot(IDH1, marker='o', linestyle='solid', linewidth=1, markersize=3, color='blue')
ax[0, 2].fill_between(range(0,9), IDH1, 20, alpha=0.2)

ax[0, 3].plot(AHNAK2, marker='o', linestyle='solid', linewidth=1, markersize=3, color='orange')
ax[0, 3].fill_between(range(0,9), AHNAK2, 20, alpha=0.2)

ax[1, 0].plot(ATRX, marker='o', linestyle='solid', linewidth=1, markersize=3, color='brown')
ax[1, 0].fill_between(range(0,9), ATRX, 20, alpha=0.2)

ax[1, 1].plot(KMT2D, marker='o', linestyle='solid', linewidth=1, markersize=3, color='red')
ax[1, 1].fill_between(range(0,9), KMT2D, 20, alpha=0.2)

ax[1, 2].plot(RYR2, marker='o', linestyle='solid', linewidth=1, markersize=3, color='darkorange')
ax[1, 2].fill_between(range(0,9), RYR2, 20, alpha=0.2)

ax[1, 3].plot(SOX1, marker='o', linestyle='solid', linewidth=1, markersize=3, color='darkgreen')
ax[1, 3].fill_between(range(0,9), SOX1, 20, alpha=0.2)

ax[2, 0].plot(SUSD2, marker='o', linestyle='solid', linewidth=1, markersize=3, color='darkred')
ax[2, 0].fill_between(range(0,9), SUSD2, 20, alpha=0.2)

ax[2, 1].plot(H3F3A, marker='o', linestyle='solid', linewidth=1, markersize=3, color='darkblue')
ax[2, 1].fill_between(range(0,9), H3F3A, 20, alpha=0.2)

ax[2, 2].plot(PIK3CA, marker='o', linestyle='solid', linewidth=1, markersize=3, color='turquoise')
ax[2, 2].fill_between(range(0,9), PIK3CA, 20, alpha=0.2)

ax[2, 3].plot(KMT2C, marker='o', linestyle='solid', linewidth=1, markersize=3, color='violet')
ax[2, 3].fill_between(range(0,9), KMT2C, 20, alpha=0.2)

plt.tight_layout()
plt.grid(True)
plt.ylim(bottom=20, top=0)
plt.xlabel("age groups")
plt.ylabel("rankings")
plt.show()
