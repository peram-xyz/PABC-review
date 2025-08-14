#!/usr/bin/env python
# coding: utf-8
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns

import mygene

def aggregate_df(og_df, typ_dict, fxns_dict, oldcolumns):

    header = pd.MultiIndex.from_product([list(typ_dict),
                                     list(fxns_dict)],
                                    names=['Cell Type','Statistic'])

    agg_df = pd.DataFrame(columns=header)

    for newcolumnname in oldcolumns:
        agg_df[newcolumnname] = og_df[oldcolumns[newcolumnname]]

    for cell_line in typ_dict:
        for fxn in fxns_dict:
            agg_df[(cell_line, fxn)] = og_df[typ_dict[cell_line]].apply(fxns_dict[fxn], axis=1)
    
    return agg_df


def plotPCA(df, title):
    expression_data = df.drop(columns=["gene_id", "gene_name", "Geneid"], errors='ignore')
    expression_data = expression_data[~(expression_data == 0).all(axis=1)] # drop rows with all zeros 
    expression_data = expression_data[~(pd.isna(expression_data)).any(axis=1)] # drop rows with any NaN 
    expression_data_T = expression_data.T


    scaler = StandardScaler()
    expression_scaled = scaler.fit_transform(expression_data_T)

    pca = PCA(n_components=2)  
    pca_result = pca.fit_transform(expression_scaled)


    pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
    pca_df["Sample"] = expression_data_T.index
    pca_df["Group"] = pca_df["Sample"].apply(lambda x: "TCGA" if x.split("_")[0] == "tpm" else
                                                        "Placenta" if x.split("_")[0] == "HE8W" or x.split("_")[0] == "HE24W" else
                                                        "MCF7" if x.split("-")[1] == "MCF7" else
                                                        "BEWO" if x.split("-")[1] == "Bewo" else
                                                        "UNDEFINED")


    plt.figure(figsize=(10, 6))
    sns.scatterplot(x='PC1', y='PC2',
                    hue="Group",  
                    data=pca_df,
                    palette="Set2",  
                   )

    plt.title(title)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

