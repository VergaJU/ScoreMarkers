#!/usr/local/bin/python
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.preprocessing as sk
import scipy.sparse as sp


class DefineLabel:
    def __init__(self, markers: str, adata: str):
        """
        :param markers: (str) Path to the csv file containing the positive and negative markers for each label
        :param adata: (str) Path to the anndata file with the cells that have to be scored
        """
        self.adata = adata
        self.markers = markers

    def read_markers(self):
        """
        Read the csv path and store the data in a pandas dataframe
        Then CHECK IF THE STRUCTURE IS CORRECT (index name, first columns)
        :return: pandas dataframe with labels and +/- markers
        """
        markers = pd.read_csv(self.markers, index_col=[0])
        if not (markers.index.name == "Label" and markers.columns[0] == "PoN"):
            print("Wrong csv format")
            return pd.DataFrame()
        else:
            return markers

    def get_markers(self):
        """
        This function split the markers dataframe in 2 different df: Positive markers and negative markers,
        then in drops the columns with index 1 to keep onlty the labels and the markers
        :return: 2 pandas dataframe with positive or negative markers
        """
        markers = self.read_markers()
        if not markers.empty:
            pos_markers = markers.loc[markers["PoN"] == "pos"]  # subset positive markers
            neg_markers = markers.loc[markers["PoN"] == "neg"]  # subset negative markers
            pos_markers = pos_markers.drop(columns="PoN")  # drop pos and gen columns
            neg_markers = neg_markers.drop(columns="PoN")
            return pos_markers, neg_markers
        else:
            print("The dataframe is empty")
            return pd.DataFrame(), pd.DataFrame()
        pass

    def get_adata(self):
        """
        This function read the anndata file, normalise it
        :return: adata with sc data (matrix .X and metadata)
        """
        adata = sc.read(self.adata)
        return adata

    def celltypes(self):
        """
        get tuple with cell labels to score
        :return: tuple
        """
        markers = self.read_markers()
        celltypes = markers.index.unique().tolist()  # get celltypes (index column unique)
        if len(celltypes) < 1:
            print("No cell labels found!")
            return list()
        return celltypes

    def score_label(self, alpha=1, beta=1):
        """
        This function give a score for each label and each cell and return a dataframe with cell barcodes as index
        and the columns with the scores for each label
        :return: pandas dataframe with cells as rows and scores ad columns, adata with the sc experiment
        :param alpha: (float) weight positive markers
        :param beta:  (float) weight negative markers
        """
        pos_markers, neg_markers = self.get_markers()
        adata = self.get_adata()
        celltypes = self.celltypes()
        cells = pd.DataFrame(adata.obs_names)  # create dataframe with cells
        adata_tmp = adata.copy()  # create temp file with normalized columns
        sc.pp.normalize_total(adata_tmp)
        #sc.pp.log1p(adata_tmp)
        sc.pp.scale(adata_tmp, zero_center=False)
        alpha = float(alpha)
        beta = float(beta)
        # mean_exp = adata_tmp.X.mean()

        for celltype in celltypes:
            cells["pos"] = 0 # for each cell type create pos and neg columns to score the markers
            cells["neg"] = 0
            count = 0 # count the markers
            for gene in pos_markers.loc[celltype]:
                if type(gene) == str:
                    try:
                        gene_counts = pd.DataFrame(adata_tmp[:, gene].X.toarray().flatten())
                        median_exp = gene_counts[gene_counts[0]>0].median()[0]
                        cells["pos"][((adata_tmp[:, gene].X > 0).toarray().flatten())] += 1 # assign pos score
                        cells["pos"][((adata_tmp[:, gene].X > median_exp).toarray().flatten())] += 1
                        #cells["pos"] += adata_tmp[:, gene].X.toarray().flatten() * 1
                        count += 1
                    except KeyError:
                        print(f"Positive Marker {gene} for cell type {celltype} not found")
                else:
                    pass
            #cells["pos"] = (cells["pos"] - cells["pos"].min()) / (cells["pos"].max() - cells["pos"].min())
            cells["pos"] = cells["pos"] / count # normalise score by the number of markers
            count = 0
            for gene in neg_markers.loc[celltype]:
                if type(gene) == str:
                    try:
                        gene_counts = pd.DataFrame(adata_tmp[:, gene].X.toarray().flatten())
                        median_exp = gene_counts[gene_counts[0]>0].median()[0]
                        cells["neg"][((adata_tmp[:, gene].X > 0).toarray().flatten())] += 1 # assign neg score
                        cells["neg"][((adata_tmp[:, gene].X > median_exp).toarray().flatten())] += 1 # assign neg score
                        #cells["neg"] += adata_tmp[:, gene].X.toarray().flatten() * 1

                        count += 1
                    except KeyError:
                        print(f"Negative Marker {gene} for cell type {celltype} not found")
                else:
                    pass
            #cells["neg"] = (cells["neg"] - cells["neg"].min()) / (cells["neg"].max() - cells["neg"].min())
            cells["neg"] = cells["neg"] / count # normalise score by the number of markers
            cells[celltype] = (alpha * cells["pos"]) - (beta * cells["neg"]) # get final score summing pos and neg by bias
        cells = cells.set_index(0) # cells as index
        del cells["pos"] # drop pos column
        del cells["neg"]
        return cells, adata
        pass

    def plots(self, cells, newfile, abs_values, threshold, thresholdvalue):
        """

        :param cells:
        :param adata:
        :param thresholdlab:
        :return:
        """
        ### prepare dataset and labels:
        cells = cells.iloc[:,:-1]
        labels = cells.columns
        cells = cells.sort_values(by=[labels[1]])
        cells = cells.reset_index()
        abs_values = np.sort(abs_values)[::-1]
        ### Scores distribution with threshold:
        plt.figure(figsize=(9, 4))
        for f in range(len(labels)):
            plt.plot(cells[labels[f]], label=labels[f])
        plt.grid(axis='y',linestyle='--', linewidth=.3)
        plt.title("Scores and threshold distribution")
        plt.xlabel("Cells")
        plt.ylabel("score")
        plt.legend()
        newimage = newfile[:-5] + "_scores_distribution.png"
        plt.savefig(newimage, dpi=300)

        ### absolute distribution scores and percentile
        threshold_label = "Threshold value (" + str(thresholdvalue) + "percentile)"
        plt.figure(figsize=(9,4))
        sns.kdeplot(abs_values, label = "abs(scores)")
        plt.axvline(threshold, color = "r", linestyle="--", label = threshold_label)
        plt.grid(axis='y',linestyle='--', linewidth=.3)
        plt.title("Absolute scores ")
        plt.xlabel("Abs(score)")
        plt.legend()
        newimage = newfile[:-5] + "_absolute_scores.png"
        plt.savefig(newimage, dpi=300)

    def get_label(self, thresholdvalue, newfile="", alpha=1, beta=1, newlabel="new_label", thresholdlab="Other", plot=True):
        """
        This function label each cell with the given label with highest score and save the new anndata file.
        TODO: chech threshold
        TODO: chech threshold
        :param newfile: (str) with the output file name
        :param alpha: (float) weight positive markers
        :param beta:  (float) weight negative markers
        :param newlabel: (str) label of the new variable
        :param thresholdvalue: (int) value of the percentile threshold to label the cells
        :param thresholdlab: (str) name of the thershold label if none of the other labels have an higher score
        :return:
        """
        cells, adata = self.score_label(alpha, beta)
        abs_values = []
        for column in cells.columns:
            abs_values += cells[column].abs().values.tolist()
        abs_values = np.array(abs_values)
        abs_values = abs_values[abs_values != 0]
        thresholdvalue = int(thresholdvalue)
        threshold = np.percentile(abs_values, thresholdvalue)
        cells.insert(loc=0, column=thresholdlab, value=threshold)  # insert threshold column before others
        cells["Label"] = cells.idxmax(axis=1)  # get label by higher score
        # newdf = self.markers  # markers file name
        # newdf = newdf[:-4] + "_results.csv"  # results name (csv)
        # cells.to_csv(newdf)
        cells_lables = cells.iloc[:, -1]  # keep index and labels column
        adata.obs[newlabel] = cells_lables  # add the new labels as metadata to the anndata object
        print(adata.obs[newlabel].value_counts())  # print result
        if newfile == "":  # if newfile is empty take old file name and append "processed.h5ad"
            newfile = self.adata
            newfile = newfile[:-5] + "_processed.h5ad"
        else:
            pass
        adata.write(newfile)  # save new anndata object
        adata.obs[newlabel].to_csv(newfile[:-5] + "_labels.csv")
        newdf = newfile[:-5] + "_results.csv"
        cells.to_csv(newdf)
        if plot:
            self.plots(cells, newfile, abs_values,threshold,thresholdvalue)
        else:
            pass

        return cells, adata



