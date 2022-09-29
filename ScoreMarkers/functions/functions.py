# ## This function will contain bla bla bla
import pandas as pd
import scanpy as sc



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
        if not(markers.index.name == "label" and markers.columns[0] == "PoN"):
            print("Wrong csv format")
            return pd.DataFrame()
        else:
            return markers

    def get_markers(self):
        """
        This function split the markers dataframe in 2 different df: Positive markers and negative markers,
        then in drops the columns with index 1 to keep onlt the labels and the markers
        :return: 2 pandas dataframe with positive or negative markers
        """
        markers = self.read_markers()
        if not markers.empty:
            pos_markers = markers.loc[markers["PoN"] == "pos"]
            neg_markers = markers.loc[markers["PoN"] == "neg"]
            pos_markers = pos_markers.drop(columns="PoN")
            neg_markers = neg_markers.drop(columns="PoN")
            return pos_markers, neg_markers
        else:
            print("The dataframe is empty")
            return pd.DataFrame(), pd.DataFrame()
        pass
    def get_adata(self):
        """
        This function read the anndata file
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
        celltypes = markers.index.unique().tolist()
        if len(celltypes) < 1:
            print("No cell labels found!")
            return list()
        return celltypes


    def score_label(self, alpha=0.8, beta=1):
        """
        This function give a score for each label and each cell and return a dataframe with cell barcodes as index
        and the columns with the scores for each label
        :return: pandas dataframe with cells as rows and scores ad columns
        """
        pos_markers, neg_markers = self.get_markers()
        adata = self.get_adata()
        celltypes = self.celltypes()
        cells = pd.DataFrame(adata.obs_names)

        for celltype in celltypes:
            cells["pos"] = 0
            cells["neg"] = 0
            count = 0
            for gene in pos_markers.loc[celltype]:
                if type(gene) == str:
                    try:
                        cells["pos"][((adata[:, gene].X > 0).toarray().flatten())] += 1
                        cells["pos"][((adata[:, gene].X == 0).toarray().flatten())] -= 0
                        count += 1
                    except KeyError:
                        print(f"Positive Marker {gene} for cell type {celltype} not found")
                else:
                    pass
            cells["pos"] = cells["pos"] / count
            count = 0
            for gene in neg_markers.loc[celltype]:
                if type(gene) == str:
                    try:
                        cells["neg"][((adata[:, gene].X > 0).toarray().flatten())] -= 1
                        cells["neg"][((adata[:, gene].X == 0).toarray().flatten())] += 0
                        count += 1
                    except KeyError:
                        print(f"Negative Marker {gene} for cell type {celltype} not found")
                else:
                    pass
            cells["neg"] = cells["neg"] / count
            cells[celltype] = (alpha * cells["pos"]) + (beta * cells["neg"])
        cells = cells.set_index(0)
        return cells
        pass



dc = DefineLabel("../")