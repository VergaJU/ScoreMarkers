import functions.functions as fun

define = fun.DefineLabel(adata="../data/nk_cells.h5ad",
                         markers="../data/nk_markers.csv")

df = define.get_label()
print(df)
