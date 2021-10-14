#!/usr/bin/env python3

import numpy as np; np.random.seed(0)
import seaborn as sns; sns.set_theme()
import pandas as pd
import glob
from matplotlib import pyplot as plt


files = glob.glob("10.07.21/*.csv")
for file in files:
    data = pd.read_csv(file,index_col="optimized_smiles")
    data.drop(["baslangic_smiles", "distance", "cikarilan_smiles"], axis=1,inplace=True)
    ax = sns.heatmap(data)
    ax.get_figure().savefig(files.split('.csv')[0] + ".png")
