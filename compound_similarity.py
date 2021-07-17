#!/usr/bin/env python
# coding: utf-8

# In[14]:


import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import SimilarityMaps
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit import Chem, DataStructs
import pandas as pd
import glob


# In[12]:


def morgan_similarity(df):
    rang = range(len(df))
    for i in rang:
        similarity = []
        cikarilan = df.iloc[i].at['cikarilan_smiles']
        for k in rang:
            optimized = df.iloc[k].at['optimized_smiles']
            mol = Chem.MolFromSmiles(cikarilan)

            if str(optimized) == 'nan':
                similarity.append(0.0)
                continue
        #     print(optimized)
            refmol = Chem.MolFromSmiles(optimized)
            fp1 = AllChem.GetMorganFingerprint(mol,2,useFeatures=True)
            fp2 = AllChem.GetMorganFingerprint(refmol,2,useFeatures=True)
            Tan = metric(fp1,fp2)
            similarity.append(Tan)
        df = pd.concat((df,pd.DataFrame(similarity,columns=[(df.iloc[i].at['cikarilan_smiles'] + '_' + metric.__name__ + '_on_Morgan')])),axis=1)
    df.to_csv(files[j],index=False)


# In[16]:


def topological_similarity(df):
    rang = range(len(df))
    for i in rang:
        similarity = []
        cikarilan = df.iloc[i].at['cikarilan_smiles']
        for k in rang:
            optimized = df.iloc[k].at['optimized_smiles']
            mol = Chem.MolFromSmiles(cikarilan)

            if str(optimized) == 'nan':
                similarity.append(0.0)
                continue
        #     print(optimized)
            refmol = Chem.MolFromSmiles(optimized)
            fp1 = FingerprintMols.FingerprintMol(mol)
            fp2 = FingerprintMols.FingerprintMol(refmol)
            Tan = DataStructs.FingerprintSimilarity(fp1,fp2, metric=metric)
            similarity.append(Tan)
        df = pd.concat((df,pd.DataFrame(similarity,columns=[(df.iloc[i].at['cikarilan_smiles'] + '_' + metric.__name__ + '_on_Topological')])),axis=1)
    df.to_csv(files[j],index=False)


# In[15]:


files = glob.glob("10.07.21/*.csv")


# In[17]:


for j in range(0,len(files)):
    for metr in [DataStructs.TanimotoSimilarity, DataStructs.DiceSimilarity, DataStructs.TanimotoSimilarity, DataStructs.CosineSimilarity, DataStructs.KulczynskiSimilarity,DataStructs.RogotGoldbergSimilarity, DataStructs.BraunBlanquetSimilarity]:
        df = pd.read_csv(files[j])
        metric = metr
        topological_similarity(df)
        
    for metr in [DataStructs.TanimotoSimilarity, DataStructs.DiceSimilarity]:
        df = pd.read_csv(files[j])
        metric = metr
        morgan_similarity(df)

