from sklearn.model_selection import KFold
from rdkit import Chem
from representation import Representation
import pandas as pd
import numpy as np


class GetData:
    
    def __init__(self, L, cell_line, descriptors='jtvae', n_fold=5, random_state=0, random_genes=False, csv_file=""):
        """
            Parameters
            -----------
            L: dictionary from L1000CDS_subset.json

            cell_line: cell_id
                params:
                    string: 'VCAP', 'MCF7', 'PC3', etc.

            descriptors: list of descriptors for chemical compounds.
                params:
                    string: 'ecfp', 'ecfp_autoencoder', 'maccs', 'topological', 'shed', 'cats2d', 'jtvae'(default)

            n_fold: number of folds
                params:
                    int: 5(default)

            random_state: random_state for Kfold
                params:
                    int: 0(default)

            random_genes: if it is true, returns random 20 genes from target values
                params:
                    bool: False(default)
                    list of random genes: [118, 919, 274, 866, 354, 253, 207, 667, 773, 563,
                                           553, 918, 934, 81, 56, 232, 892, 485, 30, 53]

            csv_file: if it is not empty, representation data used from this file
                params:
                    string: "<csv_file_path>"
        """
        
        self.L = L
        self.cell_line = cell_line
        self.descriptors = descriptors
        self.n_fold = n_fold
        self.random_state = random_state
        self.random_genes = random_genes
        self.csv_file = csv_file
        self.random_index_list = [118, 919, 274, 866, 354, 253, 207, 667, 773, 563,
                                  553, 918, 934, 81, 56, 232, 892, 485, 30, 53]

        self.X = []
        self.Y = []
        self.perts = []
        self.LmGenes = []
        self.meta_smiles = pd.read_csv('meta_SMILES.csv')

        file_path = 'LmGenes.txt'
        with open(file_path) as fp:
            line = fp.readline()
            while line:
                self.LmGenes.append(line.strip())
                line = fp.readline()

        self.representation_data = None
        if len(csv_file) != 0:
            self.representation_data = pd.read_csv(csv_file)

    def get_regression_data(self):
        unique_smiles = []
        counter = 1
        for pert_id in self.L[self.cell_line]:
            smiles = self.meta_smiles[self.meta_smiles['pert_id'] == pert_id]['SMILES'].values[0]
            mol = Chem.MolFromSmiles(smiles)
            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)

            if canonical_smiles in unique_smiles or len(canonical_smiles) > 120:
                continue
            if self.representation_data is not None:
                feature = self.representation_data[self.representation_data.index == pert_id].values[0]
            else:
                print(str(counter) + ". iteration of " + self.descriptors)
                print('len smiles: ' + str(len(canonical_smiles)))
                counter += 1
                feature = Representation().get_representation(smiles=canonical_smiles, descriptor=self.descriptors)

            unique_smiles.append(canonical_smiles)
            labels = self.L[self.cell_line][pert_id]['chdirLm']
            self.X.append(feature)
            self.Y.append(labels)
            self.perts.append(pert_id)

        x = np.asarray(self.X)
        y = np.asarray(self.Y)

        x_columns = ['SMILES']
        if self.descriptors == 'ecfp':
            for i in range(x.shape[1]-1):
                x_columns.append('ecfp_' + str(i + 1))
        elif self.descriptors == 'ecfp_autoencoder':
            for i in range(x.shape[1]-1):
                x_columns.append('ecfp_autoencoder_' + str(i + 1))
        elif self.descriptors == 'topological':
            for i in range(x.shape[1]-1):
                x_columns.append('topological_' + str(i + 1))
        elif self.descriptors == 'maccs':
            for i in range(x.shape[1]-1):
                x_columns.append('maccs_' + str(i + 1))
        elif self.descriptors == 'jtvae':
            for i in range(x.shape[1]-1):
                x_columns.append('jtvae_' + str(i + 1))
        elif self.descriptors == 'shed':
            for i in range(x.shape[1]-1):
                x_columns.append('shed_' + str(i + 1))
        elif self.descriptors == 'cats2d':
            for i in range(x.shape[1]-1):
                x_columns.append('cats2d_' + str(i + 1))

        x = pd.DataFrame(x, index=self.perts, columns=x_columns)
        y = pd.DataFrame(y, index=self.perts)
        folds = list(KFold(self.n_fold, shuffle=True, random_state=self.random_state).split(x))

        if self.random_genes:
            y_random = []
            for i in self.random_index_list:
                y_random.append(y.iloc[:, i:i + 1])
            df = y_random[0]
            for i in range(len(y_random) - 1):
                df = pd.concat([df, y_random[i + 1]], axis=1)
            y = df

        return x, y, folds

    def get_up_genes(self):
        unique_smiles = []
        counter = 1
        class_dict = {}
        for gene in self.LmGenes:
            class_dict.update({gene: 0})

        for pert_id in self.L[self.cell_line]:
            if 'upGenes' not in self.L[self.cell_line][pert_id]:
                continue

            smiles = self.meta_smiles[self.meta_smiles['pert_id'] == pert_id]['SMILES'].values[0]
            mol = Chem.MolFromSmiles(smiles)
            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)

            if canonical_smiles in unique_smiles or len(canonical_smiles) > 120:
                continue
            if self.representation_data is not None:
                feature = self.representation_data[self.representation_data.index == pert_id].values[0]
            else:
                print(str(counter) + ". iteration of " + self.descriptors)
                print('len smiles: ' + str(len(canonical_smiles)))
                counter += 1
                feature = Representation().get_representation(smiles=canonical_smiles, descriptor=self.descriptors)

            unique_smiles.append(canonical_smiles)
            up_genes = list(set(self.L[self.cell_line][pert_id]['upGenes']))
            class_dict = dict.fromkeys(class_dict, 0)

            for gene in up_genes:
                if gene in class_dict:
                    class_dict.update({gene: 1})

            labels = np.fromiter(class_dict.values(), dtype=int)
            self.X.append(feature)
            self.Y.append(labels)
            self.perts.append(pert_id)

        x = np.asarray(self.X)
        y = np.asarray(self.Y)

        x_columns = ['SMILES']
        if self.descriptors == 'ecfp':
            for i in range(x.shape[1]-1):
                x_columns.append('ecfp_' + str(i + 1))
        elif self.descriptors == 'ecfp_autoencoder':
            for i in range(x.shape[1]-1):
                x_columns.append('ecfp_autoencoder_' + str(i + 1))
        elif self.descriptors == 'topological':
            for i in range(x.shape[1]-1):
                x_columns.append('topological_' + str(i + 1))
        elif self.descriptors == 'maccs':
            for i in range(x.shape[1]-1):
                x_columns.append('maccs_' + str(i + 1))
        elif self.descriptors == 'jtvae':
            for i in range(x.shape[1]-1):
                x_columns.append('jtvae_' + str(i + 1))
        elif self.descriptors == 'shed':
            for i in range(x.shape[1]-1):
                x_columns.append('shed_' + str(i + 1))
        elif self.descriptors == 'cats2d':
            for i in range(x.shape[1]-1):
                x_columns.append('cats2d_' + str(i + 1))

        x = pd.DataFrame(x, index=self.perts, columns=x_columns)
        y = pd.DataFrame(y, index=self.perts)
        folds = list(KFold(self.n_fold, shuffle=True, random_state=self.random_state).split(x))

        if self.random_genes:
            y_random = []
            for i in self.random_index_list:
                y_random.append(y.iloc[:, i:i + 1])
            df = y_random[0]
            for i in range(len(y_random) - 1):
                df = pd.concat([df, y_random[i + 1]], axis=1)
            y = df

        return x, y, folds

    def get_down_genes(self):
        unique_smiles = []
        counter = 1
        class_dict = {}
        for gene in self.LmGenes:
            class_dict.update({gene: 0})

        for pert_id in self.L[self.cell_line]:
            if 'dnGenes' not in self.L[self.cell_line][pert_id]:
                continue

            smiles = self.meta_smiles[self.meta_smiles['pert_id'] == pert_id]['SMILES'].values[0]
            mol = Chem.MolFromSmiles(smiles)
            canonical_smiles = Chem.MolToSmiles(mol, isomericSmiles=False)

            if canonical_smiles in unique_smiles or len(canonical_smiles) > 120:
                continue
            if self.representation_data is not None:
                feature = self.representation_data[self.representation_data.index == pert_id].values[0]
            else:
                print(str(counter) + ". iteration of " + self.descriptors)
                print('len smiles: ' + str(len(canonical_smiles)))
                counter += 1
                feature = Representation().get_representation(smiles=canonical_smiles, descriptor=self.descriptors)

            unique_smiles.append(canonical_smiles)
            dn_genes = list(set(self.L[self.cell_line][pert_id]['dnGenes']))
            class_dict = dict.fromkeys(class_dict, 0)
            for gene in dn_genes:
                if gene in class_dict:
                    class_dict.update({gene: 1})

            labels = np.fromiter(class_dict.values(), dtype=int)
            self.X.append(feature)
            self.Y.append(labels)
            self.perts.append(pert_id)

        x = np.asarray(self.X)
        y = np.asarray(self.Y)

        x_columns = ['SMILES']
        if self.descriptors == 'ecfp':
            for i in range(x.shape[1]-1):
                x_columns.append('ecfp_' + str(i + 1))
        elif self.descriptors == 'ecfp_autoencoder':
            for i in range(x.shape[1]-1):
                x_columns.append('ecfp_autoencoder_' + str(i + 1))
        elif self.descriptors == 'topological':
            for i in range(x.shape[1]-1):
                x_columns.append('topological_' + str(i + 1))
        elif self.descriptors == 'maccs':
            for i in range(x.shape[1]-1):
                x_columns.append('maccs_' + str(i + 1))
        elif self.descriptors == 'jtvae':
            for i in range(x.shape[1]-1):
                x_columns.append('jtvae_' + str(i + 1))
        elif self.descriptors == 'shed':
            for i in range(x.shape[1]-1):
                x_columns.append('shed_' + str(i + 1))
        elif self.descriptors == 'cats2d':
            for i in range(x.shape[1]-1):
                x_columns.append('cats2d_' + str(i + 1))

        x = pd.DataFrame(x, index=self.perts, columns=x_columns)
        y = pd.DataFrame(y, index=self.perts)
        folds = list(KFold(self.n_fold, shuffle=True, random_state=self.random_state).split(x))

        if self.random_genes:
            y_random = []
            for i in self.random_index_list:
                y_random.append(y.iloc[:, i:i + 1])
            df = y_random[0]
            for i in range(len(y_random) - 1):
                df = pd.concat([df, y_random[i + 1]], axis=1)
            y = df

        return x, y, folds
