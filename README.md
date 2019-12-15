# ChemRep
Representations of chemical compounds

## Example Usage
GetData class has the following methods.

- get_regression_data
- get_up_genes
- get_down_genes

```
from get_data import GetData
import json

with open('L1000CDS_subset.json', 'r') as f:
    L = json.load(f)

obj = GetData(L=L, cell_line='VCAP', descriptors='jtvae', n_fold=5, random_state=42,
              random_genes=False, csv_file='VCAP_jtvae.csv')

x, y, folds = obj.get_down_genes()
x.drop(['SMILES'], axis=1, inplace=True)

for i, (trn, val) in enumerate(folds) :
    print(i+1, "fold.")
    
    trn_x = x.iloc[trn, :]
    trn_y = y.iloc[trn, :]
    val_x = x.iloc[val, :]
    val_y = y.iloc[val, :]
    
    print('Train shapes:', trn_x.shape, trn_y.shape)
    print('Test shapes:', val_x.shape, val_y.shape)

    model.fit(trn_x, trn_y)
    ...


```

Autoencoderdan tek tek representation almak uzun sürdüğü için önceden csv dosyasına kaydedilmiş
representationları kullanabilmek için GetData sınıfına `csv_file` parametresi eklendi.
Buraya `pert_id`,`SMILES` ve representation kolonları bulunan bir csv dosyası
verildiğinde bu dosyadan yararlanarak hızlı bir şekilde feature, target ve fold değerlerini GetData
sınıfındaki metodları kullanarak alabiliyoruz. Eğer bu parametre kullanılmazsa istenilen representation
tek tek çıkartılır.