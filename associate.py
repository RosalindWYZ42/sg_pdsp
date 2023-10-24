import pandas as pd
import webbrowser
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols

npatlas_ori  = 'NPAtlas_download.tsv'
pdsp_ki_ori = 'KiDatabase.csv'

npatlas = pd.read_csv(npatlas_ori, sep='\t')
pdsp_ki = pd.read_csv(pdsp_ki_ori, sep=',', encoding='ISO 8859-1')

npatlas_smiles = npatlas['compound_smiles']
pdsp_ki_smiles = pdsp_ki['SMILES']

npatlas_c_smiles = []
pdsp_ki_c_smiles = []
#for s in npatlas_smiles:
    #try:
        #cs = Chem.CanonSmiles(s)
        #npatlas_c_smiles.append(cs)
    #except:
        #None
        #Invalid SMILES
#for s in pdsp_ki_smiles:
    #try:
        #cs = Chem.CanonSmiles(s)
        #pdsp_ki_c_smiles.append(cs)
    #except:
        #None
        #Invalid SMILES

pdsp_ki_smiles_map = {}
for _, row in pdsp_ki.iterrows():
    s = row['SMILES']
    if s not in pdsp_ki_smiles_map:
        try:
            cs = Chem.CanonSmiles(s)
            pdsp_ki_c_smiles.append(cs)
            pdsp_ki_smiles_map[s] = [cs]
        except:
            None
    pdsp_ki_smiles_map[s].append(row)

associations = {'SMILES': [], 'npatlas': [], 'pdsp_ki': []}
for _, row in npatlas.iterrows():
    smiles = row['compound_smiles']
    try:
        associations['pdsp_ki'] += [pdsp_ki_smiles_map[smiles]]
        associations['SMILES'] += [smiles]
        associations['npatlas'] += [row]
    except:
        None

associations_df = pd.DataFrame.from_dict(associations)

print(associations_df)
with open('str.html','w') as f:
    associations_df.to_html(f)

filename = 'str.html'
webbrowser.open_new_tab(filename)
