from typing import List
import pandas as pd
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from socialgene.neo4j.neo4j import GraphDriver
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rich.progress import Progress
from typing import List
from rdkit import DataStructs
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
import os

#Get a 1024-dimensional vector that represents the Morgan fingerprint of the molecule.
def morgan_fingerprint(rdkitmol, radius=2, nBits=1024):
    return AllChem.GetMorganFingerprintAsBitVect(
        rdkitmol, useChirality=True, radius=radius, nBits=nBits, bitInfo={}
    )

#Construct a molecule from a SMILES string
def read_mol(input):
    m = Chem.MolFromSmiles(input, sanitize=True)
    m.UpdatePropertyCache()
    return m



def rehydrate(x):
    bv = ExplicitBitVect(1024)
    bv.FromBase64(x)
    return bv


def compare_one_to_many(
    query_bitvector: ExplicitBitVect,
    target_bitvector_list: List[ExplicitBitVect],
    target_uid_list: List[int],
    cutoff: float = 0.7,
) -> List[int]:
    """
    Compare a query bitvector to a list of target bitvectors using the BulkTanimotoSimilarity method.
    Return a list of target_uids for targets whose similarity score is greater than or equal to the cutoff.

    Args:
        query_bitvector: The query bitvector to compare against.
        target_bitvector_list: A list of target bitvectors to compare against (element order should reflect the order of target_uid_list).
        target_uid_list: A list of target uids corresponding to the target bitvectors .
        cutoff: The similarity score cutoff. Defaults to 0.7.

    Returns:
        A list of target uids for targets whose similarity score is greater than or equal to the cutoff.
    """
    return [
        target_uid_list[i]
        for i, x in enumerate(
            DataStructs.BulkTanimotoSimilarity(query_bitvector, target_bitvector_list)
        )
        if x >= cutoff
    ]




def fingerprint_from_smiles(x):
    try:
        return morgan_fingerprint(Chem.MolFromSmiles(x))
    except:
        return None

def tidy(df, SMILES='SMILES'):
    df = df[df[SMILES].notna()]
    df = df.drop_duplicates(subset=SMILES)
    df["morgan"] = df[SMILES].apply(fingerprint_from_smiles)
    df = df[~df["morgan"].isnull()]
    return df


#################################################
# Calculate Morgan Fingerprints for npatlas
#################################################
# with GraphDriver() as db:
#     res1 = db.run(
#         """
#             MATCH (n:npatlas)
#             RETURN n.uid as uid, n.smiles as smiles
#             """,
#     ).values()

# # store morgan fingerprints as base64
# with GraphDriver() as db:
#     for k, v in res1:
#         _ = db.run(
#             """
#                 WITH $input as inputs
#                 UNWIND inputs as input
#                 MATCH (n:npatlas {uid: input.uid})
#                 SET n.morgan_base64 = input.morgan_base64
#                 """,
#             input={
#                 "uid": k,
#                 "morgan_base64": morgan_fingerprint(read_mol(v)).ToBase64(),
#             },
#         )

# # store morgan fingerprints as vector
# with GraphDriver() as db:
#     for k, v in res1:
#         _ = db.run(
#             """
#                 WITH $input as inputs
#                 UNWIND inputs as input
#                 MATCH (n:npatlas {uid: input.uid})
#                 SET n.morgan_base64 = input.morgan_base64
#                 """,
#             input={"uid": k, "morgan_base64": list(morgan_fingerprint(read_mol(v)))},
#         )
#################################################
#################################################

# retrieve all npatalas compounds' morgan fingerprints from the database
with GraphDriver() as db:
    res = db.run(
        """
            MATCH (n:npatlas) where n.morgan_base64 is not null 
            RETURN n.uid as uids, n.morgan_base64 as morgan_base64
            """,
    ).to_df()
    
    
# START PLAYING WITH CSVS FROM HERE
# res ABOVE WOULD BE A PANDAS DATAFRAME WITH COLUMNS ["UIDS","MORGAN_BASE64"]

    

# convert base64 to bitvector
res["bitvec"] = res["morgan_base64"].apply(rehydrate)
del res["morgan_base64"]

df = pd.read_csv(
    'KiDatabase.csv', sep=',', encoding='ISO 8859-1')

df = os.path.join(
    os.path.expanduser('~'), 'socialgene', 'addons', 'KiDatabase.scv')

# remove rows without SMILES
df = df[df["SMILES"].notna()]
##why df[df..
# only consider unique SMILES for now
df = df.drop_duplicates(subset=["SMILES"])
# calculate morgan fingerprints for each SMILES (some SMILES are invalid)
df["morgan"] = df["SMILES"].apply(fingerprint_from_smiles)
# remove rows without morgan fingerprints
df = df[~df["morgan"].isnull()]


similar_mol = {}

with Progress() as pg:
    task = pg.add_task("Processing initial matches...", total=len(df))
    for _, row in df.iterrows():
        #use smiles as key so can refer back to ki database before drop_duplicates
        smiles=row['SMILES']
        temp = compare_one_to_many(
            row["morgan"], res["bitvec"].to_list(), res["uids"].to_list()
        )
        if temp:
            similar_mol[smiles] = temp
        pg.update(task, advance=1)

#test with csvs
#npatlas_ori = os.path.join(
    #os.path.expanduser('~'), 'socialgene', 'addons', 'NPAtlas_download.tsv')
#pdsp_ki_ori = os.path.join(
    #os.path.expanduser('~'), 'socialgene', 'addons', 'KiDatabase.csv')

npatlas_ori=r"C:\Users\18918\Documents\Github\sgpy-main\socialgene\addons\NPAtlas_download.tsv"
pdsp_ki_ori=r"C:\Users\18918\Documents\Github\sgpy-main\socialgene\addons\KiDatabase.csv"

npatlas = pd.read_csv(npatlas_ori, sep='\t')
pdsp_ki = pd.read_csv(pdsp_ki_ori, sep=',', encoding_errors='replace')

pdsp_ki = tidy(pdsp_ki)
npatlas = tidy(npatlas, SMILES='compound_smiles')

similar_mol = {}

with Progress() as pg:
    task = pg.add_task("Processing initial matches...", total=len(pdsp_ki))
    for _, row in pdsp_ki.iterrows():
        #use smiles as key so can refer back to ki database before drop_duplicates
        smiles=row['SMILES']
        temp = compare_one_to_many(
            row["morgan"], npatlas["morgan"].to_list(), npatlas["npaid"].to_list()
        )
        if temp:
            similar_mol[smiles] = temp
        pg.update(task, advance=1)



#Parse Ki databse, add nodes for unique smiles

#Add edges with similar_mol