from rdkit import RDLogger

lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

import os
from glob import glob

import pandas as pd
from FPSim2.io import create_db_file
from rdkit.Chem import MolToSmiles

from ..cyplebrity_model import CyplebrityModel

# Note: we use the model although we are computing similarity databases for exactly
# this model. However, we don't use the model's computed distances during the
# following computation.
m = CyplebrityModel()


def create_db(filepath):
    db_path = f"./{os.path.basename(filepath)}.h5"
    if os.path.exists(db_path):
        return

    # load data
    df = pd.read_csv(filepath, sep="\t")
    smiles_list = df.Smiles.tolist()

    df_predictions = m.predict(smiles_list)
    preprocessed_smiles_list = (
        df_predictions.preprocessed_mol.pipe(lambda s: s[pd.notnull(s)])
        .map(MolToSmiles)
        .pipe(lambda s: s[pd.notnull(s)])
    )

    # create fpsim db
    create_db_file(
        [(smiles, i) for i, smiles in enumerate(preprocessed_smiles_list)],
        db_path,
        "Morgan",
        {"radius": 3, "nBits": 4096},
    )


filepaths = glob("../resources/training_data/*.txt")
for filepath in filepaths:
    create_db(filepath)
