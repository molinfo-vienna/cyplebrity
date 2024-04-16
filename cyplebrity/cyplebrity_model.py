import sys
from typing import List, Tuple

import numpy as np
import pandas as pd
from FPSim2 import FPSim2Engine
from joblib import load
from molvs import Standardizer, tautomer
from nerdd_module import AbstractModel
from rdkit import Chem
from rdkit.Chem import AllChem, Mol, MolToSmiles
from rdkit.rdBase import BlockLogs

if sys.version_info < (3, 9):
    from importlib_resources import files
else:
    from importlib.resources import files

__all__ = ["CyplebrityModel"]

labels = ["1A2", "2C9", "2C19", "2D6", "3A4"]


def load_ml_nn_models(labels):
    ml_data = []
    nn_data = []
    selectors = []

    for i in range(len(labels)):
        ml_data.append(
            load(
                files("cyplebrity")
                .joinpath("resources")
                .joinpath("mm_models")
                .joinpath("p450inh_final_" + labels[i] + "_model.pkl.gz")
            )
        )
        selectors.append(
            load(
                files("cyplebrity")
                .joinpath("resources")
                .joinpath("mm_models")
                .joinpath("p450inh_final_" + labels[i] + "_selector.pkl")
            )
        )
        nn_data.append(
            FPSim2Engine(
                str(
                    files("cyplebrity")
                    .joinpath("resources")
                    .joinpath("nn_models")
                    .joinpath("p450inh_final_" + labels[i] + ".h5")
                )
            )
        )

    return selectors, ml_data, nn_data


selectors, mlm, nnm = load_ml_nn_models(labels)


def get_nearest_neighbors_scores(model, mols, enable_ad):
    """
    Method to get the distance to closest training instance

    :param model: Machine learning model
    :param smiles: list of smiles
    :param enable_ad: False/True
    :return: nearest neighbor score of model for smiles
    """
    arr_size = len(mols)

    if enable_ad:
        scores = np.empty(arr_size)
        similarity = []

        for i in range(arr_size):
            similarity = model.similarity(MolToSmiles(mols[i]), 0.0, n_workers=1)
            if len(similarity) != 0:
                scores[i] = float(similarity[0][1])
            else:
                scores[i] = 0.0
        return scores
    else:
        return np.array([-1] * arr_size)


def get_machine_learning_prediction(model, selector, morgans):
    """
    Method to get the rounded predictions of the models

    :param model: Machine learning model
    :param morgans: list Morgan2FP for query
    :return: rounded predictions of model for morgans
    """
    if len(morgans) == 0:
        return []

    return [
        np.round(i, 2) for i in model.predict_proba(selector.transform(morgans))[:, 1]
    ]


def get_morgan2_fp(mol):
    """
    Method to calculate morgan2Fps

    :param mol: rdkit molecule
    :return: rdkits GetMorganFingerprintAsBitVect(mol,3,4096) (Morgan2FP)
    """
    return AllChem.GetMorganFingerprintAsBitVect(mol, 3, 4096)


def cleanAndCheckMolecule(mol):
    """
    Performs the structure validation and standardization in one step
    returns RDKit mol or None in case of any failure
    The tauromerization is controlled by global variable
    tautomerize = True

    In case of failure one of the following global counters is increased (to be used in dataset statistics)

    m_corrupt - caused any Python exception
    m_inorganic - does not contain any Carbon
    m_tooheavy - Molecular Weight over 1000
    m_badatoms - contains atoms other than ["C","N","O","S","P","F","Cl","Br","I", "B", "Si", "Se", "H"]
    m_exotic - radical

    Keyword arguments:
    mol -- RDKit mol structure
    """
    # Those are incremented on each call to cleanAndCheckMolecule if any problem occurs
    # could possibly be used to inform the user what was wrong with a given structure, if we want to
    m_corrupt = 0
    m_inorganic = 0
    m_tooheavy = 0
    m_badatoms = 0
    m_exotic = 0
    legal_atoms = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B", "Si", "Se", "H"]
    if not mol:
        return None, [m_corrupt, m_inorganic, m_tooheavy, m_badatoms, m_exotic]
    try:
        Chem.Cleanup(mol)
    except Exception as e:
        mol = None
        m_corrupt += 1

    if not mol:
        return None, [m_corrupt, m_inorganic, m_tooheavy, m_badatoms, m_exotic]
    try:
        mol = Chem.RemoveHs(mol)
    except Exception as e:
        m_corrupt += 1
        mol = None
    if not mol:
        return None, [m_corrupt, m_inorganic, m_tooheavy, m_badatoms, m_exotic]

    wr = True
    fragatoms = 0
    molfrags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    for frag in molfrags:
        if len(frag.GetAtoms()) > fragatoms:
            fragatoms = len(frag.GetAtoms())
            molOK = frag

    for atom in molOK.GetAtoms():
        isot = atom.GetIsotope()
        if not (isot == 0):
            m_exotic += 1
            return None, [m_corrupt, m_inorganic, m_tooheavy, m_badatoms, m_exotic]
        if not (atom.GetSymbol() in legal_atoms):
            m_badatoms += 1
            return None, [m_corrupt, m_inorganic, m_tooheavy, m_badatoms, m_exotic]
    if sum(1 for atom in molOK.GetAtoms() if atom.GetAtomicNum() == 6) == 0:
        m_inorganic += 1
        return None, [m_corrupt, m_inorganic, m_tooheavy, m_badatoms, m_exotic]

    if Chem.Descriptors.MolWt(Chem.AddHs(molOK)) > 1000:
        m_tooheavy += 1
        return None, [m_corrupt, m_inorganic, m_tooheavy, m_badatoms, m_exotic]
    if Chem.Descriptors.NumRadicalElectrons(molOK) > 0:
        m_exotic += 1
        return None, [m_corrupt, m_inorganic, m_tooheavy, m_badatoms, m_exotic]
    try:
        molOK = Chem.RemoveHs(molOK)
    except Exception as e:
        m_corrupt += 1
        molOK = None

    if not molOK:
        return None, [m_corrupt, m_inorganic, m_tooheavy, m_badatoms, m_exotic]

    try:
        Chem.SanitizeMol(molOK)
    except Exception as e:
        m_corrupt += 1
        molOK = None

    if not molOK:
        return None, [m_corrupt, m_inorganic, m_tooheavy, m_badatoms, m_exotic]

    s = Standardizer()
    molOK = s.standardize(molOK)
    canon = tautomer.TautomerCanonicalizer()
    try:
        molOK = canon.canonicalize(molOK)
    except:
        molOK = None

    if not molOK:
        m_corrupt += 1
        return None, [m_corrupt, m_inorganic, m_tooheavy, m_badatoms, m_exotic]

    try:
        AllChem.Compute2DCoords(molOK)
    except Exception as e:
        m_corrupt += 1
        molOK = None

    if not molOK:
        return None, [m_corrupt, m_inorganic, m_tooheavy, m_badatoms, m_exotic]

    props = molOK.GetPropNames(includePrivate=True, includeComputed=True)
    for prop in props:
        molOK.ClearProp(prop)
    molOK.SetProp("InChI", Chem.MolToInchi(molOK))
    molOK.SetProp("Smiles", Chem.MolToSmiles(molOK))
    if wr:
        return molOK, [m_corrupt, m_inorganic, m_tooheavy, m_badatoms, m_exotic]
    else:
        return None, [m_corrupt, m_inorganic, m_tooheavy, m_badatoms, m_exotic]


# TODO: delete this function since it is not used?
def produce_confidence_warnings(similarity):
    """
    Method that adds confidence warnings to the warnings

    :param args: w (hitdexter warnings; list), similarity (similarity to training set)
    :return: w with appended confidence warnings
    """
    warningsConfidence = ["NNPH", "NNPM", "NNCH", "NNCM"]
    for i, s in enumerate(similarity):
        conf = get_confidence_for_similarity(s)
        if conf != 1:
            w.append(warningsConfidence[i])
    return w


def produce_warnings(w):
    """
    Method that adds confidence warnings to the warnings

    :param args: w (hitdexter warnings; list), similarity (similarity to training set)
    :return: w with appended confidence warnings
    """
    ew = ""
    for i, err in enumerate(w):
        if err != 0:
            if i == 1:
                ew += "corrupt,"
            if i == 2:
                ew += "inorganic,"
            if i == 3:
                ew += "tooheavy,"
            if i == 4:
                ew += "badatoms,"
            if i == 5:
                ew += "exotic,"
        else:
            ew = "No errors detected "

    return ew[:-1]


def predict(
    mols,
    enable_ad: bool = True,
):
    # calculate features

    # prepare descriptors
    morgans = [get_morgan2_fp(m) for m in mols]

    # machine learning models
    mlm_predictions = [
        get_machine_learning_prediction(mlm[i], selectors[i], morgans)
        for i in range(len(labels))
    ]

    # calculate similarities to datasets
    nnm_predictions = [
        get_nearest_neighbors_scores(nnmodel, mols, enable_ad) for nnmodel in nnm
    ]

    distances = []

    if enable_ad:
        for similarity in nnm_predictions:
            distances.append([np.round(1 - s, 2) for s in similarity])
    else:
        distances = nnm_predictions

    results = pd.DataFrame()

    for i in range(len(labels)):
        results[f"prediction_{i+1}"] = mlm_predictions[i]
        results[f"neighbor_{i+1}"] = distances[i]

    return results


class CyplebrityModel(AbstractModel):
    def __init__(self):
        super().__init__(preprocessing_pipeline="custom")

    def _preprocess_single_mol(self, mol: Mol) -> Tuple[Mol, List[str]]:
        with BlockLogs():
            preprocessed_mol, flags = cleanAndCheckMolecule(mol)
            warnings = produce_warnings(flags)
            return preprocessed_mol, [warnings]

    def _predict_mols(
        self, mols: List[Mol], applicability_domain: bool = True
    ) -> pd.DataFrame:
        return predict(mols, applicability_domain)
