import sys
from typing import Iterator, List, Optional, Tuple

import numpy as np
from FPSim2 import FPSim2Engine
from joblib import load
from molvs import Standardizer, tautomer
from nerdd_module import InvalidElementsProblem, InvalidWeightProblem, Model, Problem
from nerdd_module.polyfills import BlockLogs
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Mol, MolToSmiles

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


def _append_corrupt_molecule_problem(
    problems: List[Problem], m_corrupt
) -> List[Problem]:
    if m_corrupt > 0:
        problems.append(Problem("invalid_molecule", "Molecule could not be processed."))
    return problems


def cleanAndCheckMolecule(mol) -> Tuple[Optional[Mol], List[Problem]]:
    """
    Performs the structure validation and standardization in one step
    returns RDKit mol or None in case of any failure
    The tauromerization is controlled by global variable
    tautomerize = True

    Keyword arguments:
    mol -- RDKit mol structure
    """
    problems = []

    # Those are incremented on each call to cleanAndCheckMolecule if any problem occurs
    # could possibly be used to inform the user what was wrong with a given structure, if we want to
    m_corrupt = 0

    legal_atoms = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B", "Si", "Se", "H"]
    if not mol:
        return None, _append_corrupt_molecule_problem(problems, m_corrupt)
    try:
        Chem.Cleanup(mol)
    except Exception as e:
        mol = None
        m_corrupt += 1

    if not mol:
        return None, _append_corrupt_molecule_problem(problems, m_corrupt)
    try:
        mol = Chem.RemoveHs(mol)
    except Exception as e:
        m_corrupt += 1
        mol = None
    if not mol:
        return None, _append_corrupt_molecule_problem(problems, m_corrupt)

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
            problems.append(
                Problem(
                    "invalid_molecule",
                    f"Invalid isotope {isot} for atom {atom.GetSymbol()}",
                )
            )
            return None, _append_corrupt_molecule_problem(problems, m_corrupt)
        if atom.GetSymbol() not in legal_atoms:
            problems.append(InvalidElementsProblem([atom.GetSymbol()]))
            return None, _append_corrupt_molecule_problem(problems, m_corrupt)
    if sum(1 for atom in molOK.GetAtoms() if atom.GetAtomicNum() == 6) == 0:
        problems.append(
            Problem("invalid_molecule", "Molecule does not contain any carbon atoms.")
        )
        return None, _append_corrupt_molecule_problem(problems, m_corrupt)

    mol_wt = Descriptors.MolWt(Chem.AddHs(molOK))
    if mol_wt > 1000:
        problems.append(InvalidWeightProblem(mol_wt, 0, 1000))
        return None, _append_corrupt_molecule_problem(problems, m_corrupt)
    if Descriptors.NumRadicalElectrons(molOK) > 0:
        problems.append(
            Problem("invalid_molecule", "Molecule contains radical electrons.")
        )
        return None, _append_corrupt_molecule_problem(problems, m_corrupt)
    try:
        molOK = Chem.RemoveHs(molOK)
    except Exception as e:
        m_corrupt += 1
        molOK = None

    if not molOK:
        return None, _append_corrupt_molecule_problem(problems, m_corrupt)

    try:
        Chem.SanitizeMol(molOK)
    except Exception as e:
        m_corrupt += 1
        molOK = None

    if not molOK:
        return None, problems

    s = Standardizer()
    molOK = s.standardize(molOK)
    canon = tautomer.TautomerCanonicalizer()
    try:
        molOK = canon.canonicalize(molOK)
    except:
        molOK = None

    if not molOK:
        m_corrupt += 1
        return None, problems

    try:
        AllChem.Compute2DCoords(molOK)
    except Exception as e:
        m_corrupt += 1
        molOK = None

    if not molOK:
        return None, problems

    props = molOK.GetPropNames(includePrivate=True, includeComputed=True)
    for prop in props:
        molOK.ClearProp(prop)
    molOK.SetProp("InChI", Chem.MolToInchi(molOK))
    molOK.SetProp("Smiles", Chem.MolToSmiles(molOK))
    if wr:
        return molOK, problems
    else:
        return None, problems


def predict(
    mols,
    enable_ad: bool = True,
) -> Iterator[dict]:
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

    for j in range(len(mols)):
        yield {
            **{
                f"prediction_{i + 1}": mlm_predictions[i][j] for i in range(len(labels))
            },
            **{f"neighbor_{i + 1}": distances[i][j] for i in range(len(labels))},
        }


class CyplebrityModel(Model):
    def __init__(self):
        super().__init__(preprocessing_steps=[])

    def _preprocess(self, mol: Mol) -> Tuple[Optional[Mol], List[Problem]]:
        with BlockLogs():
            preprocessed_mol, problems = cleanAndCheckMolecule(mol)

        return preprocessed_mol, problems

    def _predict_mols(
        self, mols: List[Mol], applicability_domain: bool = True
    ) -> Iterator[dict]:
        return predict(mols, applicability_domain)
