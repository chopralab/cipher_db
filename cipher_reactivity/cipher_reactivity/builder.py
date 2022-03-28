from importlib import resources
import json
import os
import tempfile
from typing import Dict
import warnings
from cipher_reactivity.client import AskcosClient
from cipher_reactivity.sascorer import SAScorer

import graphviz

from cipher_reactivity import viz
from cipher_reactivity.docs import (
    ChemicalNode,
    Difficulty,
    ReactionNode,
    Retrosynthesis,
    SyntheticTree,
)
import tomli


ASKCOST_HOST = os.environ["ASKCOS_HOST"]
try:
    TREE_PARAMS = tomli.loads(resources.read_text("cipher_reactivity.data", "tree_params.toml"))
except ValueError:
    TREE_PARAMS = None
ASKCOS_CLIENT = AskcosClient(ASKCOST_HOST, TREE_PARAMS)

SA_SCORER = SAScorer(json.loads(resources.read_text("cipher_reactivity.data", "fpscores.json")))


def build_rxn_node(tree: Dict) -> ReactionNode:
    rn = ReactionNode()

    rn.smiles = tree["smiles"]
    rn.tforms = tree["tforms"]
    rn.tsources = tree["tsources"]
    rn.template_score = tree["template_score"]
    rn.plausiblity = tree["plausibility"]
    rn.rank = tree["rank"]
    rn.num_examples = tree["num_examples"]
    rn.necessary_reagent = tree["necessary_reagent"]
    rn.precursor_smiles = tree["precursor_smiles"]
    rn.rms_molwt = tree["rms_molwt"]
    rn.num_rings = tree["num_rings"]
    rn.scscore = tree["scscore"]
    rn.rxn_id = tree["id"]
    rn.forward_score = tree["forward_score"]
    rn.class_num = tree["class_num"]
    rn.class_name = tree["class_name"]
    rn.children = [build_chemical_node(child) for child in tree["children"]]

    return rn


def build_chemical_node(tree: Dict) -> ChemicalNode:
    cn = ChemicalNode()

    cn.smiles = tree["smiles"]
    cn.ppg = tree["ppg"]
    cn.as_product = tree["as_product"]
    cn.as_reactant = tree["as_reactant"]
    cn.terminal = tree["terminal"]
    cn.chemical_id = tree["id"]
    if not cn.terminal:
        cn.rxn = build_rxn_node(tree["children"][0])

    return cn


def build_synthetic_tree(tree: Dict) -> SyntheticTree:
    st = SyntheticTree()

    st.depth = tree["attributes"]["depth"]
    st.precursor_cost = tree["attributes"]["precursor_cost"]
    st.score = tree["attributes"]["score"]
    st.cluster_id = tree["attributes"]["cluster_id"]
    st.root = build_chemical_node(tree)

    try:
        with tempfile.NamedTemporaryFile() as fid, warnings.catch_warnings():
            warnings.simplefilter("ignore", category=graphviz.exceptions.UnknownSuffixWarning)
            g = viz.tree_to_graph(viz.clean_tree(tree))
            g.render(outfile=fid.name, format="png", cleanup=True)
            fid.seek(0)
            st.image = fid
    except (IOError, graphviz.backend.execute.CalledProcessError) as e:
        print(e)
    return st


def build_retrosynthesis(inchikey, smiles):
    retro = Retrosynthesis()

    retro.inchikey = inchikey
    retro.smiles = smiles

    trees = ASKCOS_CLIENT.get_trees(smiles)
    if trees:
        retro.trees = [build_synthetic_tree(t) for t in trees]

    retro.save()

    return retro


def build_difficulty(inchikey, smiles):
    difficulty = Difficulty()

    difficulty.inchikey = inchikey
    difficulty.smiles = smiles
    difficulty.sc_score = ASKCOS_CLIENT.sc_score(smiles)
    difficulty.sa_score = SA_SCORER(smiles)

    difficulty.save()

    return difficulty
