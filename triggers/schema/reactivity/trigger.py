from importlib import import_module
import os
import tempfile
from typing import Dict

import mongoengine as me
import pymongo as pmg
from triggers.schema.reactivity import viz

from triggers.schema.reactivity.client import AskcosClient
from triggers.schema.reactivity.docs import (
    ChemicalNode,
    Difficulty,
    ReactionNode,
    Retrosynthesis,
    SyntheticTree,
)
from triggers.schema.reactivity.sascorer import SAScorer

me.connect("cipher_aspire", host=os.environ["MONGO_URI"])
MONGO_CLIENT = pmg.MongoClient(os.environ["MONGO_URI"])
ASKCOS_CLIENT = AskcosClient(os.environ["ASKCOS_HOST"])
SA_SCORER = SAScorer(os.environ["FP_SCORES_PKL"])


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
    cn.children = [build_rxn_node(child) for child in tree["children"]]

    return cn


def build_synthetic_tree(tree: Dict) -> SyntheticTree:
    st = SyntheticTree()

    st.depth = tree["attributes"]["depth"]
    st.precursor_cost = tree["attributes"]["precursor_cost"]
    st.score = tree["attributes"]["score"]
    st.cluster_id = tree["attributes"]["cluster_id"]
    st.root = build_chemical_node(tree)

    with tempfile.NamedTemporaryFile() as fid:
        g = viz.tree_to_graph(viz.clean_tree(tree))
        g.render(outfile=fid.name, format="png", cleanup=True)
        fid.seek(0)
        st.image = fid

    return st


def update_retrosynthesis(inchikey, smi):
    retro = Retrosynthesis()

    retro.inchikey = inchikey
    retro.smi = smi
    retro.trees = [build_synthetic_tree(tree) for tree in ASKCOS_CLIENT.get_trees(smi)]


def update_difficulty(inchikey, smi):
    difficulty = Difficulty()

    difficulty.inchikey = inchikey
    difficulty.smiles = smi
    difficulty.sc_score = ASKCOS_CLIENT.sc_score(smi)
    difficulty.sa_score = SA_SCORER(smi)

    difficulty.save()

    return difficulty


def difficulty_trigger():
    feasibility = MONGO_CLIENT.cipher_aspire.feasibility
    try:
        with feasibility.watch([{"$match": {"operationType": "insert"}}]) as stream:
            for change in stream:
                pass
    except pmg.errors.PyMongoError as e:
        print(e)
