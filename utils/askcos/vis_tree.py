import json
from pprint import pprint
import os
from typing import Set

import graphviz
from rdkit import Chem
from rdkit.Chem import Draw


def clean_tree(tree):
    """clean a retroysynthetic tree json into only smiles and children"""
    if "children" in tree and "is_reaction" in tree:
        children = []
        for child in tree["children"]:
            children.append(clean_tree(child))
        return children
    elif "children" in tree and "smiles" in tree:
        if len(tree["children"]) == 1:
            children = clean_tree(tree["children"][0])
            return {"smiles": tree["smiles"], "children": children}
        elif len(tree["children"]) == 0:
            return {"smiles": tree["smiles"]}
        else:
            print("Error with format")
            return -1
    else:
        print(tree)
        print("Error with format...")
        return -1


def get_edges(tree, parent=None):
    """given a tree in the cleaned format, returns a list of edges"""
    smi = tree["smiles"]
    edges = []

    if parent is not None:
        edges = [(parent, smi)]
    if "children" in tree:
        for child in tree["children"]:
            edges += get_edges(child, parent=smi)

        return edges
    else:
        return edges


def get_leaves(edges, nodes) -> Set:
    """given a list of edges, finds nodes with no children"""
    children = set(nodes)
    for e in edges:
        parent = e[0]
        if parent in children:
            children.remove(parent)

    return children


def gen_viz(root: str, edges, splines="none"):
    """generate a graph visualization from the edges of a retrosynthetic tree"""
    unismis = set()
    d_smi_img = {}

    for e in edges:
        unismis.add(e[0])
        unismis.add(e[1])

    for smi in unismis:
        mol = Chem.MolFromSmiles(smi)
        filename = smi.replace("/", "%2F")
        Draw.MolToFile(mol, filename + ".png")
        d_smi_img[smi] = filename + ".png"

    if splines == "ortho":
        graph = graphviz.Digraph(format="png", graph_attr={"splines": splines})
    else:
        graph = graphviz.Digraph(format="png")

    parents = set()
    leaves = get_leaves(edges, unismis)
    pprint(edges, compact=False)
    print(unismis)
    print(leaves)
    for e in edges:
        # print(e)
        graph.node(
            e[0],
            label="",
            image=d_smi_img[e[0]],
            shape="rect",
            style="rounded",
            color="darkred" if e[0] == root else "darkblue",
        )
        graph.node(f"d{e[0]}", label="", width="0", height="0")

        if e[0] not in parents:
            parents.add(e[0])
            graph.edge(e[0], f"d{e[0]}", arrowhead="none")

        if e[1] in leaves:
            graph.node(
                e[1],
                label="",
                image=d_smi_img[e[1]],
                shape="rect",
                style="rounded",
                color="darkgreen",
            )
        else:
            graph.node(
                e[1],
                label="",
                image=d_smi_img[e[1]],
                shape="rect",
                style="rounded",
                color="darkblue",
            )
        graph.edge(f"d{e[0]}", e[1])



    return d_smi_img, graph


def traverse_tree(tree):
    graph = graphviz.Digraph(format="png")
    for k, v in tree.items():
        if k == "smiles":
            pass


def tree_to_image(tree, path):
    tree = clean_tree(tree)
    edges = get_edges(tree)

    # print(json.dumps(tree, indent=4))
    try:
        d_smi_img, graph = gen_viz(tree["smiles"], edges, "ortho")
        graph.render(path)
    except:
        d_smi_img, graph = gen_viz(tree["smiles"], edges)
        graph.render(path)

    for v in d_smi_img.values():
        os.remove(v)

    return graph


if __name__ == "__main__":
    tree = {
        "is_chemical": True,
        "children": [
            {
                "tforms": ["5e1f4b6e6348832850995fbf", "5e1f4b6e6348832850995f00"],
                "necessary_reagent": "",
                "id": 13,
                "is_reaction": True,
                "template_score": 0.017360973591475988,
                "num_examples": 3287,
                "plausibility": 0.992301,
                "children": [
                    {
                        "is_chemical": True,
                        "children": [
                            {
                                "tforms": ["5e1f4b6e6348832850998838"],
                                "necessary_reagent": "",
                                "id": 3,
                                "is_reaction": True,
                                "template_score": 0.06363593040428359,
                                "num_examples": 93,
                                "plausibility": 0.999670863,
                                "children": [
                                    {
                                        "is_chemical": True,
                                        "children": [],
                                        "as_reactant": 15576,
                                        "as_product": 1320,
                                        "id": 1,
                                        "smiles": "C[Mg]I",
                                        "ppg": 1.0,
                                    },
                                    {
                                        "is_chemical": True,
                                        "children": [],
                                        "as_reactant": 1036,
                                        "as_product": 560,
                                        "id": 2,
                                        "smiles": "O=C(c1ccccc1)c1ccc(Cl)cc1",
                                        "ppg": 1.0,
                                    },
                                ],
                                "smiles": "C[Mg]I.O=C(c1ccccc1)c1ccc(Cl)cc1>>CC(O)(c1ccccc1)c1ccc(Cl)cc1",
                            }
                        ],
                        "as_reactant": 8,
                        "as_product": 32,
                        "id": 6,
                        "smiles": "CC(O)(c1ccccc1)c1ccc(Cl)cc1",
                        "ppg": 0.0,
                    },
                    {
                        "is_chemical": True,
                        "children": [],
                        "as_reactant": 372,
                        "as_product": 34,
                        "id": 12,
                        "smiles": "CCN(CC)CCBr",
                        "ppg": 5.0,
                    },
                ],
                "smiles": "CC(O)(c1ccccc1)c1ccc(Cl)cc1.CCN(CC)CCBr>>CCN(CC)CCOC(C)(c1ccccc1)c1ccc(Cl)cc1",
            }
        ],
        "as_reactant": 0,
        "as_product": 2,
        "id": 9,
        "smiles": "CCN(CC)CCOC(C)(c1ccccc1)c1ccc(Cl)cc1",
        "ppg": 0.0,
    }

    graph = tree_to_image(tree, "test")
    # print(graph)
