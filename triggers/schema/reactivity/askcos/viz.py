from pathlib import Path
import tempfile
from typing import Dict, Iterable, List, Optional, Set, Tuple

import graphviz
from rdkit import Chem
from rdkit.Chem import Draw

NODE_ATTRS = dict(label="", shape="rect", style="square", weight="4", penwidth="2")
EDGE_ATTRS = dict(penwidth="2")

TMP_DIR = Path(tempfile.gettempdir()) / "tree_builder_assets"
TMP_DIR.mkdir(exist_ok=True, parents=True)

DRAW_OPTIONS = Draw.rdMolDraw2D.MolDrawOptions()
DRAW_OPTIONS.bondLineWidth = 6
DRAW_OPTIONS.fixedBondLength = 20
DRAW_OPTIONS.minFontSize = 16
DRAW_OPTIONS.useBWAtomPalette()
DRAW_OPTIONS.scaleBondWidth = True


def clean_tree(tree: Dict) -> Optional[Dict]:
    """clean a retroysynthetic tree json into only smiles and children. None if cleaning failed"""
    if "children" in tree and "is_reaction" in tree:
        return [clean_tree(child) for child in tree["children"]]

    elif "children" in tree and "smiles" in tree:
        if len(tree["children"]) == 1:
            children = clean_tree(tree["children"][0])
            return {"smiles": tree["smiles"], "children": children}
        elif len(tree["children"]) == 0:
            return {"smiles": tree["smiles"]}
        else:
            print("Error with format")
            return None
    else:
        print(tree)
        print("Error with format...")
        return None


def get_edges(tree: Dict, parent: Optional[str] = None):
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


def get_leaves(edges: List[Tuple], nodes: Iterable[str]) -> Set:
    """given a list of edges, finds nodes with no children"""
    children = set(nodes)
    for e in edges:
        parent = e[0]
        if parent in children:
            children.remove(parent)

    return children


def build_graph(root: str, edges: List[Tuple], splines="none"):
    """build a graph from the edges of a retrosynthetic tree"""
    parents, children = zip(*edges)
    d_child_parent = dict(zip(children, parents))
    parents = set(parents)

    unique_smis = {*parents, *children}
    leaves = get_leaves(edges, unique_smis)

    d_smi_img = {}
    for smi in unique_smis:
        mol = Chem.MolFromSmiles(smi)
        img_path = str(TMP_DIR / f'{smi.replace("/", "%2F")}.png')
        Draw.MolToFile(mol, img_path, (250, 250), options=DRAW_OPTIONS)
        d_smi_img[smi] = img_path

    graph = graphviz.Digraph(
        format="png", graph_attr={"splines": splines}, node_attr=NODE_ATTRS, edge_attr=EDGE_ATTRS
    )

    for smi in parents:
        graph.node(smi, image=d_smi_img[smi], color="darkred" if smi == root else "darkblue")
        graph.node(f"{smi}-branch", label="", width="0", height="0")
        graph.edge(smi, f"{smi}-branch", arrowhead="none")

    for smi in children:
        if smi not in parents:
            graph.node(
                smi, image=d_smi_img[smi], color="darkgreen" if smi in leaves else "darkblue"
            )
        parent = d_child_parent[smi]
        graph.edge(f"{parent}-branch", smi)

    return d_smi_img, graph


def tree_to_image(tree, path):
    root = tree["smiles"]
    edges = get_edges(tree)

    d_smi_img, graph = build_graph(root, edges, "ortho")
    graph.render(outfile=path, cleanup=True)

    # for img_path in d_smi_img.values():
    #     os.remove(img_path)

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

    graph = tree_to_image(clean_tree(tree), "test.png")
    # print(graph)