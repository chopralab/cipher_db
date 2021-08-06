from rdkit import Chem
from rdkit.Chem import Draw
import graphviz
import os

def clean_tree(tree):
    '''
    Clean a retroysynthetic tree json into only smiles and children
    '''
    if 'children' in tree and 'is_reaction' in tree:
        children = []
        for child in tree['children']:
            children.append(clean_tree(child))
        return children
    elif 'children' in tree and 'smiles' in tree:
        if len(tree['children']) == 1:
            children = clean_tree(tree['children'][0])
            return {'smiles': tree['smiles'], 'children': children}
        elif len(tree['children']) == 0:
            return {'smiles': tree['smiles']}
        else:
            print('Error with format')
            return -1
    else:
        print(tree)
        print('Error with format...')
        return -1

def get_edges(tree, parent=None):
    '''
    Given a tree in the cleaned format, returns a list of edges
    '''
    smi = tree['smiles']
    edges = []
    if parent is not None:
        edges = [(parent, smi)]
    if 'children' in tree:
        for child in tree['children']:
            edges += get_edges(child, parent=smi)
        return edges
    else:
        return edges

def gen_viz(edges, splines='none'):
    '''
    Given a list of edges of a retrosynthetic tree, generates a graph visualization
    '''
    unismi = set([])
    img_map = {}
    for edge in edges:
        unismi.add(edge[0])
        unismi.add(edge[1])
    for i in unismi:
        mol = Chem.MolFromSmiles(i)
        j = i.replace('/', '%2F')
        Draw.MolToFile(mol, j + '.png')
        img_map[i] = j + '.png'

    if splines == 'ortho':
        graph = graphviz.Digraph(format='png', graph_attr={'splines': splines})
    else:
        graph = graphviz.Digraph(format='png')
    parents = []
    for edge in edges:
        graph.node(edge[0], label='', image=img_map[edge[0]], shape='rect', style='rounded')
        graph.node(edge[1], label='', image=img_map[edge[1]], shape='rect', style='rounded')
        graph.node(f'd{edge[0]}', label='', width='0', height='0')
        if edge[0] not in parents:
            parents.append(edge[0])
            graph.edge(edge[0], f'd{edge[0]}', arrowhead='none')
        graph.edge(f'd{edge[0]}', edge[1])
    
    return img_map, graph

class TreeVisualizer():

    def __init__(self):
        pass

    def tree_to_file(self, tree, path):
        tree = clean_tree(tree)
        edges = get_edges(tree)
        try:
            img_map, graph = gen_viz(edges, 'ortho')
            graph.render(path)
        except:
            img_map, graph = gen_viz(edges)
            graph.render(path)
        for key,val in img_map.items():
            os.remove(val)
        return graph


if __name__ == '__main__':
    tree = {
            "is_chemical": True,
            "children": [
                {
                "tforms": [
                    "5e1f4b6e6348832850995fbf",
                    "5e1f4b6e6348832850995f00"
                ],
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
                        "tforms": [
                            "5e1f4b6e6348832850998838"
                        ],
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
                            "ppg": 1.0
                            },
                            {
                            "is_chemical": True,
                            "children": [],
                            "as_reactant": 1036,
                            "as_product": 560,
                            "id": 2,
                            "smiles": "O=C(c1ccccc1)c1ccc(Cl)cc1",
                            "ppg": 1.0
                            }
                        ],
                        "smiles": "C[Mg]I.O=C(c1ccccc1)c1ccc(Cl)cc1>>CC(O)(c1ccccc1)c1ccc(Cl)cc1"
                        }
                    ],
                    "as_reactant": 8,
                    "as_product": 32,
                    "id": 6,
                    "smiles": "CC(O)(c1ccccc1)c1ccc(Cl)cc1",
                    "ppg": 0.0
                    },
                    {
                    "is_chemical": True,
                    "children": [],
                    "as_reactant": 372,
                    "as_product": 34,
                    "id": 12,
                    "smiles": "CCN(CC)CCBr",
                    "ppg": 5.0
                    }
                ],
                "smiles": "CC(O)(c1ccccc1)c1ccc(Cl)cc1.CCN(CC)CCBr>>CCN(CC)CCOC(C)(c1ccccc1)c1ccc(Cl)cc1"
                }
            ],
            "as_reactant": 0,
            "as_product": 2,
            "id": 9,
            "smiles": "CCN(CC)CCOC(C)(c1ccccc1)c1ccc(Cl)cc1",
            "ppg": 0.0
            }
   
    tv = TreeVisualizer()
    graph = tv.tree_to_file(tree, 'test')