import pytest

from cipher_reactivity import viz

@pytest.fixture
def tree():
    return {
        "is_chemical": True,
        "as_reactant": 0,
        "as_product": 2,
        "id": 9,
        "smiles": "CCN(CC)CCOC(C)(c1ccccc1)c1ccc(Cl)cc1",
        "ppg": 0.0,
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
                                "smiles": "C[Mg]I.O=C(c1ccccc1)c1ccc(Cl)cc1>>CC(O)(c1ccccc1)c1ccc(Cl)cc1",  # noqa: E501
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
                "smiles": "CC(O)(c1ccccc1)c1ccc(Cl)cc1.CCN(CC)CCBr>>CCN(CC)CCOC(C)(c1ccccc1)c1ccc(Cl)cc1",  # noqa: E501
            }
        ],
    }

def walk_tree(cleaned_tree):
    VALID_KEYS = {"smiles", "children"}

    for key, val in cleaned_tree.items():
        assert key in VALID_KEYS
        if key == "children":
            for subtree in val:
                walk_tree(subtree)

def test_clean_tree(tree):
    cleaned_tree = viz.clean_tree(tree)
    walk_tree(cleaned_tree)

def test_clean_tree_integrity(tree):
    cleaned_tree = viz.clean_tree(tree)
    recleaned_tree = viz.clean_tree(tree)

    assert cleaned_tree == recleaned_tree