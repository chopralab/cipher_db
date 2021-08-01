import argparse
from pathlib import Path
from typing import Dict, List, Tuple

from rdkit import Chem

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

def process_reaction(
    reaction: reaction_pb2.Reaction, reactants: bool = False
) -> Tuple[str, List[str], List[str]]:
    """
    process a single reaction to obtain the products and reactants

    Parameters
    ----------
    reaction : reaction_pb2.Reaction
        the reaction to process
    reactants : bool, default=False
        whether to also process the reaction for the reactants

    Returns
    -------
    reaction_id : str
        the ORD reaction ID of the given reaction
    product_inchikeys : List[str]
        the inchikeys of the products
    reactant_inchikeys : List[str]
        the inchikeys of the reactants. Empty if reactants is False
    """
    product_mols = []
    reactant_mols = []

    for outcome in reaction.outcomes:
        for product in outcome.products:
            try:
                product_mols.append(
                    message_helpers.mol_from_compound(product)
                )
            except ValueError:
                continue
    
    if reactants:
        for reaction_input in reaction.inputs.values():
            for compound in reaction_input.components:
                # NOTE(degraff): should we be filtering solely for reagents?
                if compound.reaction_role != reaction_pb2.ReactionRole.REACTANT:
                    continue

                try:
                    reactant_mols.append(
                        message_helpers.mol_from_compound(compound)
                    )
                except ValueError:
                    continue

    product_inchikeys = [Chem.MolToInchiKey(mol) for mol in product_mols]
    reactant_inchikeys = [Chem.MolToInchiKey(mol) for mol in reactant_mols]

    return reaction.reaction_id, product_inchikeys, reactant_inchikeys

def process_file(filepath: Path, ingest_reactants: bool = False) -> List[Dict]:
    """process all the reactions in single ORD data file to obtain a list of
    table rows corresponding to each reaction

    Parameters
    ----------
    filepath : Path
        the filepath of the data file
    ingest_reactants : bool, default=False
        whether or not to ingest reactants as well

    Returns
    -------
    rows : List[Dict]
        the table rows processed from the datset in the given file. Each row 
        corresponds to a single reaction in the dataset and is a dictionary
        of the form:

        * 'reaction_id': the ORD reaction ID of the reaction
        * 'products':  a list of reaction products as inchikeys
        * 'reactants': a list of reaction reactants as inchikeys. Empty if 
            ingest_reactants is False
    """
    dataset = message_helpers.load_message(
        str(filepath), dataset_pb2.Dataset
    )
    # NOTE(degraff): maybe a single for-loop is faster than 2 comprehensions
    rows = [
        process_reaction(reaction, ingest_reactants)
        for reaction in dataset.reactions
    ]
    rows = [
        {'reaction_id': reaction_id,
         'products': products,
         'reactants': reactants}
        for reaction_id, products, reactants in rows
    ]
    return rows
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d','--ord-data-repo', type=Path,
                        default=Path.home() / 'ord/ord-data',
                        help='the path of the ORD-data repository')
    parser.add_argument('--ingest-reactants', default=False, 
                        action='store_true',
                        help='whether to ingest reactants as well as products')
    args = parser.parse_args()

    parent_dir = args.ord_data_repo / 'data'

    db_rows = []
    for data_dir in parent_dir.iterdir():
        for filepath in data_dir.iterdir():
            db_rows.extend(process_file(filepath, args.ingest_reactants))


if __name__ == '__main__':
    main()