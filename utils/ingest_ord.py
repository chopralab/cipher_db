import argparse
from itertools import chain, islice
from pathlib import Path
import pprint
from typing import (
    Dict, Iterable, Iterator, List, Optional, Tuple, TypeVar, Union
)

import ray
from rdkit import Chem, RDLogger
from tqdm import tqdm

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

T = TypeVar('T')

_TEST_MAX_REACTIONS = 40000

RDLogger.DisableLog('rdApp.*') 

try:
    ray.init('auto')
except ConnectionError:
    ray.init()

def batches(it: Iterable[T], size: int) -> Iterator[List]:
    """Batch an iterable into batches of given size, with the final
    batch potentially being smaller"""
    it = iter(it)
    return iter(lambda: list(islice(it, size)), [])

@ray.remote
def _process_reaction_batch(
    reactions: Iterable[reaction_pb2.Reaction], reactants: bool = False,
) -> List[Tuple[str, List[str], List[str]]]:
    return [
        process_reaction(reaction, reactants) for reaction in reactions
    ]

def process_reaction(
    reaction: reaction_pb2.Reaction, reactants: bool = False
) -> Tuple[str, List[str], List[str]]:
    """process a single reaction to obtain the reaction ID, a list of products, 
    and a list of reactants, both as InChIKeys

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
        the InChIKeys of the products
    reactant_inchikeys : List[str]
        the InChIKeys of the reactants. Empty if reactants is False
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

def process_ord_file(filepath: Union[str, Path],
                     ingest_reactants: bool = False,
                     max_reactions: Optional[int] = None) -> List[Dict]:
    """process all the reactions in single ORD data file to obtain a list of
    table rows corresponding to each reaction

    Parameters
    ----------
    filepath : Path
        the filepath of an ORD data file
    ingest_reactants : bool, default=False
        whether or not to ingest reactants as well
    max_reactions: Optional[int], default=None
        the maximum number of reactions to process from the dataset. If None,
        process all reactions
        
    Returns
    -------
    rows : List[Dict]
        the table rows processed from the datset in the given file. Each row 
        corresponds to a single reaction in the dataset and is a dictionary
        of the form:

        * 'reaction_id': the ORD reaction ID of the reaction
        * 'products':  a list of reaction products as InChIKeys
        * 'reactants': a list of reaction reactants as InChIKeys. Empty if 
            ingest_reactants is False
    """
    dataset = message_helpers.load_message(str(filepath), dataset_pb2.Dataset)

    BATCH_SIZE = 512
    # NOTE(degraff): maybe a single for-loop is faster than 2 comprehensions
    refs = [
        _process_reaction_batch.remote(reaction_batch, ingest_reactants)
        for reaction_batch in batches(
            dataset.reactions[:max_reactions], BATCH_SIZE
        )
    ]
    rowss = [
        ray.get(ref) for ref in tqdm(
            refs, desc='Processing reaction batches', unit='batch', smoothing=0.
        )
    ]
    # rows = [
    #     process_reaction(reaction, ingest_reactants)
    #     for reaction in tqdm(dataset.reactions[:max_reactions], unit='rxn')
    # ]
    rows = [
        {'reaction_id': reaction_id,
         'products': products,
         'reactants': reactants}
        for reaction_id, products, reactants in chain(*rowss)
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
    parser.add_argument('--test', default=False, action='store_true',
                        help='whether to run in testing mode')
    parser.add_argument('--test-filenames', nargs='+',
                        help='specific filenames of ORD data files with which to test the script.')
    parser.add_argument('--max-rxns', type=int,
                        help='the maximum number of reactions to process from each dataset. By default, process all reactions from each dataset. If the --test flag is passed, proecss _TEST_MAX_REACTIONS reactions by default.')
    args = parser.parse_args()

    db_rows = []
    if args.test:
        max_rxns = args.max_rxns or _TEST_MAX_REACTIONS
        for filename in args.test_filenames:
            db_rows.extend(process_ord_file(
                filename, args.ingest_reactants, max_rxns
            ))
        
        # pprint.pprint(db_rows[:20], compact=True)
        exit(0)

    parent_dir = args.ord_data_repo / 'data'
    for data_dir in parent_dir.iterdir():
        for filepath in data_dir.iterdir():
            db_rows.extend(process_ord_file(
                filepath, args.ingest_reactants, args.max_rxns
            ))


if __name__ == '__main__':
    main()