import argparse
from itertools import chain, islice
from pathlib import Path
import pprint
from typing import (
    Dict, Iterable, Iterator, List, Optional, Tuple, TypeVar, Union
)

import ray
from rdkit import Chem, RDLogger
import pymongo
from tqdm import tqdm

from ord_schema import message_helpers
from ord_schema.proto import dataset_pb2
from ord_schema.proto import reaction_pb2

T = TypeVar('T')
_TEST_MAX_REACTIONS = 5000
RDLogger.DisableLog('rdApp.*') 

try:
    ray.init('auto')
except ConnectionError:
    ray.init()

def get_host(username: str, password: str) -> str:
    #NOTE(degraff): there's gotta be a better way to do this
    return f'mongodb+srv://{username}:{password}@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'

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

def process_ord_file(collection: pymongo.collection.Collection,
                     filepath: Union[str, Path],
                     ingest_reactants: bool = False,
                     max_reactions: Optional[int] = None):
    """process all the reactions in single ORD data file to obtain a list of
    documents corresponding to each reaction and insert them to the collection

    Each document is of the following form:

        * 'reaction_id': the ORD reaction ID of the reaction
        * 'products':  a list of reaction products as InChIKeys
        * 'reactants': a list of reaction reactants as InChIKeys. Empty if 
            ingest_reactants is False

    Parameters
    ----------
    collection : pymongo.collection.Collection
        the collection to which data should be uploaded
    filepath : Union[str, Path]
        the filepath of the ORD data file
    ingest_reactants : bool, default=False
        whether or not to ingest reactants as well
    max_reactions: Optional[int], default=None
        the maximum number of reactions to process from the dataset. If None,
        process all reactions
    """
    dataset = message_helpers.load_message(str(filepath), dataset_pb2.Dataset)
    max_reactions = max_reactions or len(dataset.reactions)

    BATCH_SIZE = 8192
    refs = [
        _process_reaction_batch.remote(reaction_batch, ingest_reactants)
        for reaction_batch in batches(
            dataset.reactions[:max_reactions], BATCH_SIZE
        )
    ]
    for ref in tqdm(refs, desc='Processing reactions', unit_scale=BATCH_SIZE,
                    unit='reaction', leave=False):
        docs = [
            {'reaction_id': reaction_id,
             'products': products,
             'reactants': reactants}
            for reaction_id, products, reactants in ray.get(ref)
        ]
        collection.insert_many(docs)
    # rowss = [
    #     ray.get(ref) for ref in tqdm(
    #         refs, desc='Processing reaction batches', unit='batch', smoothing=0.
    #     )
    # ]
    # # rows = [
    # #     process_reaction(reaction, ingest_reactants)
    # #     for reaction in tqdm(dataset.reactions[:max_reactions], unit='rxn')
    # # ]
    # rows = [
    #     {'reaction_id': reaction_id,
    #      'products': products,
    #      'reactants': reactants}
    #     for reaction_id, products, reactants in chain(*rowss)
    # ]
    # return rows
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-u', '--username',
                        help='the MongoDB username')
    parser.add_argument('-p', '--password',
                        help='the associated MongoDB password')
    parser.add_argument('-d','--ord-data-repo', type=Path, required=True,
                        help='the path of the ORD-data repository')
    parser.add_argument('--ingest-reactants', default=False, 
                        action='store_true',
                        help='whether to ingest reactants as well as products')
    parser.add_argument('--test', default=False, action='store_true',
                        help='whether to run in testing mode')
    parser.add_argument('--test-filenames', nargs='+',
                        help='specific filenames of ORD data files with which to test the script.')
    parser.add_argument('--max-rxns', type=int,
                        help='the maximum number of reactions to process from each dataset. By default, process all reactions from each dataset. If the --test flag is passed, process _TEST_MAX_REACTIONS reactions by default.')

    args = parser.parse_args()

    host = get_host(args.username, args.password)
    client = pymongo.MongoClient(host)
    db = client.cipher_aspire

    if args.test:
        max_rxns = args.max_rxns or _TEST_MAX_REACTIONS
        for filename in args.test_filenames:
            process_ord_file(
                db.ord, filename, args.ingest_reactants, max_rxns
            )
        client.close()
    else:
        data_dir = args.ord_data_repo / 'data'
        for filepath in tqdm(data_dir.glob('*/*.pb.gz'),
                            desc='Processing files', unit='file'):
            process_ord_file(
                db.ord, filepath, args.ingest_reactants, args.max_rxns
            )

    client.close()
    exit()

if __name__ == '__main__':
    main()