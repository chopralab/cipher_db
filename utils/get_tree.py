import requests
from pprint import pprint

class TreeBuilder():
    
    def __init__(self, host: str):
        '''
        Initialize a tree building object

        Parameters:
            host: IP address of server where ASKCOS is deployed
        '''
        self.host = host

    def build_tree(self, smi: str, params: dict):
        '''
        Submits a tree building request for the given smiles and parameters, returning the json response

        Parameters:
            smi: SMILES string
            params: dictionary containing the parameters of the query

        Returns:
            resp.json(): synthetic tree or trees found for the given request in json format
        '''
        try:
            params['smiles'] = smi
            for i in range(3):
                print(f'Trying to get response for {i + 1} time')
                resp = requests.get(self.host + '/api/treebuilder/', params=params, verify=False)
                print(resp)
                if 'error' not in resp.json().keys():
                    break
        except Exception as e:
            print(f'Error getting tree for {smi}: {e}')
            return -1
        return resp.json()

if __name__ == '__main__':
    HOST = 'http://XX.XX.XX.XX' # Replace with your address with ASKCOS

    # Create the tree builder
    tb = TreeBuilder(HOST)
    smi = 'CC1=C(COC2=CC(OC)=C(CN3CCCC[C@H]3C(O)=O)C(OC)=C2)C=CC=C1C4=CC=CC=C4'

    # Parameters for the query
    params = {
        'max_depth': 9,
        'max_branching': 25,
        'expansion_time': 60,
        'max_ppg': 100,
        'template_count': 1000,
        'max_cum_prob': 0.9999,
        'chemical_property_logic': 'none',
        'max_chemprop_c': 0,
        'max_chemprop_n': 0,
        'max_chemprop_o': 0,
        'max_chemprop_h': 0,
        'chemical_popularity_logic': 'none',
        'min_chempop_reactants': 5,
        'min_chempop_products': 5,
        'filter_threshold': 0.0001,
        'return_first': 'true' # default is false
    }

    # Build the tree for the given smiles and parameters
    js = tb.build_tree(smi, params)
    # Print the tree as a json
    pprint(js)