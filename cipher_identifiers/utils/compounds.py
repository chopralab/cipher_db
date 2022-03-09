import requests
import json
from rdkit import Chem

if __name__ == "__main__":
    import sys
    sys.path.append("../../")

# Import Custom Errors
from cipher_identifiers.utils.errors import (
    InvalidRequestError,
    CompoundNotFoundError,
    InvalidSMILESError,
)

# Import all validation functions
from cipher_identifiers.docs.docs import (
    validate_smiles,
    check_inchikey_in_compounds,
    check_mid_in_models,
    check_bmid_in_biomolecules,
    check_bsid_in_binding_sites,
)

# Import ME document
from cipher_identifiers.docs.docs import Compounds


def get_ids_from_inchikey(inchikey):
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
        + inchikey
        + "/property/IUPACName,Title/json"
    )
    response = requests.get(url)
    if response:
        data = json.loads(response.text)
        cid = data["PropertyTable"]["Properties"][0]["CID"]
        iupac = data["PropertyTable"]["Properties"][0]["IUPACName"]
        name = data["PropertyTable"]["Properties"][0]["Title"]
        return int(cid), iupac, name
    else:
        raise InvalidRequestError(
            "Invalid Request - Request on InChi Key "
            + inchikey
            + " failed with status code "
            + response.status_code
        )


def id_compound_from_smiles(smiles):
    comp = Compounds()
    try:
        m = Chem.MolFromSmiles(smiles)
    except:
        raise InvalidSMILESError(
            "Invalid SMILES - The provided SMILES "
            + smiles
            + " cannot be converted into a valid molecule"
        )

    inchikey = Chem.MolToInchiKey(m)

    if Compounds.objects.with_id(inchikey) is not None:
        return False

    comp.inchikey = inchikey
    comp.smiles = Chem.MolToSmiles(m, isomericSmiles=False)
    comp.inchi = Chem.MolToInchi(m)

    try:
        comp.cid, comp.iupac, comp.name = get_ids_from_inchikey(inchikey)
    except Exception as e:
        print(e)

    comp.save()
    print(
        "Inserted Object with InChI Key {} to Compounds Collection of the database".format(inchikey)
    )
    return True


if __name__ == "__main__":
    import argparse