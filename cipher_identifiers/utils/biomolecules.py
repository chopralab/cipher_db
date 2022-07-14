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
    gen_unique_biomol_id
)

# Import ME document
from cipher_identifiers.docs.docs import Biomolecules

if __name__ == "__main__":
    import argparse
    import mongoengine as me
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--testing", action="store_true", help="Whether or not to run on the test database")
    parser.add_argument("--mid", type=str)
    parser.add_argument("--name", type=str)
    parser.add_argument("--uniprot-id", type=str)
    args = parser.parse_args()

    if args.testing:
        URI = os.environ["TESTING_URI"]
    else:
        URI = os.environ["MONGO_URI"]

    me.connect(host=URI)

    bmid = gen_unique_biomol_id()
    biomol = Biomolecules(cipher_bmid=bmid, cipher_mid=args.mid, name=args.name)
    biomol.save()
    print("Biomolecule saved with cipher_bmid: {}, name: {}, and cipher_mid: {}".format(bmid, args.name, args.mid))
