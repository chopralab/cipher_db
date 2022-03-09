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
from cipher_identifiers.docs.docs import Models

if __name__ == "__main__":
    import argparse
    import mongoengine as me
    import os

    parser = argparse.ArgumentParser()
    parser.add_argument("--testing", action="store_true", help="Whether or not to run on the test database")
    parser.add_argument("--mid", type=str, required=True, help="The Cipher MID of the inserted model")
    parser.add_argument("--source", type=str, required=True, help="The source of the inserted model")
    parser.add_argument("--parameters", type=str, default="{}", help="The parameters of the inserted model")
    args = parser.parse_args()

    if args.testing:
        URI = os.environ["TESTING_URI"]
    else:
        URI = os.environ["MONGO_URI"]

    me.connect(host=URI)

    model = Models(cipher_mid=args.mid, source=args.source, parameters=args.parameters)
    model.save()
    print("Model saved with mid: {}, source: {}, and parameters: {}".format(args.mid, args.source, args.parameters))