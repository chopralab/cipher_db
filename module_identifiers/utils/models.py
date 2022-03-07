if __name__ == "__main__":
    import sys
    sys.path.append("../../")

# Import Custom Errors
from module_identifiers.utils.errors import (
    InvalidRequestError,
    CompoundNotFoundError,
    InvalidSMILESError,
)

# Import all validation functions
from module_identifiers.docs.docs import (
    validate_smiles,
    check_inchikey_in_compounds,
    check_mid_in_models,
    check_bmid_in_biomolecules,
    check_bsid_in_binding_sites,
)

# Import ME document
from module_identifiers.docs.docs import Models

if __name__ == "__main__":
    import argparse
    import mongoengine as me

    parser = argparse.ArgumentParser()
    parser.add_argument("--testing", action="store_true", help="Whether or not to run on the test database")
    parser.add_argument("--mid", type=str, required=True, help="The Cipher MID of the inserted model")
    parser.add_argument("--source", type=str, required=True, help="The source of the inserted model")
    parser.add_argument("--parameters", type=str, default="{}", help="The parameters of the inserted model")
    args = parser.parse_args()

    LOGIN = open("../../utils/login.txt", "r")
    USERNAME = LOGIN.readline().replace("\n", "")
    PASSWORD = LOGIN.readline().replace("\n", "")
    LOGIN.close()

    if args.testing:
        URI = (
            "mongodb+srv://"
            + USERNAME
            + ":"
            + PASSWORD
            + "@aspirecluster0.hmj3q.mongodb.net/cipher_testing?retryWrites=true&w=majority"
        )
    else:
        URI = (
            "mongodb+srv://"
            + USERNAME
            + ":"
            + PASSWORD
            + "@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority"
        )

    me.connect(host=URI)

    model = Models(id=args.mid, source=args.source, parameters=args.parameters)
    model.save()
    print("Model saved with mid: {}, source: {}, and parameters: {}".format(args.mid, args.source, args.parameters))