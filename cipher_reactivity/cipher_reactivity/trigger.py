import os
from cipher_reactivity import builder

import mongoengine as me
import pymongo as pmg

MONGO_CLIENT = me.connect(host=os.environ["MONGO_URI"])


def synthesis_trigger():
    cpds = MONGO_CLIENT.cipher_aspire.compounds
    try:
        with cpds.watch([{"$match": {"operationType": "insert"}}]) as stream:
            for change in stream:
                d = change["fullDocument"]
                inchikey = d["inchikey"]
                smi = d["smiles"]

                builder.update_difficulty(inchikey, smi)
                builder.update_retrosynthesis(inchikey, smi)
    except pmg.errors.PyMongoError as e:
        print(e)


if __name__ == "__main__":
    synthesis_trigger()
