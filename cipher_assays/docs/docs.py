import mongoengine as me
import datetime
import shortuuid

from cipher_identifiers.docs.docs import validate_smiles, check_inchikey_in_compounds

class Assays(me.DynamicDocument):
    cipher_aid = me.StringField(required=True,primary_key=True)
    inchikey = me.StringField(required=True)
    smiles = me.StringField(required=True, validation=validate_smiles)
    source = me.StringField(required=True)
    receptor = me.StringField(required=True)

def gen_unique_assay_id():
    su = shortuuid.ShortUUID(alphabet="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    assigned = False
    while not assigned:
        rand_id = su.random(length=6)
        if Assays.objects.with_id(rand_id) is None:
            return rand_id
