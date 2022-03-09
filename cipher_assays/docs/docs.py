import mongoengine as me
import datetime
import shortuuid

class Assays(me.DynamicDocument):
    id = me.StringField(required=True,primary_key=True)
    inchikey = me.StringField(required=True)
    source = me.StringField(required=True)
    url = me.StringField()
    modified = me.DateTimeField(default=datetime.datetime.utcnow)


def gen_unique_assay_id():
    su = shortuuid.ShortUUID(alphabet="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    assigned = False
    while not assigned:
        rand_id = su.random(length=6)
        if Assays.objects.with_id(rand_id).count() == 0:
            return rand_id
