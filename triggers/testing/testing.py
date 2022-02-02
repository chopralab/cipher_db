import sys
sys.path.insert(0, "../schema/")
import urllib
import pymongo
from properties import pubchem, rdkit

login = open('../../utils/login.txt','r')
username = login.readline().replace('\n','')
username = urllib.parse.quote(username)
password = login.readline().replace('\n','')
password = urllib.parse.quote(password)

mongo_login = 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
client = pymongo.MongoClient(mongo_login)

testing = client['cipher_aspire']['testing']
properties = client['cipher_aspire']['properties']
compounds = client['cipher_aspire']['compounds']

rdkit_mid = "002"

for document in compounds.find():
    rdkit.insert(document["inchikey"], properties, compounds, rdkit_mid)

# rdkit.insert("ATLYLVPZNWDJBW-NHYNNZIHSA-N", testing, compounds, "002")