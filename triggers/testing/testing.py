import sys
sys.path.insert(0, "../schema/")
import urllib
import pymongo
from properties import pubchem, rdkit
from identifiers import Compounds

login = open('../../utils/login.txt','r')
username = login.readline().replace('\n','')
username = urllib.parse.quote(username)
password = login.readline().replace('\n','')
password = urllib.parse.quote(password)

mongo_login = 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_testing?retryWrites=true&w=majority'
client = pymongo.MongoClient(mongo_login)

#testing = client['cipher_aspire']['testing']
#properties = client['cipher_aspire']['properties']
#compounds = client['cipher_aspire']['compounds']

Compounds.insert("CCC(=O)N(C1CCN(CC1)CCC2=CC=CC=C2)C3=CC=CC=C3", name="Fetanyl")
