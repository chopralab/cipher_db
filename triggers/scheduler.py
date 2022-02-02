import pymongo
import urllib
import sys
from multiprocessing import Process

# Script for database trigger scheduling and resource management

login = open('../utils/login.txt','r')
username = login.readline().replace('\n','')
username = urllib.parse.quote(username)
password = login.readline().replace('\n','')
password = urllib.parse.quote(password)

mongo_login = 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
client = pymongo.MongoClient(mongo_login)

# Identifiers
compounds = client['cipher_aspire']['compounds']
substances = client['cipher_aspire']['substances']
experiments = client['cipher_aspire']['experiments']
models = client['cipher_aspire']['models']
external_databases = client['cipher_aspire']['external_databases']
biomolecules = client['cipher_aspire']['biomolecules']
binding_sites = client['cipher_aspire']['binding_sites']

# Primary Data
properties = client['cipher_aspire']['properties']
ord = client['cipher_aspire']['ord']
reactivity = client['cipher_aspire']['reactivity']
binding = client['cipher_aspire']['binding']

# Machine Learning Features
features = client['cipher_aspire']['features']

# Website Login
users = client['cipher_aspire']['users']