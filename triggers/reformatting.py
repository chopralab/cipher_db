from matplotlib import collections
import pymongo
import urllib
import sys
from multiprocessing import Process

# Script for reformatting out-of-date entries

login = open('../utils/login.txt','r')
username = login.readline().replace('\n','')
username = urllib.parse.quote(username)
password = login.readline().replace('\n','')
password = urllib.parse.quote(password)

mongo_login = 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
client = pymongo.MongoClient(mongo_login)

# Testing
testing = client['cipher_aspire']['testing']

# Identifiers
compounds = client['cipher_aspire']['compounds']


def correct_date_field(collection):
    collection.update_many(
        {},
        {"$unset": {"date_modified": 1}}
    )
    collection.update_many(
        {},
        {"$currentDate": {"modified": True}}
    )

def rename(collection, old, new):
    collection.update_many(
        {},
        {"$rename": {old: new}}
    )

if __name__ == "__main__":
    rename(compounds, "IUPACname", "iupac")