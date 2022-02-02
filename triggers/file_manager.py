import sys
sys.path.insert(0, 'schema/')
from binding import pubchem, desi_exp
import pymongo
import urllib
from os.path import exists

# Script for large file managment applications (insertion, deletion, editing, etc.)

login = open('../utils/login.txt','r')
username = login.readline().replace('\n','')
username = urllib.parse.quote(username)
password = login.readline().replace('\n','')
password = urllib.parse.quote(password)

mongo_login = 'mongodb+srv://' + username + ':' + password + '@aspirecluster0.hmj3q.mongodb.net/cipher_aspire?retryWrites=true&w=majority'
client = pymongo.MongoClient(mongo_login)

testing = client['cipher_aspire']['testing']
binding = client['cipher_aspire']['binding']

def insert_pubchem_data():
    mu_or_chembl_path = "../data/OPIOID_RECEPTORS/MU/mu_opioid_receptor_chembl_drugs.json"
    mu_or_drugbank_path = "../data/OPIOID_RECEPTORS/MU/mu_opioid_receptor_drugbank_drugs.json"

    kappa_or_chembl_path = "../data/OPIOID_RECEPTORS/KAPPA/kappa_opioid_receptor_chembl_drugs.json"
    kappa_or_drugbank_path = "../data/OPIOID_RECEPTORS/KAPPA/kappa_opioid_receptor_drugbank_drugs.json"
    kappa_or_IUPHAR_ligands_path = "../data/OPIOID_RECEPTORS/KAPPA/kappa_opioid_receptor_IUPHAR_ligands.json"

    nociceptin_or_chembl_path = "../data/OPIOID_RECEPTORS/NOCICEPTIN/nociceptin_opioid_receptor_chembl_drugs.json"
    nociceptin_or_drugbank_path = "../data/OPIOID_RECEPTORS/NOCICEPTIN/nociceptin_opioid_receptor_drugbank_drugs.json"
    nociceptin_or_IUPHAR_ligands_path = "../data/OPIOID_RECEPTORS/NOCICEPTIN/nociceptin_opioid_receptor_IUPHAR_ligands.json"

    sigma_nor_chembl_path = "../data/OPIOID_RECEPTORS/SIGMA/sigma_non_opioid_receptor_chembl_drugs.json"
    sigma_nor_drugbank_path = "../data/OPIOID_RECEPTORS/SIGMA/sigma_non_opioid_receptor_drugbank_drugs.json"
    sigma_nor_IUPHAR_ligands_path = "../data/OPIOID_RECEPTORS/SIGMA/sigma_non_opioid_receptor_IUPHAR_ligands.json"

    print(exists(mu_or_chembl_path))
    print(exists(mu_or_drugbank_path))

    print(exists(kappa_or_chembl_path))
    print(exists(kappa_or_drugbank_path))
    print(exists(kappa_or_IUPHAR_ligands_path))

    print(exists(nociceptin_or_chembl_path))
    print(exists(nociceptin_or_drugbank_path))
    print(exists(nociceptin_or_IUPHAR_ligands_path))

    print(exists(sigma_nor_chembl_path))
    print(exists(sigma_nor_drugbank_path))
    print(exists(sigma_nor_IUPHAR_ligands_path))

    pubchem.insert_from_json_file(mu_or_chembl_path, 'chembl', 'mu', binding)
    pubchem.insert_from_json_file(mu_or_drugbank_path, 'drugbank', 'mu', binding)

    pubchem.insert_from_json_file(kappa_or_chembl_path, 'chembl', 'kappa', binding)
    pubchem.insert_from_json_file(kappa_or_drugbank_path, 'drugbank', 'kappa', binding)
    pubchem.insert_from_json_file(kappa_or_IUPHAR_ligands_path, 'iuphar', 'kappa', binding)

    pubchem.insert_from_json_file(nociceptin_or_chembl_path, 'chembl', 'nociceptin', binding)
    pubchem.insert_from_json_file(nociceptin_or_drugbank_path, 'drugbank', 'nociceptin', binding)
    pubchem.insert_from_json_file(nociceptin_or_IUPHAR_ligands_path, 'iuphar', 'nociceptin', binding)

    pubchem.insert_from_json_file(sigma_nor_chembl_path, 'chembl', 'sigma', binding)
    pubchem.insert_from_json_file(sigma_nor_drugbank_path, 'drugbank', 'sigma', binding)
    pubchem.insert_from_json_file(sigma_nor_IUPHAR_ligands_path, 'iuphar', 'sigma', binding)

def insert_desi_data():
    desi_data_path = "../data/DESI/DESI_RESULTS_27Jan2022.tsv"
    test_path = "../data/DESI/test.tsv"
    print(exists(desi_data_path))
    if exists(desi_data_path):
        desi_exp.insert_from_tsv_file(test_path, testing)

def main():
    insert_desi_data()

if __name__ == "__main__":
    main()