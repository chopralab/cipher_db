# Identifiers

## Background

The identifiers module is the core system for recording information on which compounds are present in the database. All compounds referenced in any collection of the database have their record present in one or more of the collections the managed by the identifiers module. Additionally, all models and external databases used in the project have records present in one of the collections managed by the identifiers module. The collections which are managed by the module are listed below with a short background describing their purpose:

`compounds`: Stores identifying information on chemical compounds. This collection is indexed by InChi Key. 

| Variable Name | Type | Description |
| ------------- | ---- | ---------- |
| `_id` (primary key) | `str` | The unique InChi Key identifier of the compound |
| `name` | `str` | The assigned name of the compound if it is mined or provided |
| `smiles` | `str` | The SMILES string of the compound |
| `inchi` | `str` | The InChI string of the compound |
| `cid` | `str` | The PubChem CID of the compound if it is present in the database |
| `iupac` | `str` | The IUPAC name of the compound |
| `synonyms`| `List<str>` | A list of synonyms for the compound |
| `image` | `File` | The image of the compound | 
| `modified` | `DateTime` | The date and time the compound was added or last modified |


`substances`: Stores information on a compounds source (vendor, etc.)

| Variable Name | Type | Description |
| ------------- | ---- | ---------- |
| `cipher_sid` (primary key) | `str` | The unique ID for a specific substance |
| `inchikey` | `str` | The InChi Key of the substance |
| `smiles` | `str` | The SMILES string of the substance |
| `source` | `str` | The vendor (or other) source of the substance |
| `name` | `str` | The assigned name of the substance |
| `modified` | `DateTime` | The date and time the substance was added or last modified | 

`biomolecules`: Stores identifying information on proteins, pepties, etc. used in Cando signatures

| Variable Name | Type | Description |
| ------------- | ---- | ---------- |
| `cipher_bmid` (primary key) | `str`| The unique ID for a specific biomolecule |
| `name` | `str` | The assigned name of the biomolecule |
| `pbd_id` | `str` | The protein databank ID of the biomolecule (if present) |
| `chain_id` | `str` | The chain ID of the biomolecule (if applicable) |
| `uniprot_id` | `str` | The Uniprot ID of the biomolecule (if present) |
| `cipher_mid` | `str` | The computation model used to determine biomolecule structure (i.e. AlphaFold) |
| `modified` | `DateTime` | The date and time the biomolecule was added or last modified |

`binding sites`: Stores information on the bindings sites of biomolecules

| Variable Name | Type | Description |
| ------------- | ---- | ----------- |
| `cipher_bsid` | `str` | The unique ID for a specific binding site |
| `cipher_bmid` | `str` | The CIPHER BMID of the biomolecule with the given binding site |
| `template_pdb_id` | `str` | ... |
| `template_chain_id` | `str` | ... |
| `residue` | `str` | The residue on which the binding site is located |
| `coordinates` | `str` | The 3D coordinates of the binding site |
| `modified` | `str` | The date and time the binding site was added or last modified | 

`models`: Stores information on computational modules used in CIPHER

| Variable Name | Type | Description |
| ------------- | ---- | ---------- |
| `cipher_mid` | `str` | The unique ID for a specific computational model |
| `source` | `str` | The source of the model | 
| `parameters` | `str` | Any parameters which the model is run with |
| `modified` | `DateTime` | The date and time the model was added or last modified |

`external databases`: Stores information on external databases referenced in CIPHER

| Variable Name | Type | Description |
| ------------- | ---- | ---------- |
| `cipher_exdbid`| `str` | The unique ID for a specific external database used as a reference |
| `source` | `str` | The name of the external database |
| `url` | `str` | The URL of the external database (if applicable) |
| `modified` | `DateTime` | The date and time the URL was last modified |

## Installation and Setup

There are no specific setup instructions for the identifiers module. Make sure you have a Conda environment active with 