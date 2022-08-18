# Properties

## Background

This module is designed to automatically mine and calculate a series of chemical properties upon addition of a novel compound to the CIPHER database. This is accomplished using 2 methods. The first is making a series of requests to the PubChem REST API to mine hosted chemical property information from the PubChem Database if the compound is present in PubChem. The second is using the RDKit cheminformatics package to calculate a series of chemical properties and descriptors. Some of the properties overlap, however this is done intentionally as a way to ensure certain chemical properties and descriptors will be always be present regardless of whether or not the compound is present in PubChem. From here, these properties are organized into two embedded documents (one for PubChem and one for RDKit) and easily accessible via scripting or CIPHER REST API access for further use (i.e. chemical descriptors for machine learning features).

## Document Structure

`properties`

| Variable Name | Type | Description |
| ------------- | ---- | ---------- |
| `_id`| `str` | The InChIKey of the compound |
| `smiles` | `str` | The SMILES string of the compound |
| `pubchem` | `Object`| The PubChem entry for the compound |
| `rdkit`| `Object` | The RDKit entry for the compound |

`Pubchem (embedded object)`

| Variable Name | Type | Description |
| ------------- | ---- | ---------- |
| `MolecularFormula`| `str` | The molecular formula of the compound |
| `MolecularWeight` | `str` | The molecular weight of the compound |
| `XLogP` | `str` | Octanol-water partition coefficient of the compound |
| `MonoisotopicMass` | `str` | The monoisotopic mass of the compound |
| `TPSA` | `float` | The topological polar surface area of the compound (square angstroms)|
| `Complexity` | `int` | A synthetic complexity score calculated for the compound |
| `Charge` | `int` | The charge of the compound |
| `HBondDonorCount` | `int` | The number of hydrogen bond donor sites the compound has |
| `HBondAcceptorCount` | `int` | The number of hydrogen bond acceptor sites the compound has |
| `RotatableBondCount` | `int` | The number of rotatable bonds the compound has |
| `HeavyAtomCount` | `int` | The number of heavy atoms present in the compound |
| `IsotopeAtomCount` | `int` | The number of atoms which are isotopes in the compound |
| `AtomStereoCount` | `int` | The number of atoms which serve as stereoisomeric centers |
| `DefinedAtomStereoCount` | `int` | The number of defined atoms which serve as stereoisomeric centers |
| `UndefinedAtomStereoCount` | `int` | The number of undefined atoms which serve as stereoisomeric centers |
| `BondStereoCount` | `int` | The number of bonds involved in stereoisomeric centers |
| `DefinedAtomStereoCount` | `int` | The number of defined bonds which serve as stereoisomeric centers |
| `UndefinedAtomStereoCount` | `int` | The number of undefined bonds which serve as stereoisomeric centers |
| `CovalentUnitCount` | `int` | The number of covalent units as defined by PubChem |
| `Volume3D` | `float` | The 3D volume of the compound (square angstroms) |
| `XStericQuadrupole3D` | `float` | |
| `YStericQuadrupole3D` | `float` | |
| `ZStericQuadrupole3D` | `float` | |
| `FeatureCount3D` | `int` | | 
| `FeatureAcceptorCount3D` | `int` | |
| `FeatureDonorCount3D` | `int` | |
| `FeatureAnionCount3D` | `int` | |
| `FeatureCationCount3D`| `int` | |
| `FeatureRingCount3D` | `int`| |
| `FeatureHydrophobeCount3D`| `int`| |
| `ConformerModelRMSD3D` | `float`| |
| `EffectiveRotorCount3D`| `int`| |
| `ConformerCount3D` | `int` | |
| `Fingerprint2D` | `str`| |
| `modified` | `DateTime` | | 

`RDKit (embedded object)`

| Variable Name | Type | Description |
| ------------- | ---- | ---------- |
| `cipher_mid` | `str` | CIPHER model ID for the specific version of RDKit used for property calculation |
| `MolWt` | `float`| The approximate molecular weight of the compound|
| `ExactMolWt`| `float` | The exact molecular weight of the compound (precise calculation) |
| `HeavyAtomMolWt` | `float` | The average molecular weight of the molecule ignoring hydrogens |
| `MaxPartialCharge` | `float` | The maximum calculated value for the partial charge of the molecule|
| `MinPartialCharge` | `float` | The minimum calculated value for the partial charge of the molecule |
| `NumRadicalElectrons` | `int` | The number of radical electrons in the compound |
| `NumValenceElectrons` | `int` | The number of valence electrons in the compound |
| `MolLogP` | `float` | The calculated octanol-water coefficient of the compound |
| `MaxQED` | `float` | |
| `MeanQED` | `float` | |
| `NoneQED` | `float` | |
| `NHOHCount` | `int` | |
| `NOCount` | `int` | |
| `NumHAcceptors` | `int` | |
| `NumHDonors` | `int` | | 
| `NumHeteroatoms` | `int` | |
| `NumRotatableBonds` | `int`| |
| `RingCount` | `int` | |
| `modified` | `DateTime`| | 


## Installation and Setup

There are no specific setup instructions for the identifiers module. Make sure you have a Conda environment active with packages in the `cipher_db_env.yml` added to the environment. RDKit is present in that Conda environment, however if you need additional help installing RDKit, please refer to the [following guide](http://www.rdkit.org/docs/Install.html).

**Note:** In order for the PubChem REST API mining portion of the module to function, a working internet connection is required. While we never experienced any issues during testing, a timeout between PubChem REST API requests may be needed to prevent being locked out (use the [python `sleep()` function](https://docs.python.org/3/library/time.html#time.sleep)). See the [PubChem REST API volume](https://pubchemdocs.ncbi.nlm.nih.gov/programmatic-access#_RequestVolumeLimitations) rules for more information.

## Running the Properties Trigger

To run the trigger, simply run the `trigger.py` python script in the `cipher_properties` folder with the Conda environment active

## Properties Trigger Functionality

The trigger functions in the following manner
1. A compound is added to the compounds collection of the database via scripting, the website, or the REST API
2. The properties trigger watches for this addition
3. Two actions occur
    1. A script attempts to mine relevant information for the compound from PubChem (if the compound exists in PubChem)
    2. A script calculates chemical descriptors using RDKit
4. The information is added to the properties collection of the database.

## Modification

For modifications to PubChem properties follow the steps below:

1. Add the corresponding field to the `Pubchem` document class in `docs.py`
2. Add the property to the `full_property_list` variable in the `__format_request_url()` method in `properties.py` (make sure it is a listed property in the [PubChem REST API guide](https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest))

For modifications to RDKit properties follow the steps below:

1. Add the corresponding field to the `RDKit` document class in `docs.py`
2. Import the new RDKit descriptor calculation method at the top of `properties.py`
3. In the `insert_properties_from_smiles()` method in `properties.py` assign the field you created to the method output for the `RDKit()` object `rd`