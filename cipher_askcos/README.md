# ASKCOS

## Background

ASKCOS is a software package for the prediction of feasible synthetic routes towards a desired compound and associated tasks related to synthesis planning. It is used in the CIPHER database to generate feasible retrosynthetic pathways for novel non-addictive opioids added to the database. 

## ASKCOS Installation and Setup

The setup for the reactivity is done in two (independent) parts:

1. Deploy ASKCOS on some accessible IP address (private or public depending your use case) following these directions: https://github.com/ASKCOS/askcos-deploy
    - Make note of this IP address and save it somewhere
    - (optional) set up a user account on the ASKCOS server and save these credentials
2. Set up the cipher reactivity service: https://github.com/chopralab/cipher_db/tree/master/cipher_reactivity
    - install the package to your environment: `$ pip install -e .`
    - set the following environment variables in your environment:
        1. `ASKCOS_HOST` : the IP address of your ASKCOS server. This is the web service backend to which all retrosynthetic analysis jobs will be submitted
        2. (optional) `ASCKOS_USERNAME` / `ASKCOS_PASSWORD` : the username and password of the account to which tree builder jobs will be saved in the case you would like to manually inspect these results later on 

## ASKCOS Configuration

Two default configuration files are provided. They are stored inside the package and are automatically imported along with the package, so nothing needs to be performed by a client. They are located under `cipher_reactivity/data` subdirectory:

1. fpscores.json : this is a JSON file containing a database of SA scores and corresponding fingerprint hashes. I.e., it is a List<Tuple<float, int, ...>> where the 0th index of each tuple is a float (the SA score) the 1st through nth indices are all morgan fingerprint hashes that correspond to that SA score. In principle, you could supply your own database, but the SAScore was demonstrated with the database currently in the package, so switching it out would make the resulting SA scores incomparable to standard SA scores (i.e., it becomes an "apples to oranges" comparison if you change this file.)

2. tree_params.toml : a TOML config file containing default retrosynthesis tree builder job parameters. Altering this file in the code base can change the expansion time or reaction rules that are used when performing the analysis

## Running the ASKCOS Trigger

Once you have a working instance of ASKCOS deployed and the configuration is complete, running the `$trigger.py` python script should start the trigger and it will automatically watch for insertions into the main "compounds" collection of the database, automatically updating the "difficulty" and "retrosynthesis" collections with the corresponding documents.

## ASKCOS Trigger Functionality

The ASKCOS Trigger functions in the following manner upon insertion of a compound to the database:

1. A new compound is inserted into the central `compounds` collection of the database
2. This trigger is watching for that insertion and picks up the inchikey and SMILES string of that compound
3. A Difficulty document is built and saved for the corresponding compound
4. A Retrosynthesis document is built and saved for the corresponding compound
    1. An ASKCOS tree builder job is submitted to the backend ASKCOS service for the given molecule
    2. For each synthetic tree returned, build a SyntheticTree document. A single SyntheticTree corresponds to a single retrosynthetic analysis of the compound. Nodes in this tree are either chemicals or reactions. A chemical node always has 1 edge leading away from it to a reaction node, and a reaction node has `n` edges leading away from to `n` chemical nodes, one for each of the precursors in this reaction.
    3. After a tree is traversed and all corresponding documents are created (Root, ChemicalNode(s), and ReactionNode(s)), the tree is visualized in an image. This image is rendered and saved into the `image` field of the SyntheticTree document 

# Mongo DB Document Structure Used by ASKCOS (Retroysnthesis)

The `retrosynthesis` collection is comprised of `Retrosynthesis` documents. A single document contains the following fields:

| Variable Name | Type | Description |
| -------------- | ---- | ----------- |
| `inchikey`| `str`  | The inchikey of the compound (the primary key in this collection) |
| `smiles`  | `str`  | The SMILES string of the compound |
| `trees`   | `List<SyntheticTree>` | A list of SyntheticTree embedded documents, one for each individual retrosynthetic tree enumerated for this compound |
| `task_id` | `int`  | The celery task ID logged in the ASKCOS backend service|

In turn, a `SyntheticTree` is a document with the following fields:

| Variable Name | Type | Description |
| ------------- | ---- | ----------- |
| `depth`| `int` | The depth of the synthetic tree |
| `precursor_cost` | `float` | The cost of the precursors in this tree |
| `score` | `float`| The feasibility score of this tree |
| `cluster_id` | `int` | The cluster ID of the tree for tree clustering analysis |
| `root` | `ChemicalNode` | The root chemical (product) of the tree |
| `image` | `Image` | An image of the synthetic tree |

Further, a `ChemicalNode` is an embedded document with structure:

| Variable Name | Type | Description |
| ------------- | ---- | ----------- |
| `smiles` | `str` | The SMILES string of the molecule |
| `rxn`    | `ReactionNode` | The reaction template used to construct this molecule |
| `ppg` | `float` | The price per gram of the molecule (if purchaseable) |
| `as_reactant` | `int` | The number of reaction templates in the template set that feature this molecule as a reactant |
| `as_product` | `int` | The number of reaction templates in the template set that feature this molecule as a product |
| `terminal` | `boolean` | Whether this node is a leaf of the tree |
| `chemical_id` | `int` |  The ID of this chemical in the ASKCOS buyables database, if it exists |

Going even deeper, a `ReactionNode` has the following structure:

| Variable Name | Type | Description |
| ------------- | ---- | ----------- |
| `smiles` | `str` | The reaction SMILES |
| `tforms` | `list[str]` | A secondary reaction ID value within the template set |
| `tsources` | `list[str]` | The source of this reaction template set, i.e., the name of the template set |
| `template_score` | `float` | The template of this reaction from the template prioritizer model (0-1) |
| `plausibility` | `float` | The plausibility score of the reaction (0-1) |
| `rank` | `int` | The rank of this template among all possible templates to choose from in the template set |
| `num_examples` | `int` | The number of examples of this reaction in the template set |
| `necessary_reagent` | `str` | Any additional reagents required to run the reaction |
| `precursor_smiles` | `str`  | The SMILES strings (period-separated) of the reactants in this reaction |
| `rms_molwt` | `float` | The root mean square molecular weight of the product molecules |
| `num_rings` |  `int`  | The number of rings in the product |
| `scscore`   | `float` | The SC score of the product template | 
| `rxn_id`    |  `str`  | The ID of the reaction in the template set |

## Helper Code

Please check the docstrings located in Python files for additional information:

## Modification Instructions

Add instructions here on how to modify the toml file ...