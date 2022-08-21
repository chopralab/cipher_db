## RCVAE compound design module

Reduced conditional variational autoencoder for novel chemical generation on a proteomic scale. 

## Background

We used the CANDO platform to predict interaction signatures for each compound in a training set against a library of nonredundant protein structures representing the human proteome. The interaction signatures have their dimensionality reduced in an autoencoder, which models the underlying correspondence between protein structures as they behave in the proteome. The reduced signatures are then used as labels for each training compound, which a generative conditional variational autoencoder model learns to reconstruct given a target interaction signature. This pipeline allows us to novel design compounds that modulate interactions on a proteomic scale as desired to generate behaviorally novel therapeutics. More information about implementation and benchmarking may be found [here](https://www.mdpi.com/1424-8247/14/12/1277/htm). 

## Prerequisites and Setup

The RCVAE design module is written in Python using Tensorflow. In order to run the script, you must have a Conda environment active with the following packages installed:

```(yaml)
tensorflow==1.15.0
scikit-learn==1.0.2
rdkit==2020.09.1.0
keras==2.3.1
googledrivedownloader==0.4
pandas==1.4.1
```

Please also install all requiried dependencies for these packages as prompted by `pip` or `conda`.

Then, install the CVAE model architecture by cloning [this repository](https://github.com/jaechanglim/CVAE.git) into the same folder as `rcvae.py` as a dependency for the design pipeline.

Finally, the first line of `rcvae.py` should be changed to the path of your local Conda environment. Upon running the script, all trained models will be automatically gathered from Google Drive. 

## Files

- `rcvae.py`: Control file for CANDO interaction signature reduction and CVAE integration/design.
- `models/ae.h5`: Trained autoencoder model for CANDO interaction signature reduction.
- `models/model_30.ckpt-30*`: Trained CVAE model for chemical design with reduced CANDO interaction signatures.

## Usage

With the conda environment active, navigate to the `cipher_rcvae` folder and run the `./rcvae.py` python file with the following required flags:

| Argument | Flag | Description |
| ------ | --------- | ----------- |
| `TEMPLATE`| `-t --template` | Path to template objective CANDO signature. |
| `DESIGNS` | `-d --designs` | Number of compounds to generate for the given CANDO signature. |
| `OUTPUT` | `-o --output`| Output file path for generated SMILES strings. |

## Modifications

As we provide a fully trained model for conditional compound generation, generalizable to any desired compound-proteome behavior, no substantial modifications need occur with the design script or models themselves. Users may however manipulate CANDO signatures in a variety of ways to produce optimally designed compounds which we leave to future development.  