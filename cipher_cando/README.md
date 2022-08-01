# CANDO
Computational Analysis of Novel Drug Opportunities

## Background

CANDO is a unique computational drug discovery, design, and repurposing platform
that analyses drug interactions on a proteomic scale, adhering to the multitarget drug theory,
for the purposes of shotgun drug discovery and repurposing, i.e., to evaluate every drug 
for every disease. [1-13]

The platform relates small molecules based on their computed interactions with all protein 
structures, known as an interaction signature with the hypothesis being
drugs with similar interaction signatures will have similar behavior 
in biological systems and will therefore be useful against the same indications.

For this project, the interactomic interaction signature has been distilled 
to a subset of analgesic and addiction sensitive biomarkers using higher-order analytics 
and multiscale networks paired with domain expertise. This signature can be modified 
as additional experimental results are obtained herein to confirm correlation between 
analgesia/non-addictive properties and other human protein biomarkers. Our holistic consideration of 
all interatomic proteins, as opposed to a single protein, enables assessment of pain and 
addiction on a broader scale which will increase our understanding of the necessary 
interactomic profile, i.e., mechanism of action and off-target effects of non-addictive 
analgesic molecular entities, allowing for optimization to a potent and safe pain therapeutic.

## Install

You may download the source code via the releases or cloning the git repository.

**OR**

We suggest using anaconda to install the CANDO package, as this is the easiest 
and quickest way to start using our platform! 

The CANDO package relies on multiple "conda-forge" dependencies. Therefore, 
we require that you add "conda-forge" to your anaconda channels: 

`conda config --add channels conda-forge`

Create a new environment in which to install CANDO

`conda env create -n cando`

Then you can install CANDO using the following command:

`conda install -c ram-compbio cando`

## CANDO integration

Some functionality of CANDO is integrated directly into the CIPHER Database.
These features include: human proteome library, BANDOCK (BioANalytical DOCKing),
multiscale analysis and prediction of therapeutic compounds (canpredict).

### Human proteome library

We have curated a library of 8,385 protein structures from the *Homo sapien* proteome. 
This library consists of solved (5,317) and modeled (3,068) protein structures. The modeled 
structures were predicted using I-TASSER, a homology-based protein structure prediction software. [14-16]
All protein structures in this library have collated information regarding binding sites that 
were calculated using COACH, a homology-based binding site prediction software from the I-TASSER suite. [16-18]

We have also modeled an additional 9 human proteins, not previously in our library, using AlphaFold 2. [19-20] 
These proteins are associated with analgesia and/or addiction/opioid use disorder. 
The list of proteins includes:
- Dopamine D2 receptor (long)
- Dopamine D2 receptor (short)
- Dopamine D3 receptor
- Mu-type opioid receptor
- Delta-type opioid receptor
- Kappa-type opioid receptor
- Nociceptin receptor
- AMPA-selective glutamate receptor 2
- Glutamate [NMDA] receptor subunit epsilon-2

*This library of protein structures is all pre-compiled, inherent to the CANDO software, and fully 
integrated into the CIPHER DB. All compounds entered into the database, via web portal or API, will have 
interaction scores calculated for all protein structures in this library using BANDOCK.*

### BANDOCK

BANDOCK (BioANalytical DOCKing) is a bioinformatic docking protocol that compares the structures of query 
drugs to all ligands known to bind to a given site on a protein. Specifically, the COACH algorithm is used 
to elucidate potential binding sites on each query protein, which uses a consensus approach via three 
different complementary algorithms that consider substructure or sequence similarity to known binding sites 
in the PDB. COACH output includes a set of cocrystallized ligands for each potential binding site, which 
are then compared to a compound/drug of interest using chemical fingerprinting methods that binarize the 
presence or absence of particular molecular substructures. The maximum Tanimoto coefficient between the 
binary vectors of the query compound and the set of all predicted binding site ligands for a protein serve 
as a proxy for the binding strength.

CANDO is agnostic to the scoring protocol used and any other may be relatively easily plugged 
in to generate interaction matrices and processed in a similar manner. The key advantages of 
BANDOCK is that it offers fast and reliable interation scores
which enables accurate predictions across billions of interaction calculations in minutes

*This protocol is fully integrated into the CIPHER DB and will be automatically trigger upon entering a new 
compound into the database. When triggered, the interaction scores for the query compound against all human 
proteins in our library will be calculated and individually populated in the database.*

### CANPREDICT

An important part of CANDO is generating putative drug candidates for a specific disease and predicting 
indications for which added or current drugs in our library can be therapeutic. We can do this with the 
canpredict functions. [1]

Generating putative drug candidates for a specific disease is one way we may want to use the predictive 
power of our platform. For this, we use the `canpredict_compounds()` function, which uses a consensus 
method to rank putative candidates based on how many times they show up as similar to drugs that are associated 
to the disease, within some defined cutoff. The default cutoff is 10 (the most stringent from benchmarking), 
but this can be varied. In its current implementation, `canpredict_compounds()` ranks compounds based on the 
consensus count and uses the average rank of those as the tie-breaker.

Alternatively, we can define a tailored interactomic signature containing specific proteins and objective 
interaction scores which represents the optimal interactomic signature for a non-addictive analgesic. 
This signature can be used as input for the k nearest neighbors (knn) algorithm, which then identifies the 
compounds most similar to the desired interactomic signature.

*The custom interactomic signature knn protocol is integrated in CIPHER DB and is automatically trigger to run 
once a day to update an active list of the best candidate non-addictive analgesic therapeutics based on the 
empirically optimized interactomic signature.*

## References
1. Mangione W, Falls Z, Chopra G, Samudrala R. (2020) cando.py: Open Source Software for Predictive Bioanalytics of Large Scale Drug-Protein-Disease Data. Journal of chemical information and modeling (Sep), 60(9): 4131-4136. doi:10.1021/acs.jcim.0c00110
2. Mangione W, Falls Z, Melendy T, Chopra G, Samudrala R. (2020) Shotgun drug repurposing biotechnology to tackle epidemics and pandemics. Drug discovery today (Jul), 25(7): 1126-1128. doi:10.1016/j.drudis.2020.05.002
3. Schuler J, Falls Z, Mangione W, Hudson ML, Bruggemann L, Samudrala R. (2022) Evaluating the performance of drug-repurposing technologies. Drug discovery today (Jan), 27(1): 49-64.
4. Overhoff B, Falls Z, Mangione W, Samudrala R. (2021) A Deep-Learning Proteomic-Scale Approach for Drug Design. Pharmaceuticals (Basel, Switzerland) (Dec), 14(12). doi:10.3390/ph14121277
5. Fine J, Lackner R, Samudrala R, Chopra G. Computational chemoproteomics to understand the role of selected psychoactives in treating mental health indications. Scientific Reports 9, 1315, 2019.
6. Schuler J, Samudrala R. Fingerprinting CANDO: Increased accuracy with structure and ligand based shotgun drug repurposing. ACS Omega 4: 17393-17403, 2019.
7. Falls Z, Mangione W, Schuler J, Samudrala R. Exploration of interaction scoring criteria in the CANDO platform. BMC Research Notes 12: 318, 2019.
8. Mangione W, Samudrala R. Identifying protein features responsible for improved drug repurposing accuracies using the CANDO platform: Implications for drug design. Molecules 24: 167, 2019.
9. Chopra G, Samudrala R. Exploring polypharmacology in drug discovery and repurposing using the CANDO platform. Current Pharmaceutical Design 22: 3109-3123 2016.
10. Sethi G, Chopra G, Samudrala R. Multiscale modelling of relationships between protein classes and drug behavior across all diseases using the CANDO platform. Mini Reviews in Medicinal Chemistry 15: 705-717, 2015.
11. Minie M, Chopra G, Sethi G, Horst J, White G, Roy A, Hatti K, Samudrala R. CANDO and the infinite drug discovery frontier. Drug Discovery Today 19: 1353-1363, 2014.
12. Horst JA, Laurenzi A, Bernard B, Samudrala R. Computational multitarget drug discovery. Polypharmacology 263-301, 2012.
13. Jenwitheesuk E, Horst JA, Rivas K, Van Voorhis WC, Samudrala R. Novel paradigms for drug discovery: Computational multitarget screening. Trends in Pharmacological Sciences 29: 62-71, 2008.
14. W Zheng, C Zhang, Y Li, R Pearce, EW Bell, Y Zhang. Folding non-homology proteins by coupling deep-learning contact maps with I-TASSER assembly simulations. Cell Reports Methods, 1: 100014 (2021).
15. J Yang, R Yan, A Roy, D Xu, J Poisson, Y Zhang. The I-TASSER Suite: Protein structure and function prediction. Nature Methods, 12: 7-8 (2015).
16. J Yang, Y Zhang. I-TASSER server: new development for protein structure and function predictions. Nucleic Acids Research, 43: W174-W181 (2015).
17. Jianyi Yang, Ambrish Roy, and Yang Zhang. Protein-ligand binding site recognition using complementary binding-specific substructure comparison and sequence profile alignment, Bioinformatics, 29:2588-2595 (2013).
18. Jianyi Yang, Ambrish Roy, and Yang Zhang. BioLiP: a semi-manually curated database for biologically relevant ligand-protein interactions, Nucleic Acids Research, 41: D1096-D1103 (2013).
19. Jumper, J et al. Highly accurate protein structure prediction with AlphaFold. Nature (2021).
20. Varadi, M et al. AlphaFold Protein Structure Database: massively expanding the structural coverage of protein-sequence space with high-accuracy models. Nucleic Acids Research (2021).
