# Roadmap_RegulatoryVariation

**note:** Additional readmes and examples are provided in each of these directories. 

#Feature Generation:
Code, examples and information on feature generation is provided in the **genFeaturesForDiscrete** data directory here.

#Prediction Generation:
We used the code provided in https://github.com/Roy-lab/HiC-Reg to generate predictions.
L-HiC-Reg simply uses features from 1 Mb regions as input into Hi-C reg. 
We used the code in **generateLocalPredictions** (see README there) to extract the features in a 1 Mb region from the whole of the output from the feature generation code.

#Significant interaction calling: The code for calling significant interactions with the binomial method is provided in the **sigcallinter** directory.

#Scoring nodes based on significant interactions with SNPs: The code for mapping significant interactions to SNPs and scoring the genes based on their significant interactions with SNPs is provided in the **mapregionpair_snp** directory. This is what is used as input into the graph diffusion and MTGC pipeline.

#Network Generation: The code for generating the cell-type specific networks is included in the **NetworkGeneration** directory. Note that this includes the distal nd proximal networks. 

#Graph diffusion and eigenvector calculation

#Muscari
After the eigegnvectors are calculated, to do the multi-task clustering we used MUSCARI
https://github.com/Roy-lab/Muscari
