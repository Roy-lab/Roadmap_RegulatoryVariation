# genFeaturesForDiscreteData

This version of genFeatures uses pre-discretized data of chromatin marks for enhancers and promoters. It takes in continuous data to compute the mean (before looping over the folds to write out the features per E-P pair) in the window of ALL of the E-P pairs - a step called pre-computing the mean. For each feature, it clusters these means for all pairs with k-means. It reorganizes the clusters by mean so the ID corresponds with the actual values in the cluster. It then restructures the output for k-means to make it faster to write out. Only after all of that will it write out the enhancer, promoter and window features for a given pair. The reason for all of this is because we found that computing the mean and then discretizing was a better window feature than the median of discrete data. 

Additionally, for motifs it takes the count (PIQ PIV > 0.5) in E P pairs, and the count divided by the window size for the window. 

Usage is the same as previous implementations. In order to work, the motifs have to be listed in the feature file first followed by the discrete data and continuous data in the same order. Motifs are labeled with 'M', discrete data is labeled with 'C' and continuous data is labeled with 'W'.

Example feature file format, along with example inputs and outputs are provided in the 'examples' folder.
