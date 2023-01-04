# HiC-Reg: In silico prediction of high-resolution Hi-C interaction matrices
HiC-Reg is a novel approach to predict contact counts from one-dimensional regulatory signals such as epigenetic marks and regulatory protein binding. HiC-Reg provides a powerful framework to generate high-resolution profiles of contact counts that can be used to study individual locus level interactions as well as higher-order organizational units of the genome.



## Step 1: Generate pair features as input for HiC-Reg:
### 1.0 Aggregate region-level features:
Program in Scripts/aggregateSignalInRegion/

#### Usage:
```
aggregateSignal REGION hg19.fa.fai countfile outputfile
```
#### Examples:
```
aggregateSignal hg19_5kbp_chr17.txt hg19.fa.fai wgEncodeBroadHistoneGm12878CtcfStdRawDataRep1_chr17.counts wgEncodeBroadHistoneGm12878CtcfStdRawDataRep1_chr17.txt
```

#### Arguments:
- REGION: file of 5kb bin with format like hg19_5kbp_chr17.txt
- countfile: four columns count file generated from bam files
  
#### Input Files:
1. hg19_5000bp_chr17.txt
```
chr17	CONVERT	gene	0	4999	.	+	.	chr17_0_4999
chr17	CONVERT	gene	5000	9999	.	+	.	chr17_5000_9999
chr17	CONVERT	gene	10000	14999	.	+	.	chr17_10000_14999
```
2. count file: wgEncodeBroadHistoneGm12878CtcfStdRawDataRep1_chr17.counts 
```
chr17	342	389	1
chr17	1005	1056	1
chr17	1060	1111	1
chr17	1346	1397	1
```
Column1: chromosome Column2: start position Column3: end position Column4: reads (tab deliminated)


### 1.1 Generate PAIR-CONCAT or WINDOW features:
Program in Scripts/genPairFeatures/

#### Usage:
```
./genDatasetsRH HiCsparseMatrix maxdistance foldcv splittype[regionwise|pairwise] featureinputfile correlation[yes|no] outputdir prerandomize_pairs[yes|no] featype[Window|Pconcat]
```

#### Arguments:
- HiCsparseMatrix: sparse hic matrix
- maxdistance: max genomic distance for pair of regions (e.g. 1000000)
- foldcv: number of CV folds (e.g. 5)
- splittype: 
  - regionwise--split the total regions in the sparse hic matrix into N folds. 
  - pairwise--split the total pairs in the sparse hic matrix into N folds.
- featureinputfile: input file with path for each feature signal, see Gm12878_norm_featurefiles_test.txt for example.
- correlation: calculate the correlation of features in region1 and features in region2 or not. (e.g. no)
- outputdir: path for output directory.
- prerandomize_pairs: pre-randomize the pairs in the sparse hic matrix or not. (e.g. yes)
- featype: 
  - Pconcat: generate feature signal for region1 and region2. 
  - Window: generate feature signal for region1 and region2, and average feature signal for the window between these two regions. 

#### Input Files:
1. HiCsparseMatrix: Gm12878_chr17_5kb_SQRTVC_counts_pairs_100.tab
```
chr17_0_5000	chr17_5000_10000	485.3207854051
chr17_0_5000	chr17_10000_15000	212.4988640932
chr17_0_5000	chr17_15000_20000	130.1059543547
```
Column1: region1 Column2: region2 Column3: HiC count (tab deliminated)

2. featureinputfile: Gm12878_norm_featurefiles_test.txt
```
H3k4me1	Gm12878_RawData_5000bp_seqdepth_norm_H3k4me1.txt	C
H3k4me2	Gm12878_RawData_5000bp_seqdepth_norm_H3k4me2.txt	C
```
Column1: mark name Column2: feature signal file Column3: source data type, C for Chip-seq data, M for motif data (tab deliminated)


#### Example: 
```
./genDatasetsRH Gm12878_chr17_5kb_SQRTVC_counts_pairs_100.tab 1000000 5 regionwise Gm12878_norm_featurefiles_test.txt no out/ yes Window
```



## Step2: Train HiC-Reg models and make predictions
### 2.1 Training mode:
#### Train on training data and Predict on test set:
```
./regForest -t Examples/Data/Gm12878_chr17_WINDOW_train0.txt -o Examples/out/ -k1 -l10 -n20 -b Examples/Data/prior_window.txt -d Examples/Data/Gm12878_chr17_WINDOW_test0.txt
```
### 2.2 Prediction mode:
#### Predict on test sets with Models trained:
```
./regForest -t Examples/Data/Gm12878_chr17_WINDOW_train0.txt -o out/ -k1 -l10 -n20 -b Examples/Data/prior_window.txt -d Examples/Data/Gm12878_chr17_WINDOW_test0.txt -s Examples/out/regtree_node
```
#### Arguments:
- -t is the training data
- -o is the output directory
- -k is maxfactorsize (1 should be OK)
- -l is leaf size (e.g. 10)
- -n is the number of trees (e.g. 20)
- -b is the prior structure (potential regulators), see Examples/Data/prior_window.txt 
- -d is the test data
- -s is the prefix of saved models (e.g. Examples/out/regtree_node)

#### Input Files:
1. training and test data genreated by Step 1: (tab deliminated)
```
Pair	H3k4me1_E	H3k4me2_E	H3k4me1_P	H3k4me2_P	H3k4me1_W	H3k4me2_W	Distance	Count
chr17_0_5000-chr17_10000_15000	0.714904	0.435131	3.76434	1.3964	0.83982	0.418721	5000	5.36363
chr17_0_5000-chr17_15000_20000	0.714904	0.435131	1.73117	0.869143	2.30208	0.90756	10000	4.87601
chr17_0_5000-chr17_25000_30000	0.714904	0.435131	0.881374	0.708816	1.73773	0.78015	20000	4.76591
chr17_0_5000-chr17_30000_35000	0.714904	0.435131	1.49854	0.760584	1.56646	0.765883	25000	4.5534
```
header of training and test data should be the same.

2. prior_window.txt is a file to list the feature set used for training models, that is like this: (tab deliminated)
```
H3k4me1_E    Count
H3k4me1_P    Count
H3k4me1_W    Count
H3k4me2_E    Count
H3k4me2_P    Count
H3k4me2_W    Count
Distance    Count
```
Column1: feature set (header of column 2 to second last column) in training/test data, like H3k4me1_E, H3k4me1_P, H3k4me1_W, Distance.

Column2: Count (i.e. the header of the last column in training/test data)

#### Output Files:
1. testset_error.txt
(Predictions on the test data)
- Column1: Pair of two regions
- Column2: True HiC Count for this pair in the log scale 	
- Column3: Predicted HiC Count for this pair in the log scale 	
- Column4: squared prediction error
- Column5: genomic distance of the two regions

2. trainset_error.txt
(Predictions on the training data)
- Column1: Pair of two regions
- Column2: True HiC Count for this pair in the log scale 	
- Column3: Predicted HiC Count for this pair in the log scale 	
- Column4: squared prediction error
- Column5: genomic distance of the two regions

3. regtree_node_0.txt
(Saved regression tree models for tree 0, examples in Examples/out/)


## Application: Make predictions in a new cell line
### Train models using training data in Gm12878 and make predictions in K562:
```
./regForest -t Examples/Data/Gm12878_chr17_WINDOW_train0.txt -o Examples/out/ -k1 -l10 -n20 -b Examples/Data/prior_window.txt -d Examples/Data/K562_chr17_WINDOW_test0.txt
```
#### Arguments:
- -t is the training data in cell line 1 (Gm12878)
- -d is the test data in cell line 2 (K562)
- -o is the output directory
- -k is maxfactorsize (1 should be OK)
- -l is leaf size (e.g. 10)
- -n is the number of trees (e.g. 20)
- -b is the prior structure (potential regulators), see Data/prior_merge.txt.
- -s is the prefix of saved models
