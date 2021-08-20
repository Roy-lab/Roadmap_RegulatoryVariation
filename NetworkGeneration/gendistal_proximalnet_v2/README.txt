This program takes in an enhancer-promoter interaction file where the promoter is mapped to a gene, a set of accessible peaks (from DNase1) and a set of peak to motif mappings and generates an output file comprising the TFs to genes based on if the TF's motifs are in the distal or proximal region associated with a gene. The distal and proximal are also parameter configs.

Usage: genProxDistNet interactions dnase1peaks motifinstances [distal|proximal] outputnetwork
Here interactions is the EP interactions, dnase1peaks are the set of Dnase1 peals, motif instances is peak name to motif instance mapping, distal|proximal asks if we should check the first or the second region in a pair. outputnetwork is just the output network to which we write the network out.

Example usage:
./genDistalProxNet testinput/rfNames_percentile_sorted_0.9_genenames.txt testinput/idrPool.wgEncodeOpenChromDnaseH1hescRep0.narrowPeak testinput/wgEncodeOpenChromDnaseH1hescRep0_FIMO_all.txt distal testoutput/h1esc_distalnet.txt


./genDistalProxNet testinput/rfNames_percentile_sorted_0.9_genenames.txt PIQresults distal testoutput/h1esc_distalnet.txt


PIQ_Outputs/TF_Instances/pancreas/MA05261.txt

[bbaur@roy-exec-4 cellTypeSpecificNetworks]$ head PIQ_Outputs/TF_Instances/pancreas/MA05261.txt
chr1	569849	0.943308550185874	+	MA05261
chr1	874012	0.827581484524788	+	MA05261
chr1	919843	0.772964947388806	+	MA05261
chr1	936302	0.912855377008653	+	MA05261

[bbaur@roy-exec-4 cd14_primary_cells]$ head MA07191.txt 
chr1	1491555	0.720034617048897	+	MA07191	MA07191
chr1	1507550	0.704840940525588	+	MA07191	MA07191
chr1	1550139	0.9375	+	MA07191	MA07191
chr1	1758395	0.71947194719472	+	MA07191	MA07191
chr1	1798448	0.708330175092852	+	MA07191	MA07191
chr1	1814472	0.720724191063174	+	MA07191	MA07191
chr1	2083714	0.730238816010764	+	MA07191	MA07191