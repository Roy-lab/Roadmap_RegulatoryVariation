celline='pancreas'
phenotype='Schizophrenia'
matrix=sprintf('/mnt/dv/wid/projects3/Roy-enhancer-promoter/RMEC/Graph_Diffusion/Results/Graph_Diffusion_Results_NewScores/%s/%s/GraphDiff_1_UnionNorm_top5_adj.txt.txt',celline,phenotype);
	%matrix=sprintf('DSD_Reduced_Graphs/%s/%s_MatForm.txt',phenotype,cellines{i})
genes=importdata(sprintf('../../Gene_ID_Files_top5/%s/%s.txt',phenotype,celline));
strout=sprintf('/mnt/dv/wid/projects3/Roy-enhancer-promoter/RMEC/Graph_Diffusion/Results/Arboretum_Eigs_NewScores_allEigs_qnormUnion_L1/%s/%s/muscari_eigs_qnorm_test.txt',celline,phenotype)	
%strout=sprintf('/mnt/dv/wid/projects3/Roy-enhancer-promoter/RMEC/Graph_Diffusion/Results/Arboretum_Eigs_NewScores_allEigs_qnormUnion_L1/%s/%s/muscari_eigs_qnorm_test.txt',cellines,phenotype)

eigvecmat_calc_L1(matrix,genes,strout);

	%phenotype='PostbronchodilatorFEV1FVCratioinCOPD'
       %matrix=sprintf('/mnt/dv/wid/projects3/Roy-enhancer-promoter/RMEC/Graph_Diffusion/Results/Graph_Diffusion_Results_NewScores/%s/%s/GraphDiff_1_IntNorm_adj.txt.txt','h1_cells',phenotype);       
%genes=importdata(sprintf('../../Gene_ID_Files_Int/%s/%s.txt',phenotype,'h1_cells'));
%       strout=sprintf('/mnt/dv/wid/projects3/Roy-enhancer-promoter/RMEC/Graph_Diffusion/Results/Arboretum_Eigs_NewScores_allEigs_v2/%s/%s/muscari_eigs_qnorm_test.txt','h1_cells',phenotype)
 %     eigvecmat_calc(matrix,genes,strout);

