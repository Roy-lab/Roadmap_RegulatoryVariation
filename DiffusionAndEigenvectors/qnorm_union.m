cells=importdata('./orders.txt');
phenos={'Inflammatoryskindisease','Obesity-relatedtraits','Pediatricautoimmunediseases','PostbronchodilatorFEV1','PostbronchodilatorFEV1FVCratio','PostbronchodilatorFEV1FVCratioinCOPD','Restingheartrate','Schizophrenia','Type2diabetes'};
%phenos={'Bloodproteinlevels','Breastcancer','Crohnsdisease','Height','IgGglycosylation','Inflammatoryboweldisease'}
%phenos={'Autismspectrumdisorderorschizophrenia'}
phenos={'Bodymassindex'}
for i=1:length(phenos)
pheno=phenos{i}
adj_genes={};
all_genes=importdata(sprintf('/mnt/dv/wid/projects3/Roy-enhancer-promoter/RMEC/Graph_Diffusion/Results/snpgenemappings/%s/union_top1percent_uniq_top5.txt',pheno));
concat_mats=[];
lens=[];
for c=1:length(cells)
cell=cells{c}
mat=importdata(sprintf('/mnt/dv/wid/projects3/Roy-enhancer-promoter/RMEC/Graph_Diffusion/Results/Graph_Diffusion_Results_NewScores/%s/%s/GraphDiff_1_Union_top5_adj.txt',cell,pheno));
cell_genes=readtable(sprintf('/mnt/dv/wid/projects3/Roy-enhancer-promoter/RMEC/Graph_Diffusion/Results/snpgenemappings/%s/%s/idx_top1percent_top5.txt',pheno,cell),'ReadVariableNames',false,'Delimiter','\t');
cmat=zeros(length(cell_genes.Var2),length(all_genes));
for k=1:length(cell_genes.Var2)
idx=strmatch(cell_genes.Var2{k},all_genes,'exact');
cmat(:,idx)=mat(:,k);
end
concat_mats=[concat_mats;cmat];
lens=[lens,length(cell_genes.Var2)];
end

NormData = quantilenorm(concat_mats);

sizes=1:size(mat,1):(size(concat_mats,2)+size(mat,1))
for j=1:55
cell=cells{j};
cell_genes=readtable(sprintf('/mnt/dv/wid/projects3/Roy-enhancer-promoter/RMEC/Graph_Diffusion/Results/snpgenemappings/%s/%s/idx_top1percent_top5.txt',pheno,cell),'ReadVariableNames',false,'Delimiter','\t');
if j==1
c_norm=NormData(1:lens(j),:);
else
c_norm=NormData(sum(lens(1:j-1))+1:sum(lens(1:j)),:);
end
[ia,ib]=ismember(cell_genes.Var2,all_genes);
c_norm2=c_norm(:,ib);
dlmwrite(sprintf('/mnt/dv/wid/projects3/Roy-enhancer-promoter/RMEC/Graph_Diffusion/Results/Graph_Diffusion_Results_NewScores/%s/%s/GraphDiff_1_UnionNorm_top5_adj.txt.txt',cells{j},pheno),c_norm2);
end
end
