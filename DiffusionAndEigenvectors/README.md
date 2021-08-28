After performing the gene scoring (code provided in mapregionpairsnp), we next performed node and edge diffusion on the graph. 

Node diffusion: We did node diffusion using a regularized laplacian kernel (see readme in reg_laplacian).
Edge diffusion was performed with the hotnet2 pipeline: https://github.com/raphael-group/hotnet2

Eigenvectors were calculated much in the same way as the original MUSCARI pipeline (https://github.com/Roy-lab/Muscari) the only difference being we used the L1 norm instead of the L2 norm. 
