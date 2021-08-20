Compilation:
g++ *.C -I /mnt/ws/sysbio/roygroup/shared/thirdparty/gsl_install/include -L /mnt/ws/sysbio/roygroup/shared/thirdparty/gsl_install/lib  -lgsl -lgslcblas -g -o findSigIntr


Example usage:
./findSigIntr <intput file> <resolution> <radius> <scale> <output file>
./findSigIntr ../../findSigIntAnalysis/inputs/HUVEC_true_all.txt 5000 1000000 1 ../../findSigIntAnalysis/output/HUVEC_true_all.txt
./findSigIntr ../../findSigIntAnalysis/inputs/HUVEC_pred_all.txt 5000 1000000 1 ../../findSigIntAnalysis/output/HUVEC_pred_all.txt
./findSigIntr ../../findSigIntAnalysis/inputs/HUVEC_pred_all.txt 5000 1000000 7 ../../findSigIntAnalysis/output/HUVEC_pred_all_s7.txt

The second argument is the input file formatted as pair and a log count value as below 

==> ../../findSigIntAnalysis/inputs/HUVEC_pred_all.txt <==
chr1_100000000_100005000-chr1_100005000_100010000	3.207172
chr1_100000000_100005000-chr1_100010000_100015000	2.54412
chr1_100000000_100005000-chr1_100015000_100020000	2.27859
chr1_100000000_100005000-chr1_100020000_100025000	2.05304
chr1_100000000_100005000-chr1_100025000_100030000	1.7898

==> ../../findSigIntAnalysis/inputs/HUVEC_true_all.txt <==
chr1_100000000_100005000-chr1_100005000_100010000	2.99742
chr1_100000000_100005000-chr1_100010000_100015000	2.62953
chr1_100000000_100005000-chr1_100015000_100020000	1.79616
chr1_100000000_100005000-chr1_100020000_100025000	1.76211
chr1_100000000_100005000-chr1_100025000_100030000	1.81362

The second argument is the input file formatted as pair and count value (here not an integer).
The third and four arguments are the distance bin resolution (5000 bp) and selected radius (1000000 bp).
The fifth argument is a scaling factor to boost the pvalue (Set to 1 for true data, and 7 for predicted data in this case).
The last argument is the output file name.

The format of the output files is given below:

==> output/HUVEC_true_all.txt <==
Pair	Cnt	Dist	BinomQ	BinomP	GaussianQ	GaussianP
chr1_100000000_100005000-chr1_100005000_100010000	2.99742	0	0.670869	0.585381	0.941928	0.476147
chr1_10000000_10005000-chr1_10005000_10010000	3.10173	0	0.640003	0.410988	0.941928	0.392435
chr1_1000000_1005000-chr1_1005000_1010000	3.23794	0	0.640003	0.194566	0.941928	0.290708
chr1_100005000_100010000-chr1_100010000_100015000	2.89832	0	0.782474	0.749081	0.941928	0.556728

Here we have 
1. the pair id
2. the log of the counts for that pair 
3. The distancee between the two regions 
4. the q value from the Duan Binomial test for signficance 
5. the p value from the Duan binomial test 
6. the q value from the test of significance for a Gaussian model for the log counts for pairs at each distance bin
7. the p values from the Gaussian model signficance call. 

==> output/HUVEC_pred_all.txt <==
Pair	Cnt	BinomQ	BinomP	GaussianQ	GaussianP
chr1_100000000_100005000-chr1_100005000_100010000	3.20717	0	0.506955	0.100574	0.974769	0.175529
chr1_10000000_10005000-chr1_10005000_10010000	3.15	0	0.506955	0.197354	0.974769	0.235324
chr1_1000000_1005000-chr1_1005000_1010000	2.97828	0	0.510718	0.426441	0.974769	0.4652
chr1_100005000_100010000-chr1_100010000_100015000	3.18847	0	0.506955	0.143158	0.974769	0.193937

==> output/HUVEC_pred_all_s7.txt <==
Pair	Cnt	Dist	BinomQ	BinomP	GaussianQ	GaussianP
chr1_100000000_100005000-chr1_100005000_100010000	5.15308	0	0.328987	0.00226422	0.974769	0.175529
chr1_10000000_10005000-chr1_10005000_10010000	5.09591	0	0.764692	0.0206113	0.974769	0.235324
chr1_1000000_1005000-chr1_1005000_1010000	4.92419	0	0.764692	0.511536	0.974769	0.4652
chr1_100005000_100010000-chr1_100010000_100015000	5.13438	0	0.517888	0.00464457	0.974769	0.193937

Note that the output files is in a format such that if the scale input parameter  !=1, the the written out count values will be greater than the input values by log(<scale>).
