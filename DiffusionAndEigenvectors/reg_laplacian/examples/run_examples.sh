# Runs an example for a few values of lambda.
# Choose example from chain, star, degree, bg
# Overwrites existing output files!
# Uses temps file named "example_temp", "example_temp_b", "example_temp_c" 
# to create some nice output for the terminal. 

if [ $# -ne 2 ]; then
	echo "Please supply example name to run (chain, star, degree, bg) and output location"
	exit
fi

EX=$1
OUT=$2

UTIL=/mnt/ws/sysbio/roygroup/shared/programs/kernels_on_graphs/reg_laplacian/get_kernel_scores.py

if [ ! -d $OUT ]; then
	mkdir $OUT
	echo "Made directory $OUT"
fi

TEMP=example_temp
echo "Writing output to files ${EX}_\$lambda_example_output.tab for various values of lambda, also cross-validation"
echo "Also writing various things to ${TEMP}_a, b, and c"

for lam in 1 #0.1 0.5 1 5 10 50 100; do
do
	if [ -a "tiny_${EX}_kernel_${lam}.npy" ]; then
		python $UTIL tiny_${EX}_scores.tab tiny_${EX}.sif $lam tiny_${EX}_kernel_${lam}.npy --crossval ${OUT}/${EX}_${lam}_example_output.tab	
	else
		python $UTIL tiny_${EX}_scores.tab tiny_${EX}.sif $lam --crossval ${OUT}/${EX}_${lam}_example_output.tab
	fi
done

# print only full scores from each lamda and header to terminal
cut -f 1,2,3 ${OUT}/${EX}_0.1_example_output.tab | grep -v "#" > ${TEMP}	
for lam in 1 #0.5 1 5 10 50 100; do
do
	cut -f 3 ${OUT}/${EX}_${lam}_example_output.tab | grep -v "#" > ${TEMP}_b
	paste ${TEMP} ${TEMP}_b > ${TEMP}_c
	mv ${TEMP}_c ${TEMP}
done

cat ${TEMP}

rm ${TEMP}*
