USAGE:
$features - feature file for pairs generated by the genFeaturesForDiscreteData code provided in this directory.

$start - genomic position of the 1Mb region start

$end - genomic position of the 1Mb end

$num - total number of features (including distance, here it is 31)

$n - set to 1 

$out - output prefix

the output is the same as $features just for all the features from $start to $end (either or both regions in the pair must be between $start and $end

./getPairsForRegion/genPairs_Features ${features} ${start} ${end} ${num} ${end} ${n} ${out}

