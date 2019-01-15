# PBScan
Selection outlier scan with population branch statistic

Compiling:
gcc pbscan.c -lm -o pbscan

For help and options:
./pbscan

Usage example:
./pbscan -likes example.likes -pop1 list1.txt -pop2 list2.txt -pop3 list3.txt -out output -win 50 -step 51 -ms example.ms -div 2 -min 5 -maf 0.05

PBScan is run on genotype likelihoods (-likes), with -pop1 and -pop2 defining individuals belonging to the two focal populations and -pop3 defining the outgroup population. Analysis is done in 50 SNP non-overlapping sliding windows (-win 50 and -step 51), using neutral data generated with <i>ms</i> (-ms) in outlier detection. Population differentiation is estimated with <i>d</i><sub>XY</sub> (-div 2). Sites with less than 5 individuals per population (-min 5) and minor allele frequency less than 0.05 (-maf 0.05) are further excluded from the analysis. Files can be found in the example_data folder.

See the manual for more information
