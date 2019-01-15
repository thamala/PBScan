# PBScan
Selection outlier scan with population branch statistic

Compiling:
gcc pbscan.c -lm -o pbscan

For help and options:
./pbscan

Usage example (files can be found in the example_data folder):
./pbscan -likes example.likes -pop1 list1.txt -pop2 list2.txt -pop3 list3.txt -out output -win 50 -step 51 -ms example.ms -div 2 -min 5 -maf 0.05

See the manual for more information
