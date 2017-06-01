now=$(date +%Y%m%d)				##necessary for comparison of date line in vcf files
sed -i "2s/.*/##fileDate=$now/" *.vcf		##replace date line in (older) vcf output with current date. Otherwise tests will fail due to different dates.

python run_tests.py ../../.. ../../../../seqan-build/release/
