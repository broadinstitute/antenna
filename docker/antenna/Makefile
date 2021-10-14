test: FORCE
	python3 -O antenna.py --bam test/subsample2_aligned.bam --orf-bed test/orf_start.bed --output-counts output_test.csv

test_short: FORCE
	python3 antenna.py --bam test/aligned.bam --orf-bed test/orf_start.bed --output-counts output_test_short.csv
	cat output_test_short.csv | column -t -s,

test_big: FORCE
	python3 antenna.py --bam test/big_aligned.bam --orf-bed test/orf_start.bed --output-counts output_test_big.csv

profile_short: FORCE
	python3 -m cProfile antenna.py --bam test/aligned.bam --orf-bed test/orf_start.bed --output-counts output.csv

lint:
	black .

FORCE:
