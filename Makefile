seq:
	gcc studentsseq.c -o studentsseq -lm -fopenmp
par:
	gcc studentspar.c -o studentspar -lm -fopenmp
	