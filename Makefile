seq:
	gcc studentsseq.c -o studentsseq -lm
par:
	gcc studentspar.c -o studentspar -lm -fopenmp
	