all:	pzip punzip
pzip:	pzip.c
	gcc -O pzip.c -o pzip
punzip:	punzip.c
	gcc -O punzip.c -o punzip
