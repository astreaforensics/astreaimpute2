CC=gcc
CFLAGS=-gdwarf-2 -ggdb
#CFLAGS=-O2
#CFLAGS=-gdwarf-2 -g

# (c) 2024 - Ed Green
#   Astrea Forensics

OBJ=local-impute.o file-io.o sam-parse.o pileup.o hap.o
EXE=score-filter-sam astrea-impute2.14 test-pileup test-hap hap2af-table

all: $(EXE)

clean:
	rm -f $(OBJ) $(EXE) t/test_out/*

test: t/test.sh
	./t/test.sh

stats.o : stats.h stats.c
	echo "Making stats.o..."
	$(CC) $(CFLAGS) -c stats.c -o stats.o 
hap.o : hap.h hap.c
	echo "Making hap.o..."
	$(CC) $(CFLAGS) -c hap.c -lz -o hap.o

pileup.o : pileup.h pileup.c
	echo "Making pileup.o..."
	$(CC) $(CFLAGS) -c pileup.c -lz -o pileup.o

file-io.o : file-io.h file-io.c
	echo "Making file-io.o..."
	$(CC) $(CFLAGS) -c file-io.c -lz -o file-io.o

local-impute.o : local-impute.h local-impute.c file-io.o
	echo "Making local-impute.o..."
	$(CC) $(CFLAGS) -c local-impute.c -o local-impute.o -lm

sam-parse.o : sam-parse.h sam-parse.c
	echo "Making sam-parse.o..."
	$(CC) $(CFLAGS) -c sam-parse.c -o sam-parse.o -lm

score-filter-sam : score-filter-sam.c sam-parse.o
	echo "Making score-filter-sam..."
	$(CC) $(CFLAGS) sam-parse.o -o score-filter-sam score-filter-sam.c -lm

test-local-impute1 : test-local-impute1.c local-impute.o
	echo "Making test-local-impute1..."
	$(CC) $(CFLAGS) local-impute.o -o test-local-impute1 test-local-impute1.c -lm

test-pileup : test-pileup.c pileup.o file-io.o
	echo "Making test-pileup..."
	$(CC) $(CFLAGS) pileup.o file-io.o -lz -o test-pileup test-pileup.c -lm

test-hap : test-hap.c hap.o file-io.o
	echo "Making test-hap..."
	$(CC) $(CFLAGS) hap.o file-io.o -lz -o test-hap test-hap.c -lm

test-stats : test-stats.c stats.o
	echo "Making test-stats..."
	$(CC) $(CFLAGS) stats.o -o test-stats test-stats.c -lm

astrea-impute2.14 : astrea-impute2.14.c hap.o file-io.o pileup.o stats.o
	echo "Making astrea-impute2..."
	$(CC) $(CFLAGS) hap.o pileup.o file-io.o stats.o -lz -o astrea-impute2.14 astrea-impute2.14.c -lm -pthread

hap2af-table : hap2af-table.c file-io.o hap.o
	echo "Making hap2af-table..."
	$(CC) $(CFLAGS) hap.o file-io.o -lz -o hap2af-table hap2af-table.c -lm

