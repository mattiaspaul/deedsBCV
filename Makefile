HEADERS=
SOURCES=src/deedsBCV0.cpp

ifeq ($(SLOW),1)
	OPT =-O
else
	OPT =-O3 -fopenmp -mavx2 -msse4.2
endif

all: deeds

deeds: $(HEADERS) $(SOURCES) Makefile
	g++ $(SOURCES) -I src -lz -o deedsBCV -std=c++11 $(OPT)

clean:
	rm -f deedsBCV

