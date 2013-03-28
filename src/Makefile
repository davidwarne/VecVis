#!/bin/make
#
# Summary: Really simple Makefile
#
#

#Setting for using Intel compiler
#CC = icc  
#OPTS = -O2 -fno-alias -xsse4.2 -vec-report1 -DNOISE=64 -DNSIZE=512 -DNSPOT=1000
#LIBS = -limf -lGLU -lglut -lg -lGL

#Setting for gcc
CC = gcc
OPTS = -O2 
LIBS = -lm -lGLU -lglut -lg -lGL

SRCDIR = ./src
BINDIR = ./bin
SRC = $(SRCDIR)/VecVis.c \
	$(SRCDIR)/Flw/FlwReader.c \
	$(SRCDIR)/BitMap/BitMapReader.c \
	$(SRCDIR)/BitMap/BitMapWriter.c \
	$(SRCDIR)/BitMap/BitMapFile.c

EXE = $(BINDIR)/VecVis
OBJS = $(SRC:.c=.o)
INC =
PROFILE = 

.SUFFIXES: .c .o

.c.o:
	$(CC) $(OPTS) -c $< -o $@

$(EXE): $(OBJS)
	$(CC) $(OPTS) $(OBJS) -o $(EXE) $(LIBS)
	@echo Binary created!!

clean:
	set nonomatch; rm -f $(EXE)
	rm -f $(OBJS)

