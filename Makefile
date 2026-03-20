# compiler and flags
CC = gcc
CFLAGS = -Wall -lfftw3 -lm -g -O0

# name of output
TARGET = output

# name of the .xyz file
FILENAME = "hellothere.xyz"
# stadard parameter for omega
OMEGA ?= 6.0

# source
SRC = waveletanalysis_overlap_save.c

# default target
all: build run

# compile
build:
	$(CC)  $(SRC) -o $(TARGET) $(CFLAGS)

# run
run:
	./$(TARGET) $(FILENAME) $(OMEGA)

# cleanup
clean:
	rm -f $(TARGET)
