# compiler and flags
CC = gcc
CFLAGS = -Wall -lfftw3 -lm -g -O0 -DFD

# name of output
TARGET = output

# stadard parameter for omega
OMEGA ?= 15
MODE ?= FD

# source
SRC = waveletanalysis_padded.c

# default target
all: build run

# compile
build:
	$(CC)  $(SRC) -o $(TARGET) $(CFLAGS) $(addprefix -D, $(MODE))

# run
run:
	./$(TARGET) $(OMEGA)

# cleanup
clean:
	rm -f $(TARGET)
