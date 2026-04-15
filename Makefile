# compiler and flags
CC = gcc
CFLAGS = -Wall -lfftw3 -lm -g -O0

TARGET = output

FILENAME = "hellothere.xyz"

OMEGA ?= 6.0

SRC = waveletanalysis_overlap_save.c

all: build run

build:
	$(CC)  $(SRC) -o $(TARGET) $(CFLAGS)

run:
	./$(TARGET) $(FILENAME) $(OMEGA)

clean:
	rm -f $(TARGET)
