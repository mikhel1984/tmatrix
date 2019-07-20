# settings
TARGET = example
CFLAGS = -Wall
CC = gcc
LIBS = -lm

.PHONY: all test

# sources 
SRC = tmatrix.c tmatrix_calc.c tmatrix_homo.c tmatrix_vec.c
MAIN = example.c $(SRC)
TEST = tests/tests.c $(SRC)

# main program
all: $(MAIN)
	$(CC) $(CFLAGS) $(MAIN) -o $(TARGET) $(LIBS)

# build unit tests and run it
test: $(TEST)
	$(CC) $(CFLAGS) -g $(TEST) -o $@ $(LIBS)
	./test
