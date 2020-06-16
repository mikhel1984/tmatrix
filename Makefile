# settings
TARGET = example
CFLAGS = -Wall
CC = gcc
LIBS = -lm

.PHONY: all test

# sources 
SRC = lib/tmatrix.c lib/tmatrix_calc.c lib/tmatrix_homo.c lib/tmatrix_vec.c \
      lib/tmatrix_io.c
MAIN = example.c $(SRC)
TEST = tests/tests.c $(SRC)

# main program
all: $(MAIN)
	$(CC) $(CFLAGS) $(MAIN) -o $(TARGET) $(LIBS)

# build unit tests and run it
test: $(TEST)
	$(CC) $(CFLAGS) -g $(TEST) -o $@ $(LIBS)
	./test
