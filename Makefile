# settings
PROJECT = tmatrix
TARGET = $(PROJECT)_example
TEST = $(PROJECT)_test
CC = gcc
CFLAGS = -Wall
LDFLAGS = -lm -s
AR = ar
RM = rm -f 
.PHONY: all test

# sources 
OBJL = \
       src/lib/tmatrix.o \
       src/lib/tmatrix_calc.o \
       src/lib/tmatrix_homo.o \
       src/lib/tmatrix_vec.o \
       src/lib/tmatrix_io.o \
       src/lib/tmatrix_rot.o \
       src/lib/tmatrix_transform.o
OBJM = src/example.o
OBJT = src/tests/tests.o

OBJQ = src/transform.o

# main program
all: $(TARGET)

$(TARGET): $(OBJM) lib$(PROJECT).a
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

transforms: $(OBJQ) lib$(PROJECT).a
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@

test: $(TEST)

# build unit tests and run it
$(TEST): $(OBJT) lib$(PROJECT).a
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o $@
	./$@

lib$(PROJECT).a: $(OBJL)
	$(AR) rcs $@ $^

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

src/lib/tmatrix.o: src/lib/tmatrix.c
src/lib/tmatrix_calc.o: src/lib/tmatrix_calc.c
src/lib/tmatrix_homo.o: src/lib/tmatrix_homo.c
src/lib/tmatrix_vec.o: src/lib/tmatrix_vec.c
src/lib/tmatrix_io.o: src/lib/tmatrix_io.c
src/lib/tmatrix_rot.o: src/lib/tmatrix_rot.c
src/lib/tmatrix_transform.o: src/lib/tmatrix_transform.c
src/example.o: src/example.c
src/tests/tests.o: src/tests/tests.c

src/transform.o: src/transform.c

clean:
	$(RM) $(OBJL) $(OBJT) $(OBJM) lib$(PROJECT).a $(TEST) $(TARGET)
