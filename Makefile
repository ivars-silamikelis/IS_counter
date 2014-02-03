CC=gcc
CFLAGS=-Wall
DEPS = functions.h
OBJ = find_ends2.o functions.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

IS_find: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
