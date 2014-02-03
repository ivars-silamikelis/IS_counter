CC=gcc
CFLAGS=-Wall

IS_find: find_ends2.o functions.o
	$(CC) -o find_ends find_ends2.o functions.o
