# the compiler: gcc for C program, define as g++ for C++
CC = g++

# compiler flags:
#-gadds debugging information to the executable file
#-Wall turns on most, but not all, compiler warnings
CFLAGS= -g  -L/usr/lib/ -lm #-Wall 
GSLFLAGS=-lgsl -lgslcblas 

all: genPairs_Features

genPairs_Features: Framework.C
	$(CC) $(CFLAGS) -o genPairs_Features *.C $(GSLFLAGS) 

clean:
	$(RM) genPairs_Features
