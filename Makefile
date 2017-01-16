TARGET = main
LIBS = -lm
CC = g++
CFLAGS = -std=c++11 -Wall -Wextra -g -I./sais-lite-2.4.1/

.PHONY: default all clean

default: $(TARGET)
all: default

OBJECTS = sais.o main.o
HEADERS = /sais-lite-2.4.1/sais.h

sais.o: sais-lite-2.4.1/sais.c
	$(CC) $(CFLAGS) -c $< -o $@

main.o: src/main.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Wall $(LIBS) -o $@

clean:
	-rm -f *.o
	-rm -f $(TARGET)