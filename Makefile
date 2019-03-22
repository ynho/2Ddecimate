# Makefile by Yno, all rights freeds :D

EXEC=a
OPTI=-g
PKGS=SDL2

CC ?= gcc
CFLAGS = -Wall $(OPTI) `pkg-config $(PKGS) --cflags` --std=c99 -g
#LDFLAGS = `pkg-config $(PKGS) --libs` -lGL -lm
LDFLAGS = -lGL -lm -lSDL2

all: $(EXEC)

$(EXEC): decimate.o
	$(CC) $(LDFLAGS) $^ -o $@

decimate.o: decimate.c Makefile
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o

mrproper: clean
	rm -f $(EXEC)

.PHONY: all clean mrproper
