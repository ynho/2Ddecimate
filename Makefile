EXEC=a
OPTI=-g
PKGS=

CC?=gcc
#CFLAGS = -Wall $(OPTI) `pkg-config $(PKGS) --cflags` --std=c99 -g
CFLAGS=-Wall $(OPTI) --std=c99 -g
#LDFLAGS = `pkg-config $(PKGS) --libs` -lGL -lm
LDFLAGS=-lGL -lm -lSDL2

all: $(EXEC)

$(EXEC): decimate.o
	$(CC) $^ -o $@ $(LDFLAGS)

decimate.o: decimate.c Makefile
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o

mrproper: clean
	rm -f $(EXEC)

.PHONY: all clean mrproper
