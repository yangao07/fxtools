CC=			gcc
CFLAGS=		-g -Wall -O0
OBJS=		fxtools.o
PROG=		fxtools
LIB=		-lz

.SUFFIXES:.c .o

.c.o:
	$(CC) -c $(CFLAGS) $(MACRO) $< -o $@

all:$(PROG)
$(PROG):$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIB) -o $@

clean:
	rm -f *.o $(PROG)
