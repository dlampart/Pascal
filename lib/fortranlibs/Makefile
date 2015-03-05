CC=gfortran
SOURCES=mvtdstpack.f 
OBJECTS=$(SOURCES:.f=.o)
LDFLAGS=-Wall -shared
CFLAGS=-Wall -fPIC -frecursive
UNAME_S := $(shell uname)
print_vars:
	echo $(UNAME_S)

all: libmvtpack.dylib

libmvtpack.dylib: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

libmvtpack.so: $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

libmvnpack.dlyb: mvndstpack.o
	$(CC) $(LDFLAGS) mvndstpack.o -o $@

.f.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o libmvnpack.dylib libmvtpack.dylib libmvtpack.so

