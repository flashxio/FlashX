CC=g++
CXXFLAGS=-I../../include
LDFLAGS=-lnuma -lpthread -static -L../../libcache/ -lcache -lpthread -laio #-lhugetlbfs
CFLAGS=-O3

all: mem-read knem_memcpy rand-memcpy all2all-memtest

mem-read: mem-read.o mem-util.o
	$(CC) -o mem-read mem-read.o mem-util.o $(LDFLAGS)

rand-memcpy: rand-memcpy.o
	$(CC) -o rand-memcpy rand-memcpy.o $(LDFLAGS)

knem_memcpy: knem_memcpy.o
	$(CC) -o knem_memcpy knem_memcpy.o $(LDFLAGS)

all2all-memtest: all2all-memtest.o
	$(CC) -o all2all-memtest all2all-memtest.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f *~
	rm -f all2all-memtest calibrator knem_memcpy mem-read rand-memcpy stream test
