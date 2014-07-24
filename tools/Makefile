include ../Makefile.common

CXXFLAGS += -I.. -I../include
CFLAGS = -g -O3 -I.. -I../include $(TRACE_FLAGS)

LDFLAGS := -L../libsafs -lsafs $(LDFLAGS)

TARGETS = create_file memory-fill print_file
LIBFILE = ../libsafs/libsafs.a

all: $(TARGETS)

stat: stat.o
	$(CXX) -o stat stat.o $(LDFLAGS)

create_file: create_file.o
	$(CXX) -o create_file create_file.o $(LDFLAGS)

workload-gen: workload-gen.o
	$(CXX) -o workload-gen workload-gen.o $(LDFLAGS)

memory-fill: memory-fill.o
	$(CXX) -o memory-fill memory-fill.o $(LDFLAGS)

java_dump2c_dump: java_dump2c_dump.o
	$(CXX) -o java_dump2c_dump java_dump2c_dump.o $(LDFLAGS)

cache_evaluator: cache_evaluator.o $(LIBFILE)
	$(CXX) -o cache_evaluator cache_evaluator.o $(LDFLAGS)

eval_expand_SA_cache: eval_expand_SA_cache.o $(LIBFILE)
	$(CXX) -o eval_expand_SA_cache eval_expand_SA_cache.o $(LDFLAGS)

print_file: print_file.o $(LIBFILE)
	$(CXX) -o print_file print_file.o $(LDFLAGS)

clean:
	rm -f *.o
	rm -f *~
	rm -f $(TARGETS)
	rm -f *.class
