ifdef MEMCHECK
TRACE_FLAGS = -faddress-sanitizer
CC = clang
CXX = clang++
else
CC = gcc
CXX = g++
endif
TRACE_FLAGS += -fno-omit-frame-pointer # for better stack traces in error messages
TRACE_FLAGS += -fno-optimize-sibling-calls # disable tail call elimination
CXXFLAGS = -g -O3 -I.. -I../include -I../nbds.0.4.3/include/ -std=c++0x $(TRACE_FLAGS)
CFLAGS = -g -O3 -I.. -I../include $(TRACE_FLAGS)

LDFLAGS = -L../libcache -lcache -lpthread -lnuma -laio $(TRACE_FLAGS) -lprofiler -lrt

all: stat create_file create_files workload-gen memory-fill eval_expand_SA_cache verify_file print_file

stat: stat.o
	$(CXX) -o stat stat.o $(LDFLAGS)

create_file: create_file.o
	$(CXX) -o create_file create_file.o $(LDFLAGS)

create_files: create_files.o ../libcache/libcache.a
	$(CXX) -o create_files create_files.o $(LDFLAGS)

workload-gen: workload-gen.o
	$(CXX) -o workload-gen workload-gen.o $(LDFLAGS)

memory-fill: memory-fill.o
	$(CXX) -o memory-fill memory-fill.o $(LDFLAGS)

java_dump2c_dump: java_dump2c_dump.o
	$(CXX) -o java_dump2c_dump java_dump2c_dump.o ../workload.o $(LDFLAGS)

cache_evaluator: cache_evaluator.o ../libcache/libcache.a
	$(CXX) -o cache_evaluator cache_evaluator.o ../workload.o $(LDFLAGS) -L../nbds.0.4.3/ -lnbds

eval_expand_SA_cache: eval_expand_SA_cache.o ../libcache/libcache.a
	$(CXX) -o eval_expand_SA_cache eval_expand_SA_cache.o ../workload.o $(LDFLAGS)

verify_file: verify_file.o
	$(CXX) -o verify_file verify_file.o $(LDFLAGS)

print_file: print_file.o ../libcache/libcache.a
	$(CXX) -o print_file print_file.o $(LDFLAGS)

clean:
	rm *.o
