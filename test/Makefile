LDFLAGS = -lpthread -lnuma -laio $(TRACE_FLAGS) #-lprofiler
UNITTEST = SAHT_unit_test GClock_unit_test SA_expand_shrink_test
CPPFLAGS := -MD
CXXFLAGS = -I.. -g
SOURCE := $(wildcard *.c) $(wildcard *.cpp)
OBJS := $(patsubst %.c,%.o,$(patsubst %.cpp,%.o,$(SOURCE)))
DEPS := $(patsubst %.o,%.d,$(OBJS))

all: $(UNITTEST)

SAHT_unit_test: SAHT_unit_test.o ../SA_hash_table.o
	$(CXX) -o SAHT_unit_test SAHT_unit_test.o ../SA_hash_table.o $(LDFLAGS)

GClock_unit_test: GClock_unit_test.o ../gclock.o ../slab_allocator.o
	$(CXX) -o GClock_unit_test GClock_unit_test.o ../gclock.o ../slab_allocator.o $(LDFLAGS)

SA_expand_shrink_test: SA_expand_shrink_test.o ../associative_cache.o ../memory_manager.o ../slab_allocator.o ../messaging.o ../thread.o
	$(CXX) -o SA_expand_shrink_test SA_expand_shrink_test.o ../associative_cache.o ../memory_manager.o ../slab_allocator.o ../messaging.o ../thread.o $(LDFLAGS)

clean:
	rm *.o

-include $(DEPS) 
