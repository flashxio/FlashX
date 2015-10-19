project (FlashGraph)

add_library(common STATIC
	common.cpp
	config_map.cpp
	log.cpp
	mem_tracker.cpp
	slab_allocator.cpp
	thread.cpp
	NUMA_mapper.cpp
)
