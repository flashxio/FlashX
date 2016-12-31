project (FlashGraph)

add_library(safs STATIC
	aio_private.cpp
	debugger.cpp
	messaging.cpp
	read_private.cpp
	associative_cache.cpp
	direct_private.cpp
	io_interface.cpp
	native_file.cpp
	remote_access.cpp
	cache_config.cpp
	disk_read_thread.cpp
	io_request.cpp
	parameters.cpp
	safs_file.cpp
	cache.cpp
	file_mapper.cpp
	memory_manager.cpp
	part_global_cached_private.cpp
	shadow_cell.cpp
	global_cached_private.cpp
	RAID_config.cpp
	wpaio.cpp
	direct_comp_access.cpp
	comp_io_scheduler.cpp
	in_mem_io.cpp
	NUMA_mapper.cpp
	common.cpp
	config_map.cpp
	log.cpp
	mem_tracker.cpp
	slab_allocator.cpp
	thread.cpp
)
