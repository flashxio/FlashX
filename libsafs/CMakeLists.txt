project (FlashGraph)

add_library(safs STATIC
	aio_private.cpp
	debugger.cpp
	io_c_interface.cpp
	messaging.cpp
	read_private.cpp
	thread.cpp
	associative_cache.cpp
	direct_private.cpp
	io_interface.cpp
	native_file.cpp
	remote_access.cpp
	timer.cpp
	cache_config.cpp
	disk_read_thread.cpp
	io_request.cpp
	parameters.cpp
	safs_file.cpp
	virt_aio_ctx.cpp
	cache.cpp
	file_mapper.cpp
	memory_manager.cpp
	part_global_cached_private.cpp
	shadow_cell.cpp
	workload.cpp
	common.cpp
	global_cached_private.cpp
	mem_tracker.cpp
	RAID_config.cpp
	slab_allocator.cpp
	wpaio.cpp
)
