#include <assert.h>

#include <vector>

#include "io_interface.h"

static std::vector<io_interface *> io_table;

void register_io(io_interface *io)
{
	// Make sure the index hasn't been set.
	assert(io->get_io_idx() < 0);

	io_table.push_back(io);
	int idx = io_table.size() - 1;
	io->set_io_idx(idx);
}

io_interface *get_io(int idx)
{
	assert(idx >= 0);
	io_interface *io = io_table[idx];
	assert(io->get_io_idx() == idx);
	return io;
}
