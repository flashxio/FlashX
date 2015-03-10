#include "matrix_io.h"
#include "matrix_config.h"

using namespace fm;

struct comp_io
{
	bool operator()(const matrix_io &io1, const matrix_io &io2) const {
		return io1.get_top_left().get_row_idx()
			< io2.get_top_left().get_row_idx();
	}
};

void run_test()
{
	size_t row_block_size = 1024;
	size_t num_rows = random() + 1000000000L;
	size_t num_blocks = ceil(((double) num_rows) / row_block_size);
	std::vector<row_block> blocks;
	off_t off = 0;
	for (size_t i = 0; i < num_blocks; i++) {
		blocks.emplace_back(off);
		off += random() % (1024 * 1024);
	}
	blocks.emplace_back(off);

	int num_gens = 32;
	std::vector<matrix_io> ios;
	std::vector<matrix_io_generator::ptr> io_gens(num_gens);
	for (int i = 0; i < num_gens; i++) {
		row_block_mapper mapper(blocks, i, num_gens,
				matrix_conf.get_rb_io_size());
		io_gens[i] = matrix_io_generator::create(blocks,
				num_rows, 1000, 0, mapper);
	}
	while (!io_gens.empty()) {
		int idx = random() % io_gens.size();
		if (io_gens[idx]->has_next_io()) {
			if (random() % 10 == 0)
				ios.push_back(io_gens[idx]->steal_io());
			else
				ios.push_back(io_gens[idx]->get_next_io());
		}
		else
			io_gens.erase(io_gens.begin() + idx);
	}
	std::sort(ios.begin(), ios.end(), comp_io());
	printf("There are %ld ios\n", ios.size());
	off_t prev_off = ios[0].get_loc().get_offset();
	size_t prev_size = ios[0].get_size();
	off_t prev_row = ios[0].get_top_left().get_row_idx();
	size_t prev_num_rows = ios[0].get_num_rows();
	assert(prev_off == 0);
	assert(prev_row == 0);
	for (size_t i = 1; i < ios.size() - 1; i++) {
		assert(ios[i].get_loc().get_offset() == prev_off + prev_size);
		assert(ios[i].get_top_left().get_row_idx() == prev_row + prev_num_rows);
		prev_off = ios[i].get_loc().get_offset();
		prev_size = ios[i].get_size();
		prev_row = ios[i].get_top_left().get_row_idx();
		prev_num_rows = ios[i].get_num_rows();
	}
	assert(ios.back().get_loc().get_offset() + ios.back().get_size()
			== blocks.back().get_offset());
	assert(ios.back().get_top_left().get_row_idx()
			+ ios.back().get_num_rows() == num_rows);
}

int main()
{
	for (int test = 0; test < 100; test++)
		run_test();
}
