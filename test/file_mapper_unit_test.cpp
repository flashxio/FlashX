#include <assert.h>

#include <vector>

#include "../file_mapper.h"

struct extended_block_identifier
{
	struct block_identifier bid;
	off_t orig_off;
};

int main()
{
	int num_files = 16;
	std::vector<struct extended_block_identifier> locs[num_files];
	std::vector<file_info> files(num_files);
	RAID0_mapper mapper(files, 16);
	for (int i = 0; i < 100000; i++) {
		block_identifier bid;
		mapper.map(i, bid);
		struct extended_block_identifier ebid;
		ebid.bid = bid;
		ebid.orig_off = i;
		locs[bid.idx].push_back(ebid);
	}
	for (int i = 0; i < num_files; i++)
		printf("file %d: has %ld blocks\n", i, locs[i].size());

	for (int i = 0; i < num_files; i++) {
		struct extended_block_identifier prev = locs[i][0];
		for (size_t j = 1; j < locs[i].size(); j++) {
			assert(prev.orig_off + num_files == locs[i][j].orig_off);
			prev = locs[i][j];
		}
	}
}
