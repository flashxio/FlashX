#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include <vector>

#include "file_mapper.h"

struct extended_block_identifier
{
	struct block_identifier bid;
	off_t orig_off;
};

int main()
{
	int num_files = 18;
	const int BLOCK_SIZE = 16;
	std::vector<part_file_info> files(num_files);

	srandom(time(NULL));

	RAID0_mapper mapper0("", files, BLOCK_SIZE);
	std::vector<struct extended_block_identifier> locs0[num_files];
	printf("RAID0 mapper\n");
	for (int i = 0; i < 10000; i++) {
		off_t off = i * BLOCK_SIZE;
		block_identifier bid;
		mapper0.map(off, bid);
		struct extended_block_identifier ebid;
		ebid.bid = bid;
		assert(bid.idx == mapper0.map2file(off));
		ebid.orig_off = off;
		locs0[bid.idx].push_back(ebid);
	}
	for (int i = 0; i < num_files; i++)
		printf("file %d: has %ld blocks\n", i, locs0[i].size());

	for (int i = 0; i < num_files; i++) {
		struct extended_block_identifier prev = locs0[i][0];
		printf("file %d\n", i);
		for (size_t j = 1; j < locs0[i].size(); j++) {
			assert(prev.bid.off + BLOCK_SIZE == locs0[i][j].bid.off);
			assert(prev.orig_off + BLOCK_SIZE * num_files == locs0[i][j].orig_off);
			prev = locs0[i][j];
		}
	}

	RAID5_mapper mapper5("", files, BLOCK_SIZE);
	std::vector<struct extended_block_identifier> locs5[num_files];
	printf("RAID5 mapper\n");
	for (int i = 0; i < 32; i++) {
		off_t off = i * BLOCK_SIZE;
		block_identifier bid;
		mapper5.map(off, bid);
		struct extended_block_identifier ebid;
		ebid.bid = bid;
		assert(bid.idx == mapper5.map2file(off));
		ebid.orig_off = off;
		printf("off: %ld, idx in stripe: %d, stripe off: %ld\n", off, bid.idx, bid.off);
		locs5[bid.idx].push_back(ebid);
	}
	for (int i = 0; i < num_files; i++)
		printf("file %d: has %ld blocks\n", i, locs5[i].size());

	for (int i = 0; i < num_files; i++) {
		struct extended_block_identifier prev = locs5[i][0];
		printf("file %d\n", i);
		for (size_t j = 1; j < locs5[i].size(); j++) {
			assert(prev.bid.off + BLOCK_SIZE == locs5[i][j].bid.off);
			assert((locs5[i][j].orig_off / BLOCK_SIZE - prev.orig_off / BLOCK_SIZE) % num_files == num_files - 1);
			prev = locs5[i][j];
		}
	}
}
