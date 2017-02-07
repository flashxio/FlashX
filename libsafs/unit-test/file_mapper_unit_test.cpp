#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include <memory>
#include <vector>

#include "file_mapper.h"

using namespace safs;

struct extended_block_identifier
{
	struct block_identifier bid;
	off_t orig_off;
};

typedef std::vector<struct extended_block_identifier> id_vec;

void test_get_file_sizes(const file_mapper &mapper, size_t block_size)
{
	for (size_t file_size = 0; file_size < 32 * 1024; file_size++) {
		std::vector<size_t> sizes = mapper.get_size_per_disk(file_size);
		for (size_t off = 1; off < file_size; off += block_size) {
			block_identifier bid;
			mapper.map(off, bid);
			assert((size_t) bid.off < sizes[bid.idx]);
		}
	}
}

int main()
{
	int num_files = 18;
	const int BLOCK_SIZE = 16;
	std::vector<part_file_info> files(num_files);

	srandom(time(NULL));

	printf("RAID0 mapper\n");
	RAID0_mapper mapper0("", files, BLOCK_SIZE);
	test_get_file_sizes(mapper0, BLOCK_SIZE);
	std::unique_ptr<id_vec[]> locs0
		= std::unique_ptr<id_vec[]>(new id_vec[num_files]);
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

	printf("RAID5 mapper\n");
	RAID5_mapper mapper5("", files, BLOCK_SIZE);
	test_get_file_sizes(mapper5, BLOCK_SIZE);
	std::unique_ptr<id_vec[]> locs5
		= std::unique_ptr<id_vec[]>(new id_vec[num_files]);
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

	printf("hash mapper\n");
	hash_mapper mapperh("", files, BLOCK_SIZE);
	test_get_file_sizes(mapperh, BLOCK_SIZE);
}
