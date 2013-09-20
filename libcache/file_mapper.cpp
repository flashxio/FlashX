#include "file_mapper.h"

int gen_RAID_rand_start(int num_files)
{
	int r = random();
	while (r == 0) {
		r = random();
	}
	return r % num_files;
}

int RAID0_mapper::rand_start;
int RAID5_mapper::rand_start;

atomic_integer file_mapper::file_id_gen;
