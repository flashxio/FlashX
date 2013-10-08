#include "RAID_config.h"
#include "file_mapper.h"

file_mapper *RAID_config::create_file_mapper(const std::string &file_name) const
{
	std::vector<file_info> files = root_paths;
	for (unsigned i = 0; i < files.size(); i++) {
		files[i].name += std::string("/") + file_name;
	}
	switch (RAID_mapping_option) {
		case RAID0:
			return new RAID0_mapper(files, RAID_block_size);
		case RAID5:
			return new RAID5_mapper(files, RAID_block_size);
		case HASH:
			return new hash_mapper(files, RAID_block_size);
		default:
			fprintf(stderr, "wrong RAID mapping option\n");
			exit(1);
	}
}

std::set<int> RAID_config::get_node_ids() const
{
	std::set<int> node_ids;
	int num_paths = root_paths.size();
	for (int k = 0; k < num_paths; k++) {
		node_ids.insert(root_paths[k].node_id);
	}
	return node_ids;
}
