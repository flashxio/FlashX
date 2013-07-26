#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sched.h>
#include <assert.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "common.h"

extern "C" {

void bind2cpu(int cpu_id)
{
	cpu_set_t cpu_mask;
	CPU_ZERO(&cpu_mask);
	CPU_SET(cpu_id, &cpu_mask);
	int len = sizeof(cpu_mask);
	if (sched_setaffinity(gettid(), len, &cpu_mask) < 0) {
		perror("sched_setaffinity");
		exit(1);
	}
}

void bind_mem2node_id(int node_id)
{
	struct bitmask *bmp = numa_allocate_nodemask();
	numa_bitmask_setbit(bmp, node_id);
	numa_set_membind(bmp);
	numa_free_nodemask(bmp);
}

void bind2node_id(int node_id)
{
	struct bitmask *bmp = numa_allocate_nodemask();
	numa_bitmask_setbit(bmp, node_id);
	numa_bind(bmp);
	numa_free_nodemask(bmp);
}

int numa_get_mem_node()
{
	struct bitmask *bmp = numa_get_membind();
	int nbytes = numa_bitmask_nbytes(bmp);
	int num_nodes = 0;
	int node_id = -1;
	int i;
	for (i = 0; i < nbytes * 8; i++)
		if (numa_bitmask_isbitset(bmp, i)) {
			num_nodes++;
			node_id = i;
		}
	assert(num_nodes == 1);
	return node_id;
}

ssize_t get_file_size(const char *file_name)
{
	struct stat stats;
	if (stat(file_name, &stats) < 0) {
		perror("stat");
		return -1;
	}
	return stats.st_size;
}

void permute_offsets(int num, int repeats, int stride, off_t start,
		off_t offsets[])
{
	int idx = 0;
	for (int k = 0; k < repeats; k++) {
		for (int i = 0; i < num; i++) {
			offsets[idx++] = ((off_t) i) * stride + start * stride;
		}
	}
	int tot_length = idx;

	for (int i = tot_length - 1; i >= 1; i--) {
		int j = random() % tot_length;
		off_t tmp = offsets[j];
		offsets[j] = offsets[i];
		offsets[i] = tmp;
	}
}
}

bool align_check(size_t alignment)
{
	assert(alignment >= 0);
	if (alignment == 0 || alignment == 1) {
		return false;
	}
	bool aligned = true;
	while (alignment > 1) {
		// If it's not a power of 2.
		if (alignment % 2) {
			aligned = false;
			break;
		}
		alignment /= 2;
	}
	return aligned;
}

long str2size(std::string str)
{
	int len = str.length();
	long multiply = 1;
	if (str[len - 1] == 'M' || str[len - 1] == 'm') {
		multiply *= 1024 * 1024;
		str[len - 1] = 0;
	}
	else if (str[len - 1] == 'K' || str[len - 1] == 'k') {
		multiply *= 1024;
		str[len - 1] = 0;
	}
	else if (str[len - 1] == 'G' || str[len - 1] == 'g') {
		multiply *= 1024 * 1024 * 1024;
		str[len - 1] = 0;
	}
	return atol(str.c_str()) * multiply;
}

int retrieve_data_files(std::string file_file,
		std::vector<file_info> &data_files)
{
	char *line = NULL;
	size_t size = 0;
	int line_length;
	FILE *fd = fopen(file_file.c_str(), "r");
	if (fd == NULL) {
		perror("fopen");
		return 0;
	}
	while ((line_length = getline(&line, &size, fd)) > 0) {
		line[line_length - 1] = 0;
		// skip comment lines.
		if (*line == '#')
			continue;

		char *colon = strstr(line, ":");
		file_info info;
		char *name = line;
		if (colon) {
			*colon = 0;
			info.node_id = atoi(line);
			colon++;
			name = colon;
		}
		info.name = name;
		data_files.push_back(info);
		free(line);
		line = NULL;
		size = 0;
	}
	fclose(fd);
	return data_files.size();
}

sys_parameters params;
