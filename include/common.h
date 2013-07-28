#ifndef __MY_COMMON_H__
#define __MY_COMMON_H__

#include <string>
#include <vector>

#include "common_c.h"
#include "parameters.h"

enum {
	READ,
	WRITE
};

struct file_info
{
	std::string name;
	// The NUMA node id where the disk is connected to.
	int node_id;
};

/**
 * Check if the integer is a power of two.
 */
bool align_check(size_t alignment);

int retrieve_data_files(std::string file_file,
		std::vector<file_info> &data_files);

static inline std::string itoa(int n)
{
	char buf[32];
	snprintf(buf, sizeof(buf), "%d", n);
	return buf;
}

template<class T>
void rand_permute_array(T arr[], int num)
{
	for (int i = num - 1; i >= 1; i--) {
		int j = random() % i;
		T tmp = arr[j];
		arr[j] = arr[i];
		arr[i] = tmp;
	}
}

class rand_permute
{
	off_t *offset;
	long num;

public:
	/**
	 * @start: the index of the first entry.
	 */
	rand_permute(long num, int stride, long start) {
		offset = new off_t[num];
		for (int i = 0; i < num; i++) {
			offset[i] = ((off_t) i) * stride + start * stride;
		}

		rand_permute_array(offset, num);
	}

	~rand_permute() {
		delete [] offset;
	}

	off_t get_offset(long idx) const {
		return offset[idx];
	}
};

long str2size(std::string str);

class sys_parameters
{
	int RAID_block_size;
	int SA_min_cell_size;
	int test_hit_rate;
public:
	sys_parameters() {
		RAID_block_size = get_default_RAID_block_size();
		SA_min_cell_size = get_default_SA_min_cell_size();
	}

	void init(int RAID_block_size, int SA_min_cell_size, int test_hit_rate) {
		this->RAID_block_size = RAID_block_size;
		this->SA_min_cell_size = SA_min_cell_size;
		this->test_hit_rate = test_hit_rate;
	}

	int get_RAID_block_size() {
		return RAID_block_size;
	}

	int get_SA_min_cell_size() {
		return SA_min_cell_size;
	}

	int get_test_hit_rate() {
		return test_hit_rate;
	}

	static int get_default_SA_min_cell_size() {
		return 12;
	}

	static int get_default_RAID_block_size() {
		return 16;
	}
};

extern sys_parameters params;

class io_request;
void extract_pages(const io_request &req, off_t off, int npages,
		io_request &extracted);
bool inside_RAID_block(const io_request &req);

#endif
