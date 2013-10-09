#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <sched.h>
#include <assert.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <signal.h>

#include "common.h"
#include "messaging.h"

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

static void enable_debug_handler(int sig, siginfo_t *si, void *uc)
{
	enable_debug = true;
	printf("debug mode is enabled\n");
}

}

bool check_read_content(char *buf, int size, off_t off)
{
	// I assume the space in the buffer is larger than 8 bytes.
	off_t aligned_off = off & (~(sizeof(off_t) - 1));
	long data[2];
	data[0] = aligned_off / sizeof(off_t);
	data[1] = aligned_off / sizeof(off_t) + 1;
	long expected = 0;
	int copy_size = size < (int) sizeof(off_t) ? size : (int) sizeof(off_t);
	memcpy(&expected, ((char *) data) + (off - aligned_off), copy_size);
	long read_value = 0;
	memcpy(&read_value, buf, copy_size);
	if(read_value != expected)
		printf("off: %ld, size: %d, read: %ld, expect: %ld\n",
				off, size, read_value, expected);
	return read_value == expected;
}

void create_write_data(char *buf, int size, off_t off)
{
	off_t aligned_start = off & (~(sizeof(off_t) - 1));
	off_t aligned_end = (off + size) & (~(sizeof(off_t) - 1));
	long start_data = aligned_start / sizeof(off_t);
	long end_data = aligned_end / sizeof(off_t);

	/* If all data is in one 8-byte word. */
	if (aligned_start == aligned_end) {
		memcpy(buf, ((char *) &start_data) + (off - aligned_start), size);
		return;
	}

	int first_size =  (int)(sizeof(off_t) - (off - aligned_start));
	int last_size = (int) (off + size - aligned_end);

	if (first_size == sizeof(off_t))
		first_size = 0;
	if (first_size)
		memcpy(buf, ((char *) &start_data) + (off - aligned_start),
				first_size);
	for (int i = first_size; i < aligned_end - off; i += sizeof(off_t)) {
		*((long *) (buf + i)) = (off + i) / sizeof(off_t);
	}
	if (aligned_end > aligned_start
			|| (aligned_end == aligned_start && first_size == 0)) {
		if (last_size)
			memcpy(buf + (aligned_end - off), (char *) &end_data, last_size);
	}

	check_read_content(buf, size, off);
}

bool enable_debug = false;

void set_enable_debug_signal()
{
	struct sigaction sa;

	/* Establish handler for timer signal */

	sa.sa_flags = SA_SIGINFO;
	sa.sa_sigaction = enable_debug_handler;
	sigemptyset(&sa.sa_mask);
	if (sigaction(SIGUSR1, &sa, NULL) == -1) {
		perror("sigaction");
		exit(1);
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
		assert(0);
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

/**
 * Extract a request from the input request.
 * The extract request is within the range [off, off + npages * PAGE_SIZE),
 * where off is aligned with PAGE_SIZE.
 */
void extract_pages(const io_request &req, off_t off, int npages,
		io_request &extracted)
{
	off_t req_off;
	char *req_buf;
	ssize_t req_size;
	assert(req.get_num_bufs() == 1);
	assert((off & (PAGE_SIZE - 1)) == 0);
	bool check = (off >= req.get_offset() && off < req.get_offset() + req.get_size())
		|| (off + PAGE_SIZE >= req.get_offset()
				&& off + PAGE_SIZE < req.get_offset() + req.get_size())
		|| (off <= req.get_offset()
				&& off + PAGE_SIZE >= req.get_offset() + req.get_size());
	if (!check)
		fprintf(stderr, "req %lx, size: %lx, page off: %lx\n",
				req.get_offset(), req.get_size(), off);
	assert(check);
	// this is the first page in the request.
	if (off == ROUND_PAGE(req.get_offset())) {
		req_off = req.get_offset();
		req_buf = req.get_buf();
		// the remaining size in the page.
		req_size = PAGE_SIZE * npages - (req_off - off);
		if (req_size > req.get_size())
			req_size = req.get_size();
	}
	else {
		req_off = off;
		/* 
		 * We can't be sure if the request buffer is aligned
		 * with the page size.
		 */
		req_buf = req.get_buf() + (off - req.get_offset());
		ssize_t remaining = req.get_size() - (off - req.get_offset());
		req_size = remaining > PAGE_SIZE * npages ? PAGE_SIZE
			* npages : remaining;
	}
	extracted.init(req_buf, req_off, req_size, req.get_access_method(),
			req.get_io(), req.get_node_id());
}

bool inside_RAID_block(const io_request &req)
{
	int RAID_block_size = params.get_RAID_block_size() * PAGE_SIZE;
	return ROUND(req.get_offset(), RAID_block_size)
		== ROUND(req.get_offset() + req.get_size() - 1, RAID_block_size);
}
