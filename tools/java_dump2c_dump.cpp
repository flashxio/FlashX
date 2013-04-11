#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <vector>

#include "workload.h"

int main(int argc, char *argv[])
{
	char *workload_file = argv[1];
	char *output_file = argv[2];
	long length = 0;
	off_t *workloads = NULL;
	workloads = load_java_dump(workload_file, length);
	java_dump_workload *gen = new java_dump_workload(workloads, length, 0, length);
	workload_gen::set_default_access_method(READ);

	std::vector<workload_t> cworkloads;
	while(gen->has_next()) {
		workload_t w = gen->next();
		assert(w.size == PAGE_SIZE);
		cworkloads.push_back(w);
	}

	int fd = open(output_file, O_WRONLY | O_CREAT, S_IRWXU);
	if (fd < 0) {
		perror("open");
		return -1;
	}

	workload_t *w = cworkloads.data();
	for (size_t i = 0; i < cworkloads.size(); i++) {
		assert(w[i].size == PAGE_SIZE);
	}
	char *buf = (char *) w;
	size_t size = cworkloads.size() * sizeof(workload_t);
	ssize_t ret;
	do {
		ret = write(fd, buf, size);
		if (ret < 0) {
			perror("write");
			return -1;
		}
		buf += ret;
		size -= ret;
	} while (size > 0);
}
