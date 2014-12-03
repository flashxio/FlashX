#include <malloc.h>

#include "EM_vector.h"
#include "sparse_matrix.h"

using namespace fm;

class write_double_compute: public subvec_compute
{
	double *buf;
public:
	write_double_compute(double *buf) {
		this->buf = buf;
	}

	~write_double_compute() {
		free(buf);
	}

	virtual void run(char *buf, size_t size) {
	}
};

class read_double_compute: public subvec_compute
{
	off_t idx;
	size_t num;
public:
	read_double_compute(off_t idx, size_t num) {
		this->idx = idx;
		this->num = num;
	}

	virtual void run(char *buf, size_t size) {
		double *dbuf = (double *) buf;
		assert(sizeof(double) * num == size);
		for (size_t i = 0; i < num; i++)
			assert(dbuf[i] == idx + i);
	}
};

void init_buf(double *buf, off_t idx, size_t num)
{
	for (size_t i = 0; i < num; i++)
		buf[i] = idx + i;
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "test conf_file\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	EM_vector::ptr vec = EM_vector::create(1024 * 1024, sizeof(double));
	EM_vector_accessor::ptr accessor = vec->create_accessor();
	size_t num_eles = 1024 * 8;
	for (size_t i = 0; i < vec->get_size(); i += num_eles) {
		double *buf = (double *) memalign(PAGE_SIZE, num_eles * sizeof(double));
		init_buf(buf, i, num_eles);
		accessor->set_subvec((char *) buf, i, num_eles, subvec_compute::ptr(
				new write_double_compute(buf)));
	}
	accessor->wait4all();

	for (size_t i = 0; i < vec->get_size(); i += num_eles) {
		accessor->fetch_subvec(i, num_eles, subvec_compute::ptr(
				new read_double_compute(i, num_eles)));
	}
	accessor->wait4all();

	destroy_flash_matrix();
}
