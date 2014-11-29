#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "config_map.h"
#include "common.h"

#include "safs_file.h"
#include "native_file.h"
#include "io_interface.h"

#include "utils.h"

using namespace safs;
using namespace fg;

namespace fg
{
namespace utils
{
	void set_buf_cap(size_t new_cap);
}
}

size_t test_write(utils::large_writer::ptr writer)
{
	size_t tot_size = 0;
	int initializer = 0;
	size_t tot_write_bytes = 0;
	for (int i = 0; i < 100; i++) {
		size_t ints = random() % (1024 * 256);
		int *buf = (int *) malloc(ints * sizeof(int));
		for (int i = 0; i < ints; i++)
			buf[i] = initializer++;
		ssize_t ret = writer->write((char *) buf, ints * sizeof(int));
		assert(ret >= 0);
		tot_write_bytes += ret;
		free(buf);
		tot_size += ints * sizeof(int);
	}
	assert(tot_size == tot_write_bytes);
	printf("attemp to write %ld bytes\n", tot_size);
	return tot_size;
}

void test_read(utils::large_reader::ptr reader, size_t tot_size)
{
	int *buf = (int *) malloc(tot_size);
	ssize_t ret = reader->read((char *) buf, tot_size);
	assert(ret == tot_size);
	for (size_t i = 0; i < tot_size / sizeof(buf[0]); i++) {
		if (buf[i] != i)
			printf("buf[%ld]: %d, i: %ld\n", i, buf[i], i);
		assert(buf[i] == i);
	}
	free(buf);
}

int main(int argc, char *argv[])
{
	utils::set_buf_cap(1024 * 128);
	utils::large_io_creator::ptr creator
		= utils::large_io_creator::create(false, ".");

	char *tmp_file = tempnam(".", "large_writer");
	utils::large_writer::ptr writer = creator->create_writer(tmp_file);
	size_t tot_size = test_write(writer);
	writer.reset();
	native_file f(tmp_file);
	printf("%s has %ld bytes\n", tmp_file, f.get_size());
	assert(f.get_size() == ROUNDUP(tot_size, 512));
	utils::large_reader::ptr reader = creator->create_reader(tmp_file);
	test_read(reader, tot_size);
	unlink(tmp_file);

	if (argc < 2) {
		fprintf(stderr, "we need a SAFS conf file to test large_io for SAFS\n");
		exit(1);
	}

	config_map::ptr configs = config_map::create(argv[1]);
	configs->add_options("writable=1");
	init_io_system(configs);
	assert(is_safs_init());

	creator = utils::large_io_creator::create(true, ".");
	writer = creator->create_writer(tmp_file);
	assert(writer != NULL);
	tot_size = test_write(writer);
	writer.reset();
	safs_file sf(get_sys_RAID_conf(), tmp_file);
	printf("write %ld bytes and %s has %ld bytes\n", tot_size, tmp_file,
			sf.get_size());
	assert(sf.get_size() >= ROUNDUP(tot_size, 512));
	reader = creator->create_reader(tmp_file);
	test_read(reader, tot_size);
	sf.delete_file();
}
