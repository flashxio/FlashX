#include <stdio.h>
#include <assert.h>

#include "native_file.h"

using namespace safs;

int main()
{
	std::string dir_name = "/tmp/1/1/1";
	native_dir dir(dir_name);
	assert(!dir.exist());
	assert(dir.create_dir(true));

	native_dir dir2("/tmp/1");
	assert(dir2.delete_dir(true));

	dir_name = "1/1/1";
	native_dir dir3(dir_name);
	assert(!dir3.exist());
	assert(dir3.create_dir(true));

	native_dir dir4("1");
	assert(dir4.delete_dir(true));
}
