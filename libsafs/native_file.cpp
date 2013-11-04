#include <stddef.h>
#include <assert.h>

#include "native_file.h"

ssize_t native_dir::read_all_files(std::vector<std::string> &files) const
{
	DIR *d = opendir(name.c_str());
	if (d == NULL) {
		fprintf(stderr, "can't open %s: %s\n", name.c_str(),
				strerror(errno));
		return -1;
	}
	int len = offsetof(struct dirent, d_name) +
		pathconf("dirpath", _PC_NAME_MAX) + 1;
	struct dirent *entryp = (struct dirent *) malloc(len);
	struct dirent *result;
	ssize_t num = 0;
	while (true) {
		int ret = readdir_r(d, entryp, &result);
		if (ret != 0) {
			fprintf(stderr, "can't read dir %s: %s\n", name.c_str(),
					strerror(ret));
			num = -1;
			break;
		}
		// End of the directory stream.
		if (result == NULL)
			break;

		assert(result == entryp);
		files.push_back(entryp->d_name);
		num++;
	}
	free(entryp);
	return num;
}
