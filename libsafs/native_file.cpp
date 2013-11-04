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
	struct dirent *entryp;
	long path_limit = pathconf(name.c_str(), _PC_NAME_MAX);
	if (path_limit < 0) {
		perror("pathconf");
		path_limit = sizeof(entryp->d_name);
	}
	int len = offsetof(struct dirent, d_name) + path_limit + 1;
	entryp = (struct dirent *) malloc(len);
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
		std::string file = entryp->d_name;
		if (file != "." && file != "..")
			files.push_back(file);
		num++;
	}
	free(entryp);
	return num;
}
