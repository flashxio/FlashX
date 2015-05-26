#include <assert.h>

#include "sorter.h"
#include "generic_type.h"

using namespace fm;

const int vec_len = 1000000;

void test_sort()
{
	std::vector<std::vector<long> > vals(50);
	for (size_t i = 0; i < vals.size(); i++) {
		vals[i].resize(vec_len + random() % vec_len);
		for (size_t j = 0; j < vals[i].size(); j++)
			vals[i][j] = random();
		get_scalar_type<long>().get_sorter().sort((char *) vals[i].data(),
				vals[i].size(), false);
		assert(std::is_sorted(vals[i].begin(), vals[i].end()));
	}
	std::vector<std::pair<const char *, const char *> > pairs(vals.size());
	size_t tot_vals = 0;
	for (size_t i = 0; i < vals.size(); i++) {
		tot_vals += vals[i].size();
		pairs[i] = std::pair<const char *, const char *>((char *) vals[i].data(),
				(char *) (vals[i].data() + vals[i].size()));
	}
	std::vector<long> res(tot_vals);
	get_scalar_type<long>().get_sorter().merge(pairs, (char *) res.data(),
			tot_vals);
	assert(std::is_sorted(res.begin(), res.end()));
}

int main()
{
	test_sort();
}
