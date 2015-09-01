#include <assert.h>

#include "sorter.h"
#include "generic_type.h"

using namespace fm;

const int vec_len = 1000000;

void test_sort()
{
	printf("test sort\n");
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

void test_merge_with_index()
{
	printf("test merge with index\n");
	std::vector<std::vector<long> > vals(10);
	for (size_t i = 0; i < vals.size(); i++) {
		size_t len = random() % vec_len;
		vals[i].resize(len);
#pragma omp parallel for
		for (size_t j = 0; j < vals[i].size(); j++)
			vals[i][j] = random();
		std::sort(vals[i].begin(), vals[i].end());
	}

	type_sorter<long> sort;
	std::vector<std::pair<const char *, const char *> > arrs(vals.size());
	size_t num_eles = 0;
	for (size_t i = 0; i < arrs.size(); i++) {
		num_eles += vals[i].size();
		arrs[i] = std::pair<const char *, const char *>(
				(const char *) vals[i].data(),
				(const char *) (vals[i].data() + vals[i].size()));
	}
	std::vector<std::pair<int, off_t> > merge_index(num_eles);
	std::vector<long> merge_res(num_eles);
	sort.merge_with_index(arrs, (char *) merge_res.data(), num_eles,
			merge_index);

	assert(std::is_sorted(merge_res.begin(), merge_res.end()));

	std::vector<long> merge_res1(num_eles);
	sort.merge(arrs, merge_index, (char *) merge_res1.data(), merge_res1.size());
	for (size_t i = 0; i < merge_res1.size(); i++)
		assert(merge_res[i] == merge_res1[i]);
}

int main()
{
	test_merge_with_index();
	test_sort();
}
