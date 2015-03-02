#include "mem_vector.h"

using namespace fm;

mem_vector::ptr get(const type_mem_vector<int> &vec, type_mem_vector<off_t> &idxs)
{
	type_mem_vector<int>::ptr ret
		= type_mem_vector<int>::create(idxs.get_length());
#pragma omp parallel for
	for (size_t i = 0; i < idxs.get_length(); i++) {
		off_t idx = idxs.get(i);
		// Check if it's out of the range.
		if (idx < 0 && (size_t) idx >= vec.get_length()) {
			BOOST_LOG_TRIVIAL(error)
				<< boost::format("%1% is out of range") % idx;
			continue;
		}

		ret->set(i, vec.get(idx));
	}
	return std::static_pointer_cast<mem_vector>(ret);
}

/*
 * This is to measure the performance difference of random permutation
 * with and without compile-time type.
 */
void test_permute()
{
	printf("test permutation\n");
	type_mem_vector<int>::ptr vec = type_mem_vector<int>::create(1000000000);
	for (size_t i = 0; i < vec->get_length(); i++)
		vec->set(i, random());
	type_mem_vector<int>::ptr clone = type_mem_vector<int>::create(1000000000);
	clone->set_sub_vec(0, *vec);
	assert(clone->equals(*vec));

	struct timeval start, end;
	gettimeofday(&start, NULL);
	type_mem_vector<off_t>::ptr idxs = type_mem_vector<off_t>::cast(
			vec->sort_with_index());
	gettimeofday(&end, NULL);
	printf("sort takes %fseconds\n", time_diff(start, end));

	gettimeofday(&start, NULL);
	// This has compile-time type.
	mem_vector::ptr sorted1 = get(*clone, *idxs);
	gettimeofday(&end, NULL);
	printf("permute with type takes %fseconds\n", time_diff(start, end));

	gettimeofday(&start, NULL);
	// This doesn't have compile-time type.
	mem_vector::ptr sorted2 = clone->mem_vector::get(*idxs);
	gettimeofday(&end, NULL);
	printf("permute without type takes %fseconds\n", time_diff(start, end));
	assert(sorted1->equals(*sorted2));
}

int main()
{
	test_permute();
}
