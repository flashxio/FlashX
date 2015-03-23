#include "mem_vector.h"

using namespace fm;

mem_vector::ptr get(const mem_vector &vec, mem_vector &idxs)
{
	mem_vector::ptr ret = mem_vector::create(idxs.get_length(), vec.get_type());
#pragma omp parallel for
	for (size_t i = 0; i < idxs.get_length(); i++) {
		off_t idx = idxs.get<off_t>(i);
		// Check if it's out of the range.
		if (idx < 0 && (size_t) idx >= vec.get_length()) {
			BOOST_LOG_TRIVIAL(error)
				<< boost::format("%1% is out of range") % idx;
			continue;
		}

		ret->set(i, vec.get<int>(idx));
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
	mem_vector::ptr vec = mem_vector::create(1000000000, get_scalar_type<int>());
	for (size_t i = 0; i < vec->get_length(); i++)
		vec->set(i, random());
	mem_vector::ptr clone = mem_vector::cast(vec->deep_copy());
	assert(clone->equals(*vec));

	struct timeval start, end;
	gettimeofday(&start, NULL);
	mem_vector::ptr idxs = mem_vector::cast(vec->sort_with_index());
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
