#include <boost/format.hpp>

#include "log.h"

#include "mem_vec_store.h"
#include "mem_vector.h"

using namespace fm;

detail::mem_vec_store::ptr get(const detail::mem_vec_store &vec,
		detail::mem_vec_store &idxs)
{
	detail::mem_vec_store::ptr ret = detail::mem_vec_store::create(idxs.get_length(),
			vec.get_type());
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
	return ret;
}

/*
 * This is to measure the performance difference of random permutation
 * with and without compile-time type.
 */
void test_permute()
{
	printf("test permutation\n");
	detail::mem_vec_store::ptr vec = detail::mem_vec_store::create(1000000000,
			get_scalar_type<int>());
	for (size_t i = 0; i < vec->get_length(); i++)
		vec->set(i, random());
	detail::mem_vec_store::ptr clone = detail::mem_vec_store::cast(vec->deep_copy());
	mem_vector::ptr vec1 = mem_vector::create(clone);
	mem_vector::ptr vec2 = mem_vector::create(vec);
	assert(vec1->equals(*vec2));

	struct timeval start, end;
	gettimeofday(&start, NULL);
	detail::mem_vec_store::ptr idxs = detail::mem_vec_store::cast(
			vec->sort_with_index());
	gettimeofday(&end, NULL);
	printf("sort takes %fseconds\n", time_diff(start, end));

	gettimeofday(&start, NULL);
	// This has compile-time type.
	detail::mem_vec_store::ptr sorted1 = get(*clone, *idxs);
	gettimeofday(&end, NULL);
	printf("permute with type takes %fseconds\n", time_diff(start, end));

	gettimeofday(&start, NULL);
	// This doesn't have compile-time type.
	detail::mem_vec_store::ptr sorted2 = clone->detail::mem_vec_store::get(*idxs);
	gettimeofday(&end, NULL);
	printf("permute without type takes %fseconds\n", time_diff(start, end));
	vec1 = mem_vector::create(sorted1);
	vec2 = mem_vector::create(sorted2);
	assert(vec1->equals(*vec2));
}

int main()
{
	test_permute();
}
