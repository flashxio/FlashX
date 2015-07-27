#include <map>

#include "vector.h"
#include "bulk_operate.h"
#include "data_frame.h"
#include "local_vec_store.h"
#include "dense_matrix.h"

using namespace fm;
using namespace detail;

template<class T>
class count_impl: public gr_apply_operate<local_vec_store>
{
public:
	virtual void run(const void *key, const local_vec_store &vec,
			local_vec_store &output) const {
		size_t num_eles = vec.get_length();
		assert(num_eles > 0);
		const T *arr = (const T *) vec.get_raw_arr();
		T val = *arr;
		if (num_eles > 1)
			for (size_t i = 1; i < num_eles; i++)
				assert(val == arr[i]);
		assert(output.get_type() == get_scalar_type<size_t>());
		output.set<size_t>(0, num_eles);
	}

	virtual const scalar_type &get_key_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<size_t>();
	}
	virtual size_t get_num_out_eles() const {
		return 1;
	}
};

void test_groupby()
{
	printf("test groupby\n");
	smp_vec_store::ptr store = smp_vec_store::create(1000000,
			get_scalar_type<int>());
	for (size_t i = 0; i < store->get_length(); i++)
		store->set<int>(i, random() % 1000);
	vector::ptr vec = vector::create(store);
	count_impl<int> count;
	data_frame::ptr res = vec->groupby(count, true);
	printf("size: %ld\n", res->get_num_entries());

	std::map<int, size_t> ele_counts;
	const detail::smp_vec_store &vstore
		= dynamic_cast<const detail::smp_vec_store &>(vec->get_data());
	for (size_t i = 0; i < vec->get_length(); i++) {
		int val = vstore.get<int>(i);
		auto it = ele_counts.find(val);
		if (it == ele_counts.end())
			ele_counts.insert(std::pair<int, size_t>(val, 1));
		else
			it->second++;
	}

	smp_vec_store::ptr vals = smp_vec_store::cast(res->get_vec("val"));
	smp_vec_store::ptr aggs = smp_vec_store::cast(res->get_vec("agg"));
	for (size_t i = 0; i < vals->get_length(); i++) {
		int val = vals->get<int>(i);
		size_t count = aggs->get<size_t>(i);
		auto it = ele_counts.find(val);
		assert(it != ele_counts.end());
		assert(it->second == count);
	}
}

void test_append()
{
	printf("test append\n");
	size_t tot_len = 0;
	std::vector<vec_store::const_ptr> vecs(16);
	for (size_t i = 0; i < vecs.size(); i++) {
		vecs[i] = smp_vec_store::create(32, get_scalar_type<int>());
		tot_len += vecs[i]->get_length();
	}

	smp_vec_store::ptr res1 = smp_vec_store::create(0, get_scalar_type<int>());
	for (size_t i = 0; i < vecs.size(); i++) {
		off_t off = res1->get_length();
		res1->append(*vecs[i]);

		const char *sub_vec = res1->get(off);
		smp_vec_store::const_ptr orig = smp_vec_store::cast(vecs[i]);
		assert(memcmp(sub_vec, orig->get_raw_arr(),
					orig->get_length() * orig->get_entry_size()) == 0);
	}

	smp_vec_store::ptr res2 = smp_vec_store::create(0, get_scalar_type<int>());
	res2->append(vecs.begin(), vecs.end());
	assert(tot_len == res2->get_length());
	off_t off2 = 0;
	for (size_t i = 0; i < vecs.size(); i++) {
		const char *sub_vec = res2->get(off2);
		smp_vec_store::const_ptr orig = smp_vec_store::cast(vecs[i]);
		assert(memcmp(sub_vec, orig->get_raw_arr(),
					orig->get_length() * res2->get_entry_size()) == 0);
		off2 += orig->get_length();
	}
}

void test_sort()
{
	printf("test sort\n");
	smp_vec_store::ptr vec = smp_vec_store::create(1000000, get_scalar_type<int>());
	for (size_t i = 0; i < vec->get_length(); i++)
		vec->set<int>(i, random() % 1000);
	smp_vec_store::ptr clone = smp_vec_store::cast(vec->deep_copy());
	vector::ptr vec1 = vector::create(clone);
	vector::ptr vec2 = vector::create(vec);
	assert(vec1->equals(*vec2));

	smp_vec_store::ptr idxs = smp_vec_store::cast(vec->sort_with_index());
	smp_vec_store::ptr sorted = clone->smp_vec_store::get(*idxs);
	vec1 = vector::create(sorted);
	vec2 = vector::create(vec);
	assert(vec1->equals(*vec2));
}

void test_max()
{
	printf("test max\n");
	smp_vec_store::ptr vec = smp_vec_store::create(1000000, get_scalar_type<int>());
	int max = 0;
	for (size_t i = 0; i < vec->get_length(); i++) {
		int v = random() % 1000;
		vec->set<int>(i, v);
		max = std::max(max, v);
	}
	vector::ptr vec1 = vector::create(vec);
	assert(vec1->max<int>() == max);
}

void test_resize()
{
	printf("test resize\n");
	vector::ptr vec = create_vector<int>(1, 10000, 2);
	const detail::smp_vec_store &vstore
		= dynamic_cast<const detail::smp_vec_store &>(vec->get_data());
	smp_vec_store::ptr copy = smp_vec_store::cast(vec->get_data().deep_copy());
	copy->resize(100);
	size_t min_len = std::min(copy->get_length(), vec->get_length());
	for (size_t i = 0; i < min_len; i++)
		assert(vstore.get<int>(i) == copy->get<int>(i));

	copy->resize(200);
	// The semantics don't guarantee that this works, but it works with
	// the current implementation
	min_len = std::min(copy->get_length(), vec->get_length());
	for (size_t i = 0; i < min_len; i++)
		assert(vstore.get<int>(i) == copy->get<int>(i));

	copy->resize(20000);
	for (size_t i = 0; i < min_len; i++)
		assert(vstore.get<int>(i) == copy->get<int>(i));
}

void test_get_sub()
{
	printf("test get sub\n");
	vector::ptr vec = create_vector<int>(1, 10000, 2);
	smp_vec_store::ptr idxs = smp_vec_store::create(1000, get_scalar_type<off_t>());
	for (size_t i = 0; i < idxs->get_length(); i++)
		idxs->set<off_t>(i, random() % idxs->get_length());
	const detail::smp_vec_store &vstore
		= dynamic_cast<const detail::smp_vec_store &>(vec->get_data());
	smp_vec_store::ptr res = vstore.get(*idxs);
	for (size_t i = 0; i < res->get_length(); i++)
		assert(res->get<int>(i) == vstore.get<int>(idxs->get<off_t>(i)));
}

void test_copy_from()
{
	printf("test copy vector\n");
	smp_vec_store::ptr vec = smp_vec_store::create(1000000,
			get_scalar_type<int>());
	std::vector<int> stl_vec(vec->get_length());
	for (size_t i = 0; i < stl_vec.size(); i++)
		stl_vec[i] = i;
	vec->copy_from((const char *) stl_vec.data(),
			vec->get_length() * vec->get_entry_size());
	for (size_t i = 0; i < stl_vec.size(); i++)
		assert(stl_vec[i] == vec->get<int>(i));
}

void test_conv2std()
{
	printf("test convert to a std vector\n");
	vector::ptr vec = create_vector<int>(1, 10000, 2);
	std::vector<int> stl_vec = vec->conv2std<int>();
	assert(stl_vec.size() == vec->get_length());
	for (size_t i = 0; i < stl_vec.size(); i++)
		assert(stl_vec[i] == 1 + 2 * i);
}

int main()
{
	test_sort();
	test_append();
	test_groupby();
	test_max();
	test_resize();
	test_get_sub();
	test_copy_from();
	test_conv2std();
}
