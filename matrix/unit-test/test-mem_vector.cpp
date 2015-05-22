#include <map>

#include "mem_vector.h"
#include "bulk_operate.h"
#include "data_frame.h"
#include "local_vec_store.h"

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
	mem_vec_store::ptr store = mem_vec_store::create(1000000,
			get_scalar_type<int>());
	for (size_t i = 0; i < store->get_length(); i++)
		store->set<int>(i, random() % 1000);
	mem_vector::ptr vec = mem_vector::create(store);
	count_impl<int> count;
	data_frame::ptr res = vec->groupby(count, true);
	printf("size: %ld\n", res->get_num_entries());

	std::map<int, size_t> ele_counts;
	for (size_t i = 0; i < vec->get_length(); i++) {
		int val = vec->get<int>(i);
		auto it = ele_counts.find(val);
		if (it == ele_counts.end())
			ele_counts.insert(std::pair<int, size_t>(val, 1));
		else
			it->second++;
	}

	mem_vec_store::ptr vals = mem_vec_store::cast(res->get_vec("val"));
	mem_vec_store::ptr aggs = mem_vec_store::cast(res->get_vec("agg"));
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
	mem_vec_store::ptr res = mem_vec_store::create(16, get_scalar_type<int>());
	size_t tot_len = res->get_length();
	std::vector<vec_store::const_ptr> vecs(16);
	for (size_t i = 0; i < vecs.size(); i++) {
		vecs[i] = mem_vec_store::create(32, get_scalar_type<int>());
		tot_len += vecs[i]->get_length();
	}
	res->append(vecs.begin(), vecs.end());
	assert(tot_len == res->get_length());
}

void test_sort()
{
	printf("test sort\n");
	mem_vec_store::ptr vec = mem_vec_store::create(1000000, get_scalar_type<int>());
	for (size_t i = 0; i < vec->get_length(); i++)
		vec->set<int>(i, random() % 1000);
	mem_vec_store::ptr clone = mem_vec_store::cast(vec->deep_copy());
	mem_vector::ptr vec1 = mem_vector::create(clone);
	mem_vector::ptr vec2 = mem_vector::create(vec);
	assert(vec1->equals(*vec2));

	mem_vec_store::ptr idxs = mem_vec_store::cast(vec->sort_with_index());
	mem_vec_store::ptr sorted = clone->mem_vec_store::get(*idxs);
	vec1 = mem_vector::create(sorted);
	vec2 = mem_vector::create(vec);
	assert(vec1->equals(*vec2));
}

void test_max()
{
	printf("test max\n");
	mem_vec_store::ptr vec = mem_vec_store::create(1000000, get_scalar_type<int>());
	int max = 0;
	for (size_t i = 0; i < vec->get_length(); i++) {
		int v = random() % 1000;
		vec->set<int>(i, v);
		max = std::max(max, v);
	}
	mem_vector::ptr vec1 = mem_vector::create(vec);
	assert(vec1->max<int>() == max);
}

void test_resize()
{
	printf("test resize\n");
	mem_vector::ptr vec = mem_vector::cast(create_vector<int>(1, 10000, 2));
	mem_vec_store::ptr copy = mem_vec_store::cast(vec->get_data().deep_copy());
	copy->resize(100);
	size_t min_len = std::min(copy->get_length(), vec->get_length());
	for (size_t i = 0; i < min_len; i++)
		assert(vec->get<int>(i) == copy->get<int>(i));

	copy->resize(200);
	// The semantics don't guarantee that this works, but it works with
	// the current implementation
	min_len = std::min(copy->get_length(), vec->get_length());
	for (size_t i = 0; i < min_len; i++)
		assert(vec->get<int>(i) == copy->get<int>(i));

	copy->resize(20000);
	for (size_t i = 0; i < min_len; i++)
		assert(vec->get<int>(i) == copy->get<int>(i));
}

void test_get_sub()
{
	printf("test get sub\n");
	mem_vector::ptr vec = mem_vector::cast(create_vector<int>(1, 10000, 2));
	mem_vec_store::ptr idxs = mem_vec_store::create(1000, get_scalar_type<off_t>());
	for (size_t i = 0; i < idxs->get_length(); i++)
		idxs->set<off_t>(i, random() % idxs->get_length());
	mem_vec_store::ptr res = mem_vec_store::cast(vec->get_raw_store())->get(*idxs);
	for (size_t i = 0; i < res->get_length(); i++)
		assert(res->get<int>(i) == vec->get<int>(idxs->get<off_t>(i)));
}

int main()
{
	test_sort();
	test_append();
	test_groupby();
	test_max();
	test_resize();
	test_get_sub();
}
