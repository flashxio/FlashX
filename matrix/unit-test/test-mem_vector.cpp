#include "mem_vector.h"
#include "bulk_operate.h"
#include "data_frame.h"

using namespace fm;

template<class T>
class count_impl: public agg_operate
{
public:
	virtual void run(size_t num_eles, const void *left_arr,
			void *output) const {
		assert(num_eles > 0);
		const T *arr = (const T *) left_arr;
		T val = *arr;
		if (num_eles > 1)
			for (size_t i = 1; i < num_eles; i++)
				assert(val == arr[i]);
		*(size_t *) output = num_eles;
	}

	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<T>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<size_t>();
	}
};

void test_groupby()
{
	printf("test groupby\n");
	mem_vector::ptr vec = mem_vector::create(1000000, get_scalar_type<int>());
	for (size_t i = 0; i < vec->get_length(); i++)
		vec->set<int>(i, random() % 1000);
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

	mem_vector::ptr vals = mem_vector::cast(res->get_vec("val"));
	mem_vector::ptr aggs = mem_vector::cast(res->get_vec("agg"));
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
	mem_vector::ptr res = mem_vector::create(16, get_scalar_type<int>());
	size_t tot_len = res->get_length();
	std::vector<vector::ptr> vecs(16);
	for (size_t i = 0; i < vecs.size(); i++) {
		vecs[i] = std::static_pointer_cast<vector>(
				mem_vector::create(32, get_scalar_type<int>()));
		tot_len += vecs[i]->get_length();
	}
	res->append(vecs.begin(), vecs.end());
	assert(tot_len == res->get_length());
}

void test_sort()
{
	printf("test sort\n");
	mem_vector::ptr vec = mem_vector::create(1000000, get_scalar_type<int>());
	for (size_t i = 0; i < vec->get_length(); i++)
		vec->set(i, random() % 1000);
	mem_vector::ptr clone = mem_vector::cast(vec->deep_copy());
	assert(clone->equals(*vec));

	mem_vector::ptr idxs = mem_vector::cast(vec->sort_with_index());
	mem_vector::ptr sorted = clone->mem_vector::get(*idxs);
	assert(sorted->equals(*vec));
}

int main()
{
	test_sort();
	test_append();
	test_groupby();
}
