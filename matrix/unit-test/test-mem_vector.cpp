#include "mem_vector.h"
#include "bulk_operate.h"
#include "data_frame.h"

using namespace fm;

template<class T>
class find_next_impl: public agg_operate
{
public:
	virtual void run(size_t num_eles, const void *left_arr,
			void *output) const {
		const T *curr = (const T *) left_arr;
		T val = *curr;
		size_t loc = 1;
		for (; loc < num_eles && curr[loc] == val; loc++);
		*(size_t *) output = loc;
	}

	virtual size_t input_entry_size() const {
		return sizeof(T);
	}
	virtual size_t output_entry_size() const {
		return sizeof(size_t);
	}
};

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

	virtual size_t input_entry_size() const {
		return sizeof(T);
	}
	virtual size_t output_entry_size() const {
		return sizeof(size_t);
	}
};

template<class T>
class vec_creator_impl: public vec_creator
{
public:
	virtual size_t get_entry_size() const {
		return sizeof(T);
	}

	virtual vector::ptr create(size_t length) const {
		return type_mem_vector<T>::create(length);
	}
};

void test_groupby()
{
	printf("test groupby\n");
	type_mem_vector<int>::ptr vec = type_mem_vector<int>::create(1000000);
	for (size_t i = 0; i < vec->get_length(); i++)
		vec->set(i, random() % 1000);
	find_next_impl<int> find_next;
	count_impl<int> count;
	vec_creator_impl<size_t> create;
	data_frame::ptr res = vec->groupby(find_next, count, create);
	printf("size: %ld\n", res->get_num_entries());

	std::map<int, size_t> ele_counts;
	for (size_t i = 0; i < vec->get_length(); i++) {
		int val = vec->get(i);
		auto it = ele_counts.find(val);
		if (it == ele_counts.end())
			ele_counts.insert(std::pair<int, size_t>(val, 1));
		else
			it->second++;
	}

	type_mem_vector<int>::ptr vals = type_mem_vector<int>::cast(
			res->get_vec("val"));
	type_mem_vector<size_t>::ptr aggs = type_mem_vector<size_t>::cast(
			res->get_vec("agg"));
	for (size_t i = 0; i < vals->get_length(); i++) {
		int val = vals->get(i);
		size_t count = aggs->get(i);
		auto it = ele_counts.find(val);
		assert(it != ele_counts.end());
		assert(it->second == count);
	}
}

void test_append()
{
	printf("test append\n");
	type_mem_vector<int>::ptr res = type_mem_vector<int>::create(16);
	size_t tot_len = res->get_length();
	std::vector<vector::ptr> vecs(16);
	for (size_t i = 0; i < vecs.size(); i++) {
		vecs[i] = std::static_pointer_cast<vector>(
				type_mem_vector<int>::create(32));
		tot_len += vecs[i]->get_length();
	}
	res->append(vecs.begin(), vecs.end());
	assert(tot_len == res->get_length());
}

int main()
{
	test_append();
	test_groupby();
}
