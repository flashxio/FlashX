#include "mem_data_frame.h"
#include "bulk_operate.h"
#include "mem_vector.h"
#include "vector_vector.h"

using namespace fm;

class sum_apply_operate: public gr_apply_operate<data_frame>
{
public:
	void run(const void *key, const data_frame &val, mem_vector &out) const;

	const scalar_type &get_key_type() const {
		static scalar_type_impl<int> t;
		return t;
	}

	const scalar_type &get_output_type() const {
		static scalar_type_impl<long> t;
		return t;
	}

	size_t get_num_out_eles() const {
		return 1;
	}
};

void sum_apply_operate::run(const void *key, const data_frame &val,
		mem_vector &out) const
{
	assert(val.get_num_vecs() == 2);
	assert(out.is_type<long>());
	type_mem_vector<int>::ptr vec = type_mem_vector<int>::cast(val.get_vec(1));
	long sum = 0;
	for (size_t i = 0; i < vec->get_length(); i++)
		sum += vec->get(i);
	type_mem_vector<long> &t_out = (type_mem_vector<long> &) out;
	out.resize(1);
	t_out.set(0, sum);
}

void test_groupby()
{
	mem_data_frame::ptr df = mem_data_frame::create();
	size_t length = 1000000;
	type_mem_vector<int>::ptr vec1 = type_mem_vector<int>::create(length);
	for (size_t i = 0; i < vec1->get_length(); i++)
		vec1->set(i, random() % 1000);
	type_mem_vector<int>::ptr vec2 = type_mem_vector<int>::create(length);
	for (size_t i = 0; i < vec2->get_length(); i++)
		vec2->set(i, random() % 1000);

	std::map<int, long> map;
	for (size_t i = 0; i < length; i++) {
		int v1 = vec1->get(i);
		int v2 = vec2->get(i);
		auto it = map.find(v1);
		// New value
		if (it == map.end())
			map.insert(std::pair<int, long>(v1, v2));
		else
			it->second += v2;
	}

	df->add_vec("vec1", vec1);
	df->add_vec("vec2", vec2);
	sum_apply_operate op;
	vector_vector::ptr vv = df->groupby("vec1", op);
	type_mem_vector<long>::ptr v = type_mem_vector<long>::cast(vv->cat());
	assert(map.size() == v->get_length());
	size_t idx = 0;
	for (auto it = map.begin(); it != map.end(); it++, idx++) {
		assert(it->second == v->get(idx));
	}
}

int main()
{
	test_groupby();
}
