#include <map>

#include "mem_data_frame.h"
#include "bulk_operate.h"
#include "mem_vector.h"
#include "vector_vector.h"
#include "local_vec_store.h"

using namespace fm;

class sum_apply_operate: public gr_apply_operate<sub_data_frame>
{
public:
	void run(const void *key, const sub_data_frame &val,
			local_vec_store &out) const;

	const scalar_type &get_key_type() const {
		return get_scalar_type<int>();
	}

	const scalar_type &get_output_type() const {
		return get_scalar_type<long>();
	}

	size_t get_num_out_eles() const {
		return 1;
	}
};

void sum_apply_operate::run(const void *key, const sub_data_frame &val,
		local_vec_store &out) const
{
	assert(val.size() == 2);
	assert(out.is_type<long>());
	const local_vec_store &vec = *val[1];
	long sum = 0;
	for (size_t i = 0; i < vec.get_length(); i++)
		sum += vec.get<int>(i);
	out.resize(1);
	out.set<long>(0, sum);
}

class copy_apply_operate: public gr_apply_operate<sub_data_frame>
{
public:
	void run(const void *key, const sub_data_frame &val,
			local_vec_store &out) const;

	const scalar_type &get_key_type() const {
		return get_scalar_type<int>();
	}

	const scalar_type &get_output_type() const {
		return get_scalar_type<int>();
	}

	size_t get_num_out_eles() const {
		return 1;
	}
};

void copy_apply_operate::run(const void *key, const sub_data_frame &val,
		local_vec_store &out) const
{
	assert(val.size() == 2);
	assert(out.is_type<int>());
	const local_vec_store &vec = *val[1];
	out.resize(vec.get_length());
	for (size_t i = 0; i < vec.get_length(); i++)
		out.set<int>(i, vec.get<int>(i));
}

void test_groupby()
{
	mem_data_frame::ptr df = mem_data_frame::create();
	size_t length = 1000000;
	detail::mem_vec_store::ptr vec1 = detail::mem_vec_store::create(length,
			get_scalar_type<int>());
	for (size_t i = 0; i < vec1->get_length(); i++)
		vec1->set<int>(i, random() % 1000);
	detail::mem_vec_store::ptr vec2 = detail::mem_vec_store::create(length,
			get_scalar_type<int>());
	for (size_t i = 0; i < vec2->get_length(); i++)
		vec2->set<int>(i, random() % 1000);

	std::map<int, long> map;
	for (size_t i = 0; i < length; i++) {
		int v1 = vec1->get<int>(i);
		int v2 = vec2->get<int>(i);
		auto it = map.find(v1);
		// New value
		if (it == map.end())
			map.insert(std::pair<int, long>(v1, v2));
		else
			it->second += v2;
	}
	df->add_vec("vec1", vec1);
	df->add_vec("vec2", vec2);

	sum_apply_operate sum_op;
	vector_vector::ptr vv = df->groupby("vec1", sum_op);
	mem_vector::ptr v = mem_vector::cast(vv->cat());
	assert(map.size() == v->get_length());
	size_t idx = 0;
	for (auto it = map.begin(); it != map.end(); it++, idx++) {
		assert(it->second == v->get<long>(idx));
	}

	copy_apply_operate copy_op;
	vv = df->groupby("vec1", copy_op);
	mem_vector::ptr v1 = mem_vector::cast(vv->cat());
	assert(vec1->get_length() == v1->get_length());
}

int main()
{
	test_groupby();
}
