#include "mem_vector.h"
#include "bulk_operate.h"
#include "data_frame.h"
#include "factor.h"
#include "vector_vector.h"

using namespace fm;

class adj_apply: public gr_apply_operate<mem_vector>
{
public:
	virtual void run(const void *key, const mem_vector &val,
			mem_vector &vec) const {
		vec.resize(val.get_length());
		vec.set_sub_vec(0, val);
	}

	virtual const scalar_type &get_key_type() const {
		return get_scalar_type<int>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<int>();
	}
	virtual size_t get_num_out_eles() const {
		return 0;
	}
};

class set_label_operate: public type_set_operate<factor_value_t>
{
	factor f;
public:
	set_label_operate(const factor &_f): f(_f) {
	}

	virtual void set(factor_value_t *arr, size_t num_eles, off_t row_idx,
			off_t col_idx) const {
		assert(col_idx == 0);
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = (row_idx + i) % f.get_num_levels();
	}
};

class part_apply_operate: public gr_apply_operate<sub_vector_vector>
{
public:
	virtual void run(const void *key, const sub_vector_vector &val,
			mem_vector &vec) const {
	}
	virtual const scalar_type &get_key_type() const {
		return get_scalar_type<factor_value_t>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<size_t>();
	}
	virtual size_t get_num_out_eles() const {
		return 0;
	}
};

void test_groupby()
{
	printf("test groupby\n");
	type_mem_vector<int>::ptr vec = type_mem_vector<int>::create(1000000);
	for (size_t i = 0; i < vec->get_length(); i++)
		vec->set(i, random() % 1000);
	data_frame::ptr res = vec->groupby(adj_apply(), false);
	printf("size: %ld\n", res->get_num_entries());

	vector_vector::ptr vv = vector_vector::cast(res->get_vec("agg"));

	factor f(50);
	factor_vector::ptr labels = factor_vector::create(f, res->get_num_entries());
	labels->set_data(set_label_operate(f));
	labels->sort();
	vv = vv->groupby(*labels, part_apply_operate());
	printf("There are %ld vectors\n", vv->get_num_vecs());
}

int main()
{
	test_groupby();
}
