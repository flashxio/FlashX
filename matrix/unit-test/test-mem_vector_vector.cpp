#include <boost/foreach.hpp>

#include "mem_vector.h"
#include "bulk_operate.h"
#include "data_frame.h"
#include "factor.h"
#include "mem_vector_vector.h"

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

vector::ptr create_mem_vec(size_t len)
{
	type_mem_vector<int>::ptr v = type_mem_vector<int>::create(len);
	for (size_t i = 0; i < len; i++)
		v->set(i, random());
	return v;
}

vector_vector::ptr create_mem_vv(size_t num_vecs, size_t max_vec_len)
{
	std::vector<vector::ptr> vecs(num_vecs);
	for (size_t i = 0; i < vecs.size(); i++)
		vecs[i] = create_mem_vec(random() % max_vec_len);

	type_mem_vector_vector<int>::ptr vv = type_mem_vector_vector<int>::create();
	vv->append(vecs.begin(), vecs.end());
	return vv;
}

void verify_data(const char *buf1, const char *buf2, size_t len)
{
	for (size_t i = 0; i < len; i++)
		assert(buf1[i] == buf2[i]);
}

void test_append_vecs()
{
	std::vector<vector::ptr> vecs(1000);
	for (size_t i = 0; i < vecs.size(); i++)
		vecs[i] = create_mem_vec(random() % 1000);

	type_mem_vector_vector<int>::ptr vv1 = type_mem_vector_vector<int>::create();
	BOOST_FOREACH(vector::ptr v, vecs)
		vv1->append(*v);
	assert(vv1->get_num_vecs() == vecs.size());
	for (size_t i = 0; i < vecs.size(); i++) {
		assert(vv1->get_length(i) == vecs[i]->get_length());
		verify_data(vv1->get_raw_arr(i), mem_vector::cast(vecs[i])->get_raw_arr(),
				vecs[i]->get_length() * vecs[i]->get_entry_size());
	}

	type_mem_vector_vector<int>::ptr vv2 = type_mem_vector_vector<int>::create();
	vv2->append(vecs.begin(), vecs.end());
	assert(vv2->get_num_vecs() == vecs.size());
	for (size_t i = 0; i < vecs.size(); i++) {
		assert(vv2->get_length(i) == vecs[i]->get_length());
		verify_data(vv2->get_raw_arr(i), mem_vector::cast(vecs[i])->get_raw_arr(),
				vecs[i]->get_length() * vecs[i]->get_entry_size());
	}
}

void test_append_vvs()
{
	std::vector<vector::ptr> vvs(1000);
	std::vector<size_t> vec_lens;
	std::vector<const char *> vec_data;
	for (size_t i = 0; i <vvs.size(); i++) {
		vector_vector::ptr vv = create_mem_vv(100, 1000);
		vvs[i] = vv;
		for (size_t j = 0; j < vv->get_num_vecs(); j++) {
			vec_lens.push_back(vv->get_length(j));
			vec_data.push_back(vv->get_raw_arr(j));
		}
	}

	type_mem_vector_vector<int>::ptr vv1 = type_mem_vector_vector<int>::create();
	BOOST_FOREACH(vector::ptr vv, vvs)
		vv1->append(*vv);
	assert(vv1->get_num_vecs() == vec_lens.size());
	for (size_t i = 0; i < vec_lens.size(); i++) {
		assert(vv1->get_length(i) == vec_lens[i]);
		verify_data(vv1->get_raw_arr(i), vec_data[i],
				vv1->get_length(i) * vv1->get_entry_size());
	}

	type_mem_vector_vector<int>::ptr vv2 = type_mem_vector_vector<int>::create();
	vv2->append(vvs.begin(), vvs.end());
	assert(vv2->get_num_vecs() == vec_lens.size());
	for (size_t i = 0; i < vec_lens.size(); i++) {
		assert(vv2->get_length(i) == vec_lens[i]);
		verify_data(vv2->get_raw_arr(i), vec_data[i],
				vv2->get_length(i) * vv2->get_entry_size());
	}
}

int main()
{
	test_groupby();
	test_append_vecs();
	test_append_vvs();
}
