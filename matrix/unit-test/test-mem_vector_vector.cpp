#include <boost/foreach.hpp>

#include "mem_vector.h"
#include "bulk_operate.h"
#include "data_frame.h"
#include "factor.h"
#include "mem_vector_vector.h"
#include "local_vec_store.h"

using namespace fm;

class adj_apply: public gr_apply_operate<local_vec_store>
{
public:
	virtual void run(const void *key, const local_vec_store &val,
			local_vec_store &vec) const {
		vec.resize(val.get_length());
		memcpy(vec.get_raw_arr(), val.get_raw_arr(),
				val.get_length() * val.get_entry_size());
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

class set_label_operate: public type_set_vec_operate<factor_value_t>
{
	factor f;
	size_t num_same_label;
public:
	set_label_operate(const factor &_f, size_t tot_num_eles): f(_f) {
		num_same_label = tot_num_eles / f.get_num_levels();
	}

	virtual void set(factor_value_t *arr, size_t num_eles, off_t start_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = (start_idx + i) / num_same_label;
	}
};

class part_apply_operate: public gr_apply_operate<sub_vector_vector>
{
public:
	virtual void run(const void *key, const sub_vector_vector &val,
			local_vec_store &vec) const {
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
	detail::smp_vec_store::ptr store = detail::smp_vec_store::create(1000000,
			get_scalar_type<int>());
#pragma omp parallel for
	for (size_t i = 0; i < store->get_length(); i++)
		store->set(i, random() % 1000);
	printf("set the vector store\n");
	mem_vector::ptr vec = mem_vector::create(store);
	data_frame::ptr res = vec->groupby(adj_apply(), false);
	printf("size: %ld\n", res->get_num_entries());

	mem_vector_vector::ptr vv = mem_vector_vector::create(
			detail::mem_vv_store::cast(res->get_vec("agg")));

	factor f(50);
	factor_vector::ptr labels = factor_vector::create(f, res->get_num_entries(),
			set_label_operate(f, res->get_num_entries()));
	labels->sort();
	vector_vector::ptr gr_res = vv->groupby(*labels, part_apply_operate());
	printf("There are %ld vectors\n", gr_res->get_num_vecs());
}

detail::vec_store::ptr create_mem_vec(size_t len)
{
	detail::smp_vec_store::ptr v = detail::smp_vec_store::create(
			len, get_scalar_type<int>());
	for (size_t i = 0; i < len; i++)
		v->set<int>(i, random() % 10000);
	return v;
}

detail::mem_vv_store::ptr create_mem_vv(size_t num_vecs, size_t max_vec_len)
{
	std::vector<detail::vec_store::const_ptr> vecs(num_vecs);
	for (size_t i = 0; i < vecs.size(); i++)
		vecs[i] = create_mem_vec(random() % max_vec_len);

	detail::mem_vv_store::ptr vv = detail::mem_vv_store::create(get_scalar_type<int>());
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
	printf("test appending vectors\n");
	std::vector<detail::vec_store::const_ptr> vecs(1000);
	for (size_t i = 0; i < vecs.size(); i++)
		vecs[i] = create_mem_vec(random() % 1000);

	detail::mem_vv_store::ptr vv1 = detail::mem_vv_store::create(get_scalar_type<int>());
	BOOST_FOREACH(detail::vec_store::const_ptr v, vecs)
		vv1->append(*v);
	assert(vv1->get_num_vecs() == vecs.size());
	for (size_t i = 0; i < vecs.size(); i++) {
		assert(vv1->get_length(i) == vecs[i]->get_length());
		verify_data(vv1->get_raw_arr(i), detail::mem_vec_store::cast(vecs[i])->get_raw_arr(),
				vecs[i]->get_length() * vecs[i]->get_entry_size());
	}

	detail::mem_vv_store::ptr vv2 = detail::mem_vv_store::create(get_scalar_type<int>());
	vv2->append(vecs.begin(), vecs.end());
	assert(vv2->get_num_vecs() == vecs.size());
	for (size_t i = 0; i < vecs.size(); i++) {
		assert(vv2->get_length(i) == vecs[i]->get_length());
		verify_data(vv2->get_raw_arr(i), detail::mem_vec_store::cast(vecs[i])->get_raw_arr(),
				vecs[i]->get_length() * vecs[i]->get_entry_size());
	}
}

void test_append_vvs()
{
	printf("test append vector vectors\n");
	std::vector<detail::vec_store::const_ptr> vvs(1000);
	std::vector<size_t> vec_lens;
	std::vector<const char *> vec_data;
	for (size_t i = 0; i <vvs.size(); i++) {
		detail::mem_vv_store::ptr vv = create_mem_vv(100, 1000);
		vvs[i] = vv;
		for (size_t j = 0; j < vv->get_num_vecs(); j++) {
			vec_lens.push_back(vv->get_length(j));
			vec_data.push_back(vv->get_raw_arr(j));
		}
	}

	detail::mem_vv_store::ptr vv1 = detail::mem_vv_store::create(get_scalar_type<int>());
	BOOST_FOREACH(detail::vec_store::const_ptr vv, vvs)
		vv1->append(*vv);
	assert(vv1->get_num_vecs() == vec_lens.size());
	for (size_t i = 0; i < vec_lens.size(); i++) {
		assert(vv1->get_length(i) == vec_lens[i]);
		verify_data(vv1->get_raw_arr(i), vec_data[i],
				vv1->get_length(i) * vv1->get_entry_size());
	}

	detail::mem_vv_store::ptr vv2 = detail::mem_vv_store::create(get_scalar_type<int>());
	vv2->append(vvs.begin(), vvs.end());
	assert(vv2->get_num_vecs() == vec_lens.size());
	for (size_t i = 0; i < vec_lens.size(); i++) {
		assert(vv2->get_length(i) == vec_lens[i]);
		verify_data(vv2->get_raw_arr(i), vec_data[i],
				vv2->get_length(i) * vv2->get_entry_size());
	}
}

class time2_apply_operate: public arr_apply_operate
{
public:
	time2_apply_operate(): arr_apply_operate(0) {
	}
	virtual void run(const local_vec_store &in,
			local_vec_store &out) const {
		out.resize(in.get_length());
		assert(out.get_type() == get_scalar_type<long>());
		for (size_t i = 0; i < out.get_length(); i++) {
			out.set<long>(i, in.get<int>(i) * 2);
			assert(out.get<long>(i) == in.get<int>(i) * 2);
		}
	}
	virtual const scalar_type &get_input_type() const {
		return get_scalar_type<int>();
	}
	virtual const scalar_type &get_output_type() const {
		return get_scalar_type<long>();
	}
};

void test_flatten()
{
	printf("test flatten a sub vector vector\n");
	detail::mem_vv_store::ptr vv_store = create_mem_vv(100, 1000);
	mem_vector_vector::ptr vv = mem_vector_vector::create(vv_store);
	mem_vector::ptr vec = mem_vector::cast(vv->cat());
	assert(memcmp(vec->get_raw_arr(), vv->get_raw_arr(0),
				vec->get_length() * vec->get_entry_size()) == 0);

	detail::mem_vv_store::const_ptr sub_vv = vv_store->get_sub_vec_vec(10, 20);
	mem_vector::ptr sub_vec = mem_vector::create(
			detail::mem_vec_store::cast(sub_vv->cat()));
	off_t sub_off = 0;
	for (int i = 0; i < 10; i++)
		sub_off += vv->get_length(i);
	off_t sub_len = 0;
	for (int i = 10; i < 30; i++)
		sub_len += vv->get_length(i);
	assert(sub_vec->get_length() == sub_len);
	for (size_t i = 0; i < sub_vec->get_length(); i++)
		assert(sub_vec->get<int>(i) == vec->get<int>(i + sub_off));
	assert(memcmp(sub_vec->get_raw_arr(), vv->get_raw_arr(10),
				sub_vec->get_length() * sub_vec->get_entry_size()) == 0);
}

void test_apply()
{
	printf("test apply to each vector\n");
	mem_vector_vector::ptr vv = mem_vector_vector::create(create_mem_vv(100, 1000));
	vector_vector::ptr vv2 = vv->serial_apply(time2_apply_operate());
	assert(vv2->get_length() == vv->get_length());
	assert(vv->get_type() == get_scalar_type<int>());
	assert(vv2->get_type() == get_scalar_type<long>());
	for (size_t i = 0; i < vv->get_length(); i++) {
		assert(vv->get_length(i) == vv2->get_length(i));
		int *vec = (int *) vv->get_raw_arr(i);
		long *vec2 = (long *) vv2->get_raw_arr(i);
		for (size_t k = 0; k < vv->get_length(i); k++)
			assert(vec[k] * 2 == vec2[k]);
	}

	mem_vector::ptr vec = mem_vector::cast(vv->cat());
	mem_vector::ptr vec2 = mem_vector::cast(vv2->cat());
	assert(vec->get_length() == vec2->get_length());
	for (size_t i = 0; i < vec->get_length(); i++)
		assert(vec->get<int>(i) * 2 == vec2->get<long>(i));

	vv2 = vv->apply(time2_apply_operate());
	vec2 = mem_vector::cast(vv2->cat());
	printf("vec len: %ld, vec2 len: %ld\n", vec->get_length(), vec2->get_length());
	assert(vec->get_length() == vec2->get_length());
	for (size_t i = 0; i < vec->get_length(); i++)
		assert(vec->get<int>(i) * 2 == vec2->get<long>(i));
}

int main()
{
	test_groupby();
	test_append_vecs();
	test_append_vvs();
	test_flatten();
	test_apply();
}
