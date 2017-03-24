#include <boost/foreach.hpp>
#include <map>

#include "vector.h"
#include "bulk_operate.h"
#include "data_frame.h"
#include "factor.h"
#include "vector_vector.h"
#include "local_vec_store.h"
#include "local_vv_store.h"
#include "mem_vv_store.h"

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

class part_apply_operate: public gr_apply_operate<local_vv_store>
{
	const scalar_type &val_type;
public:
	part_apply_operate(const scalar_type &type): val_type(type) {
	}

	virtual void run(const void *key, const local_vv_store &val,
			local_vec_store &vec) const {
		assert(val.get_type() == val_type);
		assert(vec.get_type() == get_scalar_type<size_t>());
		vec.resize(1);
		size_t num = 0;
		for (size_t i = 0; i < val.get_length(); i++)
			num += val.get_length(i);
		vec.set<size_t>(0, num);
	}
	virtual const scalar_type &get_key_type() const {
		return get_scalar_type<factor_value_t>();
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
	detail::smp_vec_store::ptr store = detail::smp_vec_store::create(1000000,
			get_scalar_type<int>());
#pragma omp parallel for
	for (size_t i = 0; i < store->get_length(); i++)
		store->set<int>(i, random() % 1000);
	printf("set the vector store\n");
	vector::ptr vec = vector::create(store);
	data_frame::ptr res = vec->groupby(adj_apply(), false);
	printf("size: %ld\n", res->get_num_entries());

	vector_vector::ptr vv = vector_vector::create(
			detail::mem_vv_store::cast(res->get_vec("agg")));

	factor f(50);
	factor_vector::ptr labels = factor_vector::create(f, res->get_num_entries(),
			-1, true, set_label_operate(f, res->get_num_entries()));
	labels->sort();
	vector_vector::ptr gr_res = vv->groupby(*labels,
			part_apply_operate(vv->get_type()));
	const detail::mem_vv_store &res_store
		= dynamic_cast<const detail::mem_vv_store &>(gr_res->get_data());
	printf("There are %ld vectors\n", gr_res->get_num_vecs());

	std::map<factor_value_t, size_t> label_map;
	const detail::smp_vec_store &label_store
		= dynamic_cast<const detail::smp_vec_store &>(labels->get_data());
	for (size_t i = 0; i < labels->get_length(); i++) {
		factor_value_t label = label_store.get<factor_value_t>(i);
		if (label_map.find(label) == label_map.end())
			label_map.insert(std::pair<factor_value_t, size_t>(label, 0));
		label_map[label]++;
	}
	assert(gr_res->get_num_vecs() == label_map.size());

	size_t vec_idx = 0;
	size_t label_idx = 0;
	for (auto it = label_map.begin(); it != label_map.end(); it++) {
		factor_value_t label = it->first;
		size_t num_label_vecs = it->second;
		size_t num_eles = 0;
		for (size_t i = 0; i < num_label_vecs; i++) {
			assert(label_store.get<factor_value_t>(vec_idx) == label);
			num_eles += vv->get_length(vec_idx);
			vec_idx++;
		}
		assert(num_eles == *(size_t *) res_store.get_raw_arr(label_idx));
		label_idx++;
	}
	assert(vec_idx == vv->get_num_vecs());
	assert(label_idx == gr_res->get_num_vecs());
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

void verify_data(const int *buf1, const int *buf2, size_t len)
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
		verify_data((const int *) vv1->get_raw_arr(i),
				(const int *) detail::mem_vec_store::cast(vecs[i])->get_raw_arr(),
				vecs[i]->get_length());
	}

	detail::mem_vv_store::ptr vv2 = detail::mem_vv_store::create(get_scalar_type<int>());
	vv2->append(vecs.begin(), vecs.end());
	assert(vv2->get_num_vecs() == vecs.size());
	for (size_t i = 0; i < vecs.size(); i++) {
		assert(vv2->get_length(i) == vecs[i]->get_length());
		verify_data((const int *) vv2->get_raw_arr(i),
				(const int *) detail::mem_vec_store::cast(vecs[i])->get_raw_arr(),
				vecs[i]->get_length());
	}
}

void test_append_vvs()
{
	printf("test append vector vectors\n");
	std::vector<detail::vec_store::const_ptr> vvs(100);
	std::vector<size_t> vec_lens;
	std::vector<const char *> vec_data;
	for (size_t i = 0; i <vvs.size(); i++) {
		detail::mem_vv_store::ptr vv = create_mem_vv(10, 1000);
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
		verify_data((const int *) vv1->get_raw_arr(i),
				(const int *) vec_data[i], vv1->get_length(i));
	}

	detail::mem_vv_store::ptr vv2 = detail::mem_vv_store::create(get_scalar_type<int>());
	vv2->append(vvs.begin(), vvs.end());
	assert(vv2->get_num_vecs() == vec_lens.size());
	for (size_t i = 0; i < vec_lens.size(); i++) {
		assert(vv2->get_length(i) == vec_lens[i]);
		verify_data((const int *) vv2->get_raw_arr(i),
				(const int *) vec_data[i], vv2->get_length(i));
	}
}

class time2_apply_operate: public arr_apply_operate
{
public:
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
	virtual size_t get_num_out_eles(size_t num_input) const {
		return 0;
	}
};

void test_flatten()
{
	printf("test flatten a sub vector vector\n");
	detail::mem_vv_store::ptr vv_store = create_mem_vv(100, 1000);
	vector_vector::ptr vv = vector_vector::create(vv_store);
	vector::ptr vec = vv->cat();
	assert(memcmp(dynamic_cast<const detail::mem_vec_store &>(
					vec->get_data()).get_raw_arr(),
				vv_store->get_raw_arr(0),
				vec->get_length() * vec->get_entry_size()) == 0);
}

void test_apply()
{
	printf("test apply to each vector\n");
	detail::mem_vv_store::ptr vv_store = create_mem_vv(100, 1000);
	vector_vector::ptr vv = vector_vector::create(vv_store);

	vector_vector::ptr vv2 = vector_vector::cast(vv->apply(time2_apply_operate()));
	assert(vv->get_type() == get_scalar_type<int>());
	assert(vv2->get_type() == get_scalar_type<long>());
	for (size_t i = 0; i < vv->get_length(); i++) {
		assert(vv->get_length(i) == vv2->get_length(i));
		int *vec = (int *) vv_store->get_raw_arr(i);
		long *vec2 = (long *) dynamic_cast<const detail::mem_vv_store &>(
				vv2->get_data()).get_raw_arr(i);
		for (size_t k = 0; k < vv->get_length(i); k++)
			assert(vec[k] * 2 == vec2[k]);
	}

	vector::ptr vec = vv->cat();
	vector::ptr vec2 = vv2->cat();
	printf("vec len: %ld, vec2 len: %ld\n", vec->get_length(), vec2->get_length());
	assert(vec->get_length() == vec2->get_length());
	const detail::smp_vec_store &vstore1
		= dynamic_cast<const detail::smp_vec_store &>(vec->get_data());
	const detail::smp_vec_store &vstore2
		= dynamic_cast<const detail::smp_vec_store &>(vec2->get_data());
	for (size_t i = 0; i < vec->get_length(); i++)
		assert(vstore1.get<int>(i) * 2 == vstore2.get<long>(i));
}

void verify_portion(local_vv_store::const_ptr lvv, detail::mem_vv_store::ptr vv_store)
{
	size_t global_start = lvv->get_global_start();
	for (size_t i = 0; i < lvv->get_length(); i++) {
		const char *lvec = lvv->get_raw_arr(i);
		const char *vec = vv_store->get_raw_arr(global_start + i);
		printf("%ld: %ld, %ld: %ld\n", i, lvv->get_length(i), global_start + i,
				vv_store->get_length(global_start + i));
		assert(lvv->get_length(i) == vv_store->get_length(global_start + i));
		assert(memcmp(lvec, vec,
					lvv->get_length(i) * lvv->get_type().get_size()) == 0);
		lvec = lvv->get_sub_arr(i, i + 1);
		assert(memcmp(lvec, vec,
					lvv->get_length(i) * lvv->get_type().get_size()) == 0);
	}
}

void test_portion()
{
	printf("test getting portions\n");
	detail::mem_vv_store::ptr vv_store = create_mem_vv(100, 1000);
	local_vv_store::ptr lvv = local_vv_store::cast(vv_store->get_portion(9, 7));
	assert(lvv->get_length() == 7);
	assert(lvv->get_global_start() == 9);
	verify_portion(lvv, vv_store);

	printf("test sub vector of a portion\n");
	lvv->expose_sub_vec(2, 3);
	assert(lvv->get_global_start() == 9 + 2);
	assert(lvv->get_length() == 3);
	verify_portion(lvv, vv_store);

	printf("test reset of a portion\n");
	lvv->reset_expose();
	assert(lvv->get_global_start() == 9);
	assert(lvv->get_length() == 7);
	verify_portion(lvv, vv_store);

	printf("test get_portion of a portion\n");
	local_vv_store::ptr llvv = local_vv_store::cast(lvv->get_portion(2, 3));
	assert(llvv->get_length() == 3);
	assert(llvv->get_global_start() == 11);
	verify_portion(llvv, vv_store);

	printf("test a const portion\n");
	detail::mem_vv_store::const_ptr cvv_store = vv_store;
	local_vv_store::const_ptr clvv = local_vv_store::cast(cvv_store->get_portion(9, 7));
	verify_portion(clvv, vv_store);
}

int main()
{
	test_portion();
	test_groupby();
	test_append_vecs();
	test_append_vvs();
	test_flatten();
	test_apply();
}
