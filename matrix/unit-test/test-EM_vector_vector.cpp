#include "vector_vector.h"
#include "data_frame.h"
#include "EM_vv_store.h"
#include "local_vv_store.h"
#include "local_vec_store.h"
#include "factor.h"
#include "sparse_matrix.h"

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
public:
	virtual void run(const void *key, const local_vv_store &val,
			local_vec_store &vec) const {
		assert(val.get_type() == get_scalar_type<factor_value_t>());
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

vector_vector::ptr create_vv()
{
	detail::smp_vec_store::ptr store = detail::smp_vec_store::create(8000000,
			get_scalar_type<int>());
#pragma omp parallel for
	for (size_t i = 0; i < store->get_length(); i++)
		store->set<int>(i, random() % 1000000);
	vector::ptr vec = vector::create(store);
	data_frame::ptr res = vec->groupby(adj_apply(), false);

	std::vector<detail::vec_store::const_ptr> vvs(1);
	vvs[0] = res->get_vec("agg");
	detail::EM_vv_store::ptr em_vv = detail::EM_vv_store::create(store->get_type());
	em_vv->append(vvs.begin(), vvs.end());
	return vector_vector::create(em_vv);
}

void test_groupby()
{
	printf("test groupby\n");
	vector_vector::ptr vv = create_vv();

	factor f(50);
	factor_vector::ptr labels = factor_vector::create(f, vv->get_num_vecs(), -1,
			false, set_label_operate(f, vv->get_num_vecs()));
	assert(!labels->is_in_mem());
	labels->sort();
	vector_vector::ptr gr_res = vv->groupby(*labels, part_apply_operate());
	printf("There are %ld vectors\n", gr_res->get_num_vecs());

	std::map<factor_value_t, size_t> label_map;
	local_vec_store::const_ptr llabels = labels->get_data().get_portion(0,
			labels->get_length());
	for (size_t i = 0; i < llabels->get_length(); i++) {
		factor_value_t label = llabels->get<factor_value_t>(i);
		if (label_map.find(label) == label_map.end())
			label_map.insert(std::pair<factor_value_t, size_t>(label, 0));
		label_map[label]++;
	}
	assert(gr_res->get_num_vecs() == label_map.size());

	size_t vec_idx = 0;
	size_t label_idx = 0;
	local_vv_store::const_ptr lres = local_vv_store::cast(
			gr_res->get_data().get_portion(0, gr_res->get_num_vecs()));
	for (auto it = label_map.begin(); it != label_map.end(); it++) {
		factor_value_t label = it->first;
		size_t num_label_vecs = it->second;
		size_t num_eles = 0;
		for (size_t i = 0; i < num_label_vecs; i++) {
			assert(llabels->get<factor_value_t>(vec_idx) == label);
			num_eles += vv->get_length(vec_idx);
			vec_idx++;
		}
		assert(num_eles == *(size_t *) lres->get_raw_arr(label_idx));
		label_idx++;
	}
	assert(vec_idx == vv->get_num_vecs());
	assert(label_idx == gr_res->get_num_vecs());
}

int main(int argc, char *argv[])
{
	if (argc < 2) {
		fprintf(stderr, "test conf_file\n");
		exit(1);
	}

	std::string conf_file = argv[1];
	config_map::ptr configs = config_map::create(conf_file);
	init_flash_matrix(configs);

	test_groupby();

	destroy_flash_matrix();
}
