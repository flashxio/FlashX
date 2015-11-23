#include <map>

#include "data_frame.h"
#include "bulk_operate.h"
#include "vector.h"
#include "vector_vector.h"
#include "local_vec_store.h"
#include "EM_vector.h"
#include "sparse_matrix.h"

using namespace fm;

class rand_set_vec: public type_set_vec_operate<int>
{
public:
	virtual void set(int *arr, size_t num_eles, off_t start_idx) const {
		for (size_t i = 0; i < num_eles; i++)
			arr[i] = random() % 1000;
	}
};

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

void test_in_mem_groupby()
{
	printf("test in-mem groupby\n");
	data_frame::ptr df = data_frame::create();
	size_t length = 1000000;
	detail::smp_vec_store::ptr vec1 = detail::smp_vec_store::create(length,
			get_scalar_type<int>());
	vec1->set_data(rand_set_vec());
	detail::smp_vec_store::ptr vec2 = detail::smp_vec_store::create(length,
			get_scalar_type<int>());
	vec2->set_data(rand_set_vec());

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
	vector::ptr v = vv->cat();
	const detail::smp_vec_store &vstore
		= dynamic_cast<const detail::smp_vec_store &>(v->get_data());
	assert(map.size() == v->get_length());
	size_t idx = 0;
	for (auto it = map.begin(); it != map.end(); it++, idx++) {
		assert(it->second == vstore.get<long>(idx));
	}

	copy_apply_operate copy_op;
	vv = df->groupby("vec1", copy_op);
	vector::ptr v1 = vv->cat();
	assert(vec1->get_length() == v1->get_length());
}

void test_EM_groupby()
{
	printf("test EM groupby\n");
	data_frame::ptr df = data_frame::create();
	size_t length = 4000000;
	detail::vec_store::ptr vec1 = detail::vec_store::create(length,
			get_scalar_type<int>(), -1, false);
	vec1->set_data(rand_set_vec());
	detail::vec_store::ptr vec2 = detail::vec_store::create(length,
			get_scalar_type<int>(), -1, false);
	vec2->set_data(rand_set_vec());
	df->add_vec("vec1", vec1);
	df->add_vec("vec2", vec2);

	sum_apply_operate sum_op;
	vector_vector::ptr vv = df->groupby("vec1", sum_op);

	copy_apply_operate copy_op;
	vv = df->groupby("vec1", copy_op);
	vector::ptr v1 = vv->cat();
	assert(vec1->get_length() == v1->get_length());
}

void test_in_mem_sort()
{
	printf("test in-mem sort\n");
	detail::smp_vec_store::ptr vec1 = detail::smp_vec_store::create(10000,
			get_scalar_type<int>());
	vec1->set_data(rand_set_vec());
	detail::smp_vec_store::ptr vec2 = detail::smp_vec_store::cast(vec1->deep_copy());
	for (size_t i = 0; i < vec2->get_length(); i++)
		vec2->set<int>(i, vec2->get<int>(i) * 2);
	std::vector<named_vec_t> vecs(2);
	vecs[0].first = "1";
	vecs[0].second = vec1;
	vecs[1].first = "2";
	vecs[1].second = vec2;
	data_frame::ptr df = data_frame::create(vecs);
	data_frame::const_ptr sorted = df->sort("2");
	assert(sorted->get_vec_name(0) == "1");
	assert(sorted->get_vec_name(1) == "2");
	assert(sorted->get_vec(0)->is_sorted());
	assert(sorted->get_vec(1)->is_sorted());
}

void test_EM_sort()
{
	printf("test EM sort\n");
	detail::vec_store::ptr vec1 = detail::EM_vec_store::create(1024 * 1024,
			get_scalar_type<int>());
	vec1->set_data(rand_set_vec());
	// TODO we should use a vector with different values to test it.
	detail::vec_store::ptr vec2 = vec1->deep_copy();
	std::vector<named_vec_t> vecs(2);
	vecs[0].first = "1";
	vecs[0].second = vec1;
	vecs[1].first = "2";
	vecs[1].second = vec2;
	data_frame::ptr df = data_frame::create(vecs);
	data_frame::const_ptr sorted = df->sort("2");
	assert(sorted->get_vec_name(0) == "1");
	assert(sorted->get_vec_name(1) == "2");
	assert(sorted->get_vec(0)->is_sorted());
	assert(sorted->get_vec(1)->is_sorted());
}

data_frame::ptr create_data_frame(bool in_mem)
{
	detail::vec_store::ptr vec1 = detail::vec_store::create(1024 * 1024,
			get_scalar_type<int>(), -1, in_mem);
	vec1->set_data(rand_set_vec());
	detail::vec_store::ptr vec2 = detail::vec_store::create(1024 * 1024,
			get_scalar_type<int>(), -1, in_mem);
	vec2->set_data(rand_set_vec());
	std::vector<named_vec_t> vecs(2);
	vecs[0].first = "1";
	vecs[0].second = vec1;
	vecs[1].first = "2";
	vecs[1].second = vec2;
	return data_frame::create(vecs);
}

void test_merge()
{
	printf("test merging multiple data frames\n");
	std::vector<data_frame::const_ptr> dfs;
	size_t tot_num_entries = 0;
	for (size_t i = 0; i < 5; i++) {
		data_frame::ptr df = create_data_frame(true);
		tot_num_entries += df->get_num_entries();
		dfs.push_back(df);
	}
	data_frame::ptr res = merge_data_frame(dfs, true);
	assert(res->get_num_vecs() == dfs.front()->get_num_vecs());
	assert(res->get_num_entries() == tot_num_entries);

	res = merge_data_frame(dfs, false);
	assert(res->get_num_vecs() == dfs.front()->get_num_vecs());
	assert(res->get_num_entries() == tot_num_entries);
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

	test_merge();
	test_in_mem_groupby();
	test_EM_groupby();
	test_in_mem_sort();
	test_EM_sort();

	destroy_flash_matrix();
}
