#include "generic_hashtable.h"

using namespace fm;

int main()
{
	std::vector<int> keys(1000);
	for (size_t i = 0; i < keys.size(); i++)
		keys[i] = random() % 100;
	std::vector<int> vals(1000);
	for (size_t i = 0; i < vals.size(); i++)
		vals[i] = 1;

	std::unordered_map<int, int> counts;
	for (size_t i = 0; i < keys.size(); i++)
		counts.insert(std::pair<int, int>(keys[i], 0));
	for (size_t i = 0; i < keys.size(); i++)
		counts[keys[i]]++;

	bulk_operate::const_ptr add
		= bulk_operate::conv2ptr(get_scalar_type<int>().get_basic_ops().get_add());
	agg_operate::const_ptr sum = agg_operate::create(add);

	generic_hashtable_impl<int, sizeof(int)> ht(get_scalar_type<int>());
	generic_hashtable_impl<int, sizeof(int)> ht2(get_scalar_type<int>());
	ht.insert(keys.size(), keys.data(), vals.data(), *sum);
	data_frame::ptr df = ht.conv2df();
	detail::smp_vec_store::const_ptr key_vec = detail::smp_vec_store::cast(
			df->get_vec(0));
	detail::smp_vec_store::const_ptr val_vec = detail::smp_vec_store::cast(
			df->get_vec(1));
	for (size_t i = 0; i < key_vec->get_length(); i++) {
		int key = key_vec->get<int>(i);
		int val = val_vec->get<int>(i);
		auto ret = counts.find(key);
		assert(ret != counts.end());
		assert(ret->second == val);
	}

	ht2.insert(keys.size(), keys.data(), vals.data(), *sum);
	ht.merge(ht2, *sum);
	df = ht.conv2df();
	key_vec = detail::smp_vec_store::cast(df->get_vec(0));
	val_vec = detail::smp_vec_store::cast(df->get_vec(1));
	for (size_t i = 0; i < key_vec->get_length(); i++) {
		int key = key_vec->get<int>(i);
		int val = val_vec->get<int>(i);
		auto ret = counts.find(key);
		assert(ret != counts.end());
		assert(ret->second * 2 == val);
	}

	return 0;
}
