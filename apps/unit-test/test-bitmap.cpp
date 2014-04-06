#include <algorithm>
#include <set>

#include "bitmap.h"

const int max_bits = 1024 * 1024 * 128;

int main()
{
	int num = 100000;
	bitmap map1(max_bits, 0);
	std::set<size_t> elements;

	for (int test = 0; test < 10; test++) {
		for (int i = 0; i < num; i++) {
			int v = random() % max_bits;
			elements.insert(v);
			map1.set(v);
		}
		printf("There are %ld elements\n", elements.size());

		printf("There are %ld set bits\n", map1.get_num_set_bits());
		std::vector<size_t> bitmap_results;
		map1.get_set_bits(bitmap_results);

		std::vector<size_t> bitmap_results1;
		size_t num_bits = map1.get_num_bits();
		const int NUM_CHECK_BITS = 8192;
		for (size_t idx = 0; idx < num_bits; idx += NUM_CHECK_BITS)
			map1.get_set_bits(idx, idx + NUM_CHECK_BITS, bitmap_results1);

		assert(bitmap_results1.size() == elements.size());
		assert(bitmap_results.size() == elements.size());
		std::vector<size_t>::iterator v_it = bitmap_results.begin();
		std::vector<size_t>::iterator v_it1 = bitmap_results1.begin();
		std::set<size_t>::iterator s_it = elements.begin();
		for (; v_it != bitmap_results.end() && v_it1 != bitmap_results1.end()
				&& s_it != elements.end(); v_it++, v_it1++, s_it++) {
			assert(*v_it == *s_it);
			assert(*v_it1 == *s_it);
		}
		assert(v_it == bitmap_results.end());
		assert(v_it1 == bitmap_results1.end());
		assert(s_it == elements.end());
		map1.clear();
		elements.clear();

		bitmap tmp(random() % max_bits, 0);
		tmp.set_all();
		std::vector<size_t> result2;
		tmp.get_set_bits(result2);
	}
}
