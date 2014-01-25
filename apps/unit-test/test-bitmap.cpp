#include <algorithm>
#include <set>

#include "bitmap.h"

const int max_bits = 1024 * 1024 * 128;

int main()
{
	int num = 100000;
	bitmap map1(max_bits);
	bitmap map2(max_bits);
	std::set<size_t> elements;

	for (int test = 0; test < 10; test++) {
		for (int i = 0; i < num; i++) {
			int v = random() % max_bits;
			elements.insert(v);
			map1.set(v);
		}
		for (int i = 0; i < num; i++) {
			int v = random() % max_bits;
			elements.insert(v);
			map2.set(v);
		}

		map1.merge(map2);
		std::vector<size_t> bitmap_results;
		map1.get_set_bits(bitmap_results);
		assert(bitmap_results.size() == elements.size());
		std::vector<size_t>::iterator v_it = bitmap_results.begin();
		std::set<size_t>::iterator s_it = elements.begin();
		for (; v_it != bitmap_results.end() && s_it != elements.end(); v_it++, s_it++)
			assert(*v_it == *s_it);
		map1.clear();
		map2.clear();
		elements.clear();
	}
}
