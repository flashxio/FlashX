#include <stdlib.h>
#include <assert.h>

#include <utility>
#include <vector>
#include <algorithm>

#include "hilbert_curve.h"

typedef std::pair<off_t, off_t> coordinate_t;

struct coordinate_order
{
	coordinate_t coo;
	size_t order;

	bool operator<(const coordinate_order &o) const {
		return this->order < o.order;
	}
};

std::vector<coordinate_t> gen_coos(size_t n)
{
	std::vector<coordinate_order> orders;
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			coordinate_order o;
			o.coo.first = i;
			o.coo.second = j;
			o.order = hilbert_xy2d(n, i, j);
			orders.push_back(o);
		}
	}
	std::sort(orders.begin(), orders.end());

	std::vector<coordinate_t> coos(orders.size());
	for (size_t i = 0; i < coos.size(); i++)
		coos[i] = orders[i].coo;
	return coos;
}

void verify_coos(std::vector<coordinate_t>::iterator begin,
		std::vector<coordinate_t>::iterator end, size_t n)
{
	if (n == 1)
		return;
	assert(end - begin == n * n);
	assert(n % 2 == 0);
	off_t min_x = std::numeric_limits<off_t>::max();
	off_t min_y = std::numeric_limits<off_t>::max();
	for (std::vector<coordinate_t>::iterator it = begin; it != end; it++) {
		coordinate_t coo = *it;
		min_x = std::min(min_x, coo.first);
		min_y = std::min(min_y, coo.second);
	}
	for (std::vector<coordinate_t>::iterator it = begin; it != end; it++) {
		coordinate_t coo = *it;
		assert(coo.first <= min_x + n);
		assert(coo.second <= min_y + n);
	}
	size_t half_len = n * n / 4;
	verify_coos(begin, begin + half_len, n / 2);
	verify_coos(begin + half_len, begin + half_len * 2, n / 2);
	verify_coos(begin + half_len * 2, begin + half_len * 3, n / 2);
	verify_coos(begin + half_len * 3, begin + half_len * 4, n / 2);
	assert(begin + half_len * 4 == end);
}

void test_size(size_t n)
{
	std::vector<coordinate_t> coos = gen_coos(n);
	assert(coos.size() == n * n);
	verify_coos(coos.begin(), coos.end(), n);
}

int main()
{
	for (size_t i = 1; i <= 64; i *= 2)
		test_size(i);
}
