#ifndef __MY_STAT_H__
#define __MY_STAT_H__

/**
 * Copyright 2014 Da Zheng
 *
 * This file is part of SA-GraphLib.
 *
 * SA-GraphLib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SA-GraphLib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SA-GraphLib.  If not, see <http://www.gnu.org/licenses/>.
 */

class log_hist_bucket
{
	size_t lower_bound;
	size_t upper_bound;
	size_t count;
public:
	log_hist_bucket() {
		count = 0;
		lower_bound = 0;
		upper_bound = INT_MAX;
	}

	log_hist_bucket(int idx) {
		count = 0;
		lower_bound = pow(10, idx);
		if (idx == 0)
			lower_bound = 0;
		upper_bound = pow(10, idx + 1);
	}

	size_t get_lower_bound() const {
		return lower_bound;
	}

	size_t get_upper_bound() const {
		return upper_bound;
	}

	size_t get_count() const {
		return count;
	}

	void inc_count(size_t num) {
		count += num;
	}
};

class log_histogram
{
	std::vector<log_hist_bucket> buckets;
public:
	log_histogram(int num_buckets): buckets(num_buckets) {
		for (int i = 0; i < num_buckets; i++)
			buckets[i] = log_hist_bucket(i);
	}

	log_hist_bucket &find_bucket(size_t v) {
		for (size_t i = 0; i < buckets.size(); i++)
			if (buckets[i].get_lower_bound() <= v
					&& v < buckets[i].get_upper_bound())
				return buckets[i];
		printf("can't find a bucket for %ld. All buckets cover [%ld, %ld).\n",
				v, buckets.front().get_lower_bound(),
				buckets.back().get_upper_bound());
		assert(0);
	}

	void add_value(size_t v) {
		find_bucket(v).inc_count(1);
	}

	int get_num_buckets() const {
		return buckets.size();
	}

	const log_hist_bucket &get_bucket(int idx) const {
		return buckets[idx];
	}

	void print(FILE *f) const {
		int non_empty = get_num_buckets() - 1;
		for (; non_empty > 0 && get_bucket(non_empty).get_count() == 0; non_empty--);
		for (int i = 0; i <= non_empty; i++) {
			fprintf(f, "[%ld, %ld): %ld\n", get_bucket(i).get_lower_bound(),
					get_bucket(i).get_upper_bound(),
					get_bucket(i).get_count());
		}
	}
};

#endif
