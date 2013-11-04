/**
 * Copyright 2013 Da Zheng
 *
 * This file is part of SAFSlib.
 *
 * SAFSlib is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SAFSlib is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SAFSlib.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "shadow_cell.h"

#ifdef USE_SHADOW_PAGE

void clock_shadow_cell::add(shadow_page pg)
{
	if (!queue.is_full()) {
		queue.push_back(pg);
		return;
	}
	/*
	 * We need to evict a page from the set.
	 * Find the first page whose reference bit isn't set.
	 */
	bool inserted = false;
	do {
		for (int i = 0; i < queue.size(); i++) {
			last_idx = (last_idx + 1) % queue.size();
			shadow_page old = queue.get(last_idx);
			/* 
			 * The page has been referenced recently,
			 * we should spare it.
			 */
			if (old.referenced()) {
				queue.get(last_idx).set_referenced(false);
				continue;
			}
			queue.set(pg, last_idx);
			inserted = true;
			break;
		}
		/* 
		 * If we can't insert the page in the for loop above,
		 * we need to go through the for loop again.
		 * But for the second time, we will definitely
		 * insert the page.
		 */
	} while (!inserted);
}

shadow_page clock_shadow_cell::search(off_t off)
{
	for (int i = 0; i < queue.size(); i++) {
		shadow_page pg = queue.get(i);
		if (pg.get_offset() == off) {
			queue.get(i).set_referenced(true);
			return pg;
		}
	}
	return shadow_page();
}

void clock_shadow_cell::scale_down_hits()
{
	for (int i = 0; i < queue.size(); i++) {
		queue.get(i).set_hits(queue.get(i).get_hits() / 2);
	}
}

shadow_page LRU_shadow_cell::search(off_t off)
{
	for (int i = 0; i < queue.size(); i++) {
		shadow_page pg = queue.get(i);
		if (pg.get_offset() == off) {
			queue.remove(i);
			queue.push_back(pg);
			return pg;
		}
	}
	return shadow_page();
}

void LRU_shadow_cell::scale_down_hits()
{
	for (int i = 0; i < queue.size(); i++) {
		queue.get(i).set_hits(queue.get(i).get_hits() / 2);
	}
}

/*
 * remove the idx'th element in the queue.
 * idx is the logical position in the queue,
 * instead of the physical index in the buffer.
 */
template<class T, int SIZE>
void embedded_queue<T, SIZE>::remove(int idx)
{
	assert(idx < num);
	/* the first element in the queue. */
	if (idx == 0) {
		pop_front();
	}
	/* the last element in the queue. */
	else if (idx == num - 1){
		num--;
	}
	/*
	 * in the middle.
	 * now we need to move data.
	 */
	else {
		T tmp[num];
		T *p = tmp;
		/* if the end of the queue is physically behind the start */
		if (start + num <= SIZE) {
			/* copy elements in front of the removed element. */
			memcpy(p, &buf[start], sizeof(T) * idx);
			p += idx;
			/* copy elements behind the removed element. */
			memcpy(p, &buf[start + idx + 1], sizeof(T) * (num - idx - 1));
		}
		/* 
		 * the removed element is between the first element
		 * and the end of the buffer.
		 */
		else if (idx + start < SIZE) {
			/* copy elements in front of the removed element. */
			memcpy(p, &buf[start], sizeof(T) * idx);
			p += idx;
			/*
			 * copy elements behind the removed element
			 * and before the end of the buffer.
			 */
			memcpy(p, &buf[start + idx + 1], sizeof(T) * (SIZE - start - idx - 1));
			p += (SIZE - start - idx - 1);
			/* copy the remaining elements in the beginning of the buffer. */
			memcpy(p, buf, sizeof(T) * (num - (SIZE - start)));
		}
		/*
		 * the removed element is between the beginning of the buffer
		 * and the last element.
		 */
		else {
			/* copy elements between the first element and the end of the buffer. */
			memcpy(p, &buf[start], sizeof(T) * (SIZE - start));
			p += (SIZE - start);
			/* copy elements between the beginning of the buffer and the removed element. */
			idx = (idx + start) % SIZE;
			memcpy(p, buf, sizeof(T) * idx);
			p += idx;
			/* copy elements after the removed element and before the last element */
			memcpy(p, &buf[idx + 1], sizeof(T) * ((start + num) % SIZE - idx - 1));
		}
		memcpy(buf, tmp, sizeof(T) * (num - 1));
		start = 0;
		num--;
	}
}

template class embedded_queue<shadow_page, NUM_SHADOW_PAGES>;

#endif

