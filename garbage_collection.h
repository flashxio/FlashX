#ifndef __GARBAGE_COLLECTION_H__
#define __GARBAGE_COLLECTION_H__

class collectable
{
public:
	virtual void set_used(bool used) = 0;
	virtual bool is_used() = 0;
};

template <class T>
class garbage_collection
{
	T *buf;
	int size;
	int start;	// The index of the first used entry
	/*
	 * The number of used entries.
	 * It doesn't indicate the real number of used entries.
	 * Instead, combining with `start', it guarantees to
	 * give us unused entries.
	 * Some entries in the middle might already be set free,
	 * we still can't change it until the first used entry
	 * (pointed by `start') is set free.
	 */
	int num_used;
public:
	garbage_collection(int size) {
		this->size = size;
		buf = new T[size];
		start = 0;
		num_used = 0;
	}

	~garbage_collection() {
		delete [] buf;
	}

	int num_free_entries() {
		return size - num_used;
	}

	T *allocate_obj(int num) {
		if (num_used + num <= size) {
			T *ret = &buf[(start + num_used) % size];
			for (int i = 0; i < num; i++) {
				buf[(start + num_used + i) % size].set_used(true);
			}
			num_used += num;
			return ret;
		}
		else {
			// TODO I need to wait
			return NULL;
		}
	}

	void collect(T *entry) {
		entry->set_used(false);
		if (&buf[start] == entry) {
			/* move to the entry that is still used */
			int tmp = start;
			start = (start + 1) % size;
			num_used--;
			for (; start != tmp && num_used > 0; start = (start + 1) % size) {
				if (buf[start].is_used())
					break;
				num_used--;
			}
		}
		printf("start: %d, num used: %d\n", start, num_used);
	}
};

#endif
