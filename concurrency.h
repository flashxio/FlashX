#ifndef __CONCURRENCY_H__
#define __CONCURRENCY_H__

class atomic_integer
{
	volatile int v;
public:
	atomic_integer() {
		v = 0;
	}

	void inc(int by) {
		__sync_fetch_and_add(&v, by);
	}

	void dec(int by) {
		__sync_fetch_and_sub(&v, by);
	}

	int get() {
		return v;
	}
};

class seq_lock
{
	volatile unsigned long count;
	pthread_spinlock_t lock;
public:
	seq_lock() {
		count = 0;
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
	}

	~seq_lock() {
		pthread_spin_destroy(&lock);
	}

	void read_lock(unsigned long &count) {
		/*
		 * If the count is odd, it means another thread is changing
		 * the data structure the lock tries to protect. 
		 */
		do {
			count = this->count;
		} while (count & 1);
	}

	bool read_unlock(unsigned long count) {
		return this->count == count;
	}

	void write_lock() {
		pthread_spin_lock(&lock);
		__sync_fetch_and_add(&count, 1);
	}

	void write_unlock() {
		__sync_fetch_and_add(&count, 1);
		pthread_spin_unlock(&lock);
	}
};

class atomic_flags
{
	volatile int flags;
public:
	atomic_flags() {
		flags = 0;
	}

	void set_flags(int flag) {
		__sync_fetch_and_or(&flags, 0x1 << flag);
	}

	void clear_flags(int flag) {
		__sync_fetch_and_and(&flags, ~(0x1 << flag));
	}

	bool test_flags(int flag) {
		return flags & (0x1 << flag);
	}

	bool test_and_set_flags(int flag) {
		int orig = __sync_fetch_and_or(&flags, 0x1 << flag);
		return orig & (0x1 << flag);
	}
};

#endif
