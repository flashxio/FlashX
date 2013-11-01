#ifndef __TRACE_LOGGER_H__
#define __TRACE_LOGGER_H__

#include <stdio.h>
#include <pthread.h>

#include <string>
#include <vector>

#include "thread.h"
#include "container.h"

class log_thread: public thread
{
	FILE *f;
	std::string trace_file;
	thread_safe_FIFO_queue<std::vector<io_request> *> queue;
public:
	log_thread(const std::string &trace_file): thread(
			"trace_log_thread", 0), queue(0, 1024) {
		f = fopen(trace_file.c_str(), "w");
	}

	void add(std::vector<io_request> *reqs) {
		int ret = queue.add(&reqs, 1);
		if (ret < 1)
			fprintf(stderr, "can't add traced requests to the queue");
		activate();
	}

	void run() {
		while (!queue.is_empty()) {
			std::vector<io_request> *reqs = queue.pop_front();
			for (size_t i = 0; i < reqs->size(); i++) {
				io_request req = reqs->at(i);
				fprintf(f, ",%ld,%ld,R,%ld\n",
						req.get_offset(), req.get_size(), req.get_size());
			}
		}
	}

	void close() {
		fclose(f);
		stop();
		join();
	}
};

const size_t MAX_LOG_BUF = 1024 * 1024;

class trace_logger
{
	log_thread *thread;
	std::vector<io_request> *logged_reqs;
	pthread_spinlock_t lock;
public:
	trace_logger(const std::string &trace_file) {
		logged_reqs = new std::vector<io_request>();
		pthread_spin_init(&lock, PTHREAD_PROCESS_PRIVATE);
		thread = new log_thread(trace_file);
		thread->start();
	}

	void log(io_request reqs[], int num) {
		pthread_spin_lock(&lock);
		assert(logged_reqs);
		logged_reqs->insert(logged_reqs->end(), reqs, reqs + num);
		if (logged_reqs->size() >= MAX_LOG_BUF) {
			thread->add(logged_reqs);
			logged_reqs = new std::vector<io_request>();
		}
		pthread_spin_unlock(&lock);
	}

	void close() {
		pthread_spin_lock(&lock);
		if (!logged_reqs->empty())
			thread->add(logged_reqs);
		logged_reqs = NULL;
		pthread_spin_unlock(&lock);
		thread->close();
	}
};

#endif
