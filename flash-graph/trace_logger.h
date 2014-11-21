#ifndef __TRACE_LOGGER_H__
#define __TRACE_LOGGER_H__

/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashGraph.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

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
	thread_safe_FIFO_queue<std::vector<safs::request_range> *> queue;
public:
	log_thread(const std::string &trace_file): thread(
			"trace_log_thread", 0), queue("log_queue", 0, 1024, INT_MAX) {
		f = fopen(trace_file.c_str(), "w");
	}

	void add(std::vector<safs::request_range> *reqs) {
		int ret = queue.add(&reqs, 1);
		if (ret < 1)
			fprintf(stderr, "can't add traced requests to the queue\n");
		activate();
	}

	void run() {
		while (!queue.is_empty()) {
			std::vector<safs::request_range> *reqs = queue.pop_front();
			for (size_t i = 0; i < reqs->size(); i++) {
				safs::request_range req = reqs->at(i);
				fprintf(f, ",%ld,%ld,R,%ld\n",
						req.get_loc().get_offset(), req.get_size(), req.get_size());
			}
			delete reqs;
		}
	}

	void close() {
		fclose(f);
		stop();
		join();
	}
};

const size_t MAX_LOG_BUF = 1024 * 32;

static void destroy_queue(void *p)
{
	std::vector<safs::request_range> *q = (std::vector<safs::request_range> *) p;
	delete q;
}

class trace_logger
{
	log_thread *thread;
	pthread_key_t queue_key;

	std::vector<safs::request_range> *get_per_thread_queue() {
		std::vector<safs::request_range> *p
			= (std::vector<safs::request_range> *) pthread_getspecific(queue_key);
		if (p == NULL) {
			p = new std::vector<safs::request_range>();
			pthread_setspecific(queue_key, p);
		}
		return p;
	}
public:
	typedef std::shared_ptr<trace_logger> ptr;

	trace_logger(const std::string &trace_file) {
		thread = new log_thread(trace_file);
		thread->start();
		pthread_key_create(&queue_key, destroy_queue);
	}

	~trace_logger() {
		close();
	}

	void log(safs::request_range reqs[], int num) {
		std::vector<safs::request_range> *q = get_per_thread_queue();
		q->insert(q->end(), reqs, reqs + num);
		if (q->size() >= MAX_LOG_BUF) {
			thread->add(q);
			pthread_setspecific(queue_key, new std::vector<safs::request_range>());
		}
	}

	void close() {
		// TODO we don't add all requests to the logging thread.
		thread->close();
	}
};

#endif
