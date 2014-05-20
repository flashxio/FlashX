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

#include "messaging.h"

int multicast_msg_sender::flush()
{
	if (buf.is_empty()) {
		assert(mmsg == NULL);
		return 0;
	}
	this->mmsg = NULL;
	this->num_dests = 0;
	dest_list.clear();
	queue->add(&buf, 1);
	if (buf.get_num_objs() > 1)
		printf("there are %d objs in the msg\n", buf.get_num_objs());
	// We have to make sure all messages have been sent to the queue.
	assert(buf.is_empty());
	message tmp(alloc.get());
	buf = tmp;
	return 1;
}

bool multicast_msg_sender::add_dest(local_vid_t id)
{
	int ret = buf.inc_msg_size(sizeof(id.id));
	if (ret == 0) {
		flush();

		multicast_message *mmsg_template
			= (multicast_message *) mmsg_temp_buf;
		vertex_message *p = (vertex_message *) buf.add(*mmsg_template);
		// We just add the buffer. We should be able to add the new message.
		assert(p);
		this->mmsg = multicast_message::convert2multicast(p);
		dest_list = this->mmsg->get_dest_list();
		buf.inc_msg_size(sizeof(id.id));
	}
	num_dests++;
	dest_list.add_dest(id);
	return true;
}

int multicast_msg_sender::add_dests(local_vid_t ids[], int num)
{
	int orig_num = num;
	while (num > 0) {
		int num_allowed = min(buf.get_remaining_size(),
				num * sizeof(ids[0].id)) / sizeof(ids[0].id);
		if (num_allowed == 0) {
			flush();

			multicast_message *mmsg_template
				= (multicast_message *) mmsg_temp_buf;
			vertex_message *p = (vertex_message *) buf.add(*mmsg_template);
			// We just add the buffer. We should be able to add the new message.
			assert(p);
			this->mmsg = multicast_message::convert2multicast(p);
			dest_list = this->mmsg->get_dest_list();

			num_allowed = min(buf.get_remaining_size(),
					num * sizeof(ids[0].id)) / sizeof(ids[0].id);
			assert(num_allowed > 0);
		}
		buf.inc_msg_size(num_allowed * sizeof(ids[0].id));
		num_dests += num_allowed;
		dest_list.add_dests(ids, num_allowed);
		ids += num_allowed;
		num -= num_allowed;
		assert(num >= 0);
	}
	return orig_num;
}
