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
	message tmp(alloc);
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
