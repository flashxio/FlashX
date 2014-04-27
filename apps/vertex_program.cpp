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

#include "thread.h"

#include "vertex_program.h"
#include "messaging.h"
#include "worker_thread.h"

void vertex_program::multicast_msg(vertex_id_t ids[], int num,
		vertex_message &msg)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	assert(t == curr);
	t->multicast_msg(ids, num, msg);
}

void vertex_program::multicast_msg(edge_seq_iterator &it, vertex_message &msg)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	assert(curr == t);
	t->multicast_msg(it, msg);
}

void vertex_program::send_msg(vertex_id_t dest, vertex_message &msg)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	assert(t == curr);
	t->send_msg(dest, msg);
}

void vertex_program::activate_vertices(vertex_id_t ids[], int num)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	assert(curr == t);
	t->send_activation(ids, num);
}

void vertex_program::activate_vertices(edge_seq_iterator &it)
{
	worker_thread *curr = (worker_thread *) thread::get_curr_thread();
	assert(curr == t);
	t->send_activation(it);
}
