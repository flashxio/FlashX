#include "message_processor.h"
#include "messaging.h"
#include "graph_engine.h"
#include "worker_thread.h"
#include "steal_state.h"

message_processor::message_processor(graph_engine &_graph,
		worker_thread &_owner): graph(_graph),
	owner(_owner), msg_q(_owner.get_node_id(), "graph_msg_queue", 16, INT_MAX),
	stolenv_msgs(_owner.get_node_id(), PAGE_SIZE, true)
{
	steal_state = std::unique_ptr<steal_state_t>(new steal_state_t(graph, owner));
}

void message_processor::buf_msg(vertex_message &vmsg)
{
	if (stolenv_msgs.is_empty() || stolenv_msgs.back().add(vmsg) == NULL) {
		message msg(owner.get_msg_allocator());
		stolenv_msgs.push_back(msg);
		// We have to make sure this is successful.
		assert(stolenv_msgs.back().add(vmsg));
	}
}

/**
 * This class converts a multicast message to a p2p message.
 */
class multicast_p2p_converter
{
	vertex_id_t dest;
	multicast_message &mmsg;
public:
	multicast_p2p_converter(vertex_id_t dest,
			multicast_message &_mmsg): mmsg(_mmsg) {
		this->dest = dest;
	}

	int get_serialized_size() const {
		return mmsg.get_serialized_size() - mmsg.get_num_dests() * sizeof(vertex_id_t);
	}

	int serialize(char *buf, int size) const {
		int serialized_size = this->get_serialized_size();
		assert(serialized_size <= size);
		memcpy(buf, &mmsg, serialized_size);
		vertex_message *vmsg = (vertex_message *) buf;
		*vmsg = vertex_message(serialized_size, mmsg.is_activate());
		vmsg->set_dest(dest);
		return serialized_size;
	}
};

void message_processor::buf_mmsg(vertex_id_t id, multicast_message &mmsg)
{
	multicast_p2p_converter converter(id, mmsg);
	if (stolenv_msgs.is_empty() || stolenv_msgs.back().add(converter) == NULL) {
		message msg(owner.get_msg_allocator());
		stolenv_msgs.push_back(msg);
		// We have to make sure this is successful.
		assert(stolenv_msgs.back().add(converter));
	}
}

void message_processor::process_multicast_msg(multicast_message &mmsg,
		bool check_steal)
{
	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	vertex_program &curr_vprog = t->get_vertex_program();

	int num_dests = mmsg.get_num_dests();
	multicast_dest_list dest_list = mmsg.get_dest_list();

	if (!check_steal) {
		curr_vprog.run_on_multicast_message(graph, mmsg);
		for (int i = 0; i < num_dests; i++) {
			vertex_id_t id = dest_list.get_dest(i);
			if (mmsg.is_activate())
				owner.activate_vertex(id);
		}
		return;
	}

	for (int i = 0; i < num_dests; i++) {
		vertex_id_t id = dest_list.get_dest(i);
		// TODO now the size is the entire message. Now the message
		// is considered as non-empty.
		if (!mmsg.is_empty()) {
			if (check_steal && steal_state->is_stolen(id)) {
				buf_mmsg(id, mmsg);
			}
			else {
				compute_vertex &info = graph.get_vertex(id);
				curr_vprog.run_on_message(graph, info, mmsg);
			}
		}
		if (mmsg.is_activate())
			owner.activate_vertex(id);
	}
}

void message_processor::process_msg(message &msg, bool check_steal)
{
	worker_thread *t = (worker_thread *) thread::get_curr_thread();
	vertex_program &curr_vprog = t->get_vertex_program();

	const int VMSG_BUF_SIZE = 128;
	vertex_message *v_msgs[VMSG_BUF_SIZE];
	while (!msg.is_empty()) {
		int num = msg.get_next(v_msgs, VMSG_BUF_SIZE);
		assert(num > 0);
		// If we aren't in the mode of load balancing and these aren't
		// multicast messages, we can use the fast path.
		// We only need to check the first message. All messages are
		// of the same type.
		if (!check_steal && !v_msgs[0]->is_multicast()) {
			curr_vprog.run_on_messages(graph,
					(const vertex_message **) v_msgs, num);
			for (int i = 0; i < num; i++) {
				vertex_id_t id = v_msgs[i]->get_dest();
				if (v_msgs[i]->is_activate())
					owner.activate_vertex(id);
			}
			continue;
		}

		if (v_msgs[0]->is_multicast()) {
			for (int i = 0; i < num; i++)
				process_multicast_msg(*multicast_message::cast2multicast(
							v_msgs[i]), check_steal);
			continue;
		}

		// If we are here, we are in the load balancing mode.
		assert(check_steal);
		for (int i = 0; i < num; i++) {
			vertex_id_t id = v_msgs[i]->get_dest();
			if (!v_msgs[i]->is_empty()) {
				if (steal_state->is_stolen(id)) {
					buf_msg(*v_msgs[i]);
				}
				else {
					compute_vertex &info = graph.get_vertex(id);
					curr_vprog.run_on_message(graph, info, *v_msgs[i]);
				}
			}
			if (v_msgs[i]->is_activate())
				owner.activate_vertex(id);
		}
	}
}

void message_processor::process_msgs()
{
	if (steal_state->get_num_returned() > 0 && !stolenv_msgs.is_empty()) {
		// TODO we might have to make sure that a lot of messages can be
		// processed. Otherwise, we are wasting time.
		// The vertices have been processed so no other threads will steal
		// them and process them.
		// TODO will it be still true if we allow a vertex to be activated
		// multiple times in the same iteration?
		stack_array<message> msgs(stolenv_msgs.get_num_entries());
		int num_fetched = stolenv_msgs.fetch(msgs.data(),
				stolenv_msgs.get_num_entries());
		for (int i = 0; i < num_fetched; i++)
			process_msg(msgs[i], true);
	}

	const int MSG_BUF_SIZE = 16;
	message msgs[MSG_BUF_SIZE];
	steal_state->guard_msg_processing();
	bool check_steal = steal_state->steal_mode_enabled();
	while (!msg_q.is_empty()) {
		int num_fetched = msg_q.fetch(msgs, MSG_BUF_SIZE);
		for (int i = 0; i < num_fetched; i++)
			process_msg(msgs[i], check_steal);
	}
	steal_state->unguard_msg_processing();
}

void message_processor::steal_vertices(compute_vertex *vertices[], int num)
{
	steal_state->steal_vertices(vertices, num);
}

void message_processor::reset()
{
	steal_state->reset();
	assert(msg_q.is_empty());
	assert(stolenv_msgs.is_empty());
}

void message_processor::return_vertices(vertex_id_t ids[], int num)
{
	steal_state->return_vertices(ids, num);
}

void steal_state_t::steal_vertices(compute_vertex *vertices[], int num)
{
	prepare_steal.fetch_add(1);
	int num_locals = 0;
	for (int i = 0; i < num; i++) {
		int part_id;
		off_t off;
		graph.get_partitioner()->map2loc(vertices[i]->get_id(), part_id, off);
		// It's possible that a vertex that the current thread tries to
		// steal doesn't belong to the owner thread.
		if (part_id == worker_id) {
			stolen_bitmap.set(off);
			num_locals++;
		}
	}
	// TODO I need a memory barrier here.
	// If the guard is odd, it means the owner thread is processing
	// messages. Wait for it to finish processing messages.
	// The thread can't proceed regardless of the state of the owner
	// thread if the owner thread is processing messages.
	while (guard.load() % 2 > 0);
	num_stolen += num_locals;
}

void steal_state_t::return_vertices(vertex_id_t ids[], int num)
{
	num_returned += num;
	for (int i = 0; i < num; i++) {
		int part_id;
		off_t off;
		graph.get_partitioner()->map2loc(ids[i], part_id, off);
		// The vertices have to be returned to the owner thread.
		assert(worker_id == part_id);
		stolen_bitmap.clear(off);
	}
}
