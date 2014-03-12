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
	int num_dests = mmsg.get_num_dests();
	multicast_dest_list dest_list = mmsg.get_dest_list();
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
				const vertex_message *msgs[] = {&mmsg};
				info.run_on_messages(graph, msgs, 1);
			}
		}
		if (mmsg.is_activate())
			owner.activate_vertex(id);
	}
}

void message_processor::process_msg(message &msg, bool check_steal)
{
	const int VMSG_BUF_SIZE = 128;
	vertex_message *v_msgs[VMSG_BUF_SIZE];
	while (!msg.is_empty()) {
		int num = msg.get_next(v_msgs, VMSG_BUF_SIZE);
		for (int i = 0; i < num; i++) {
			if (v_msgs[i]->is_multicast()) {
				process_multicast_msg(*multicast_message::cast2multicast(
							v_msgs[i]), check_steal);
				continue;
			}
			vertex_id_t id = v_msgs[i]->get_dest();
			if (!v_msgs[i]->is_empty()) {
				if (check_steal && steal_state->is_stolen(id)) {
					buf_msg(*v_msgs[i]);
				}
				else {
					compute_vertex &info = graph.get_vertex(id);
					info.run_on_messages(graph,
							(const vertex_message **) &v_msgs[i], 1);
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

void message_processor::steal_vertices(vertex_id_t ids[], int num)
{
	steal_state->steal_vertices(ids, num);
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
