#include "disk_read_thread.h"
#include "parameters.h"
#include "container.cpp"

/* just call the callback of the initiator. */
class initiator_callback: public callback
{
public:

	int invoke(io_request *rq) {
		io_interface *io = rq->get_io();
		/* 
		 * after a request is processed,
		 * we need to notify the initiator thread.
		 * It's possible the initiator thread is itself,
		 * we need to stop the infinite loop.
		 */
		if (io->get_callback())
			io->get_callback()->invoke(rq);
		return 0;
	}
};

disk_read_thread::disk_read_thread(const logical_file_partition &partition,
		long size, int node_id): queue(partition.get_file_name(0),
			IO_QUEUE_SIZE, IO_QUEUE_SIZE) {
	aio = new async_io(partition, size, AIO_DEPTH_PER_FILE, node_id);
	aio->set_callback(new initiator_callback());
	this->node_id = node_id;
	num_accesses = 0;

	int ret = pthread_create(&id, NULL, process_requests, (void *) this);
	if (ret) {
		perror("pthread_create");
		exit(1);
	}
}

void disk_read_thread::run() {
	numa_run_on_node(node_id);
	printf("disk read thread runs on node %d\n", node_id);
	aio->init();
	io_request reqs[MAX_FETCH_REQS];
	while (true) {
		/* 
		 * this is the only thread that fetch requests from the queue.
		 * If there are no incoming requests and there are pending IOs,
		 * let's complete the pending IOs first.
		 */
		while (queue.is_empty() && aio->num_pending_IOs() > 0)
			aio->wait4complete();
		int num = queue.fetch(reqs, MAX_FETCH_REQS);
		num_accesses += num;
		aio->access(reqs, num);
	}
	// TODO I need to call cleanup() of aio.
}

void *process_requests(void *arg)
{
	disk_read_thread *thread = (disk_read_thread *) arg;
	printf("disk_read_thread: pid: %d, tid: %ld\n", getpid(), gettid());
	thread->run();
	return NULL;
}
