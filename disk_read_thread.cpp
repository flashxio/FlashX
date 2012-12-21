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

	int invoke(multibuf_io_request *rq) {
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

disk_read_thread::disk_read_thread(const char *name,
		long size): queue(name, IO_QUEUE_SIZE) {
	aio = new async_io(&name, 1, size, AIO_DEPTH_PER_FILE);
	aio->set_callback(new initiator_callback());

	int ret = pthread_create(&id, NULL, process_requests, (void *) this);
	if (ret) {
		perror("pthread_create");
		exit(1);
	}
}

void disk_read_thread::run() {
	aio->init();
	io_request reqs[MAX_FETCH_REQS];
	while (true) {

		/* 
		 * this is the only thread that fetch requests
		 * from the queue.
		 * get all requests in the queue
		 */
		int num = queue.fetch(reqs, MAX_FETCH_REQS);

		aio->access(reqs, num);
	}
	// TODO I need to call cleanup() of aio.
}

#include <sys/types.h>
#include <sys/syscall.h>
#define gettid() syscall(__NR_gettid)

void *process_requests(void *arg)
{
	disk_read_thread *thread = (disk_read_thread *) arg;
	printf("disk_read_thread: pid: %d, tid: %ld\n", getpid(), gettid());
	thread->run();
	return NULL;
}
