#include <linux/kthread.h>
#include <linux/cpu.h>
#include <linux/module.h>	/* Needed by all modules */
#include <linux/kernel.h>	/* Needed for KERN_INFO */
#include <linux/dmaengine.h>
#include <linux/slab.h>

static int version_id = 0;

struct dmachan_gather_thread_data
{
	struct dma_chan **all_cpu_local_chans;
	wait_queue_head_t wq;
	atomic_t completed_nthreads;
};

/**
 * A vector of DMA channels on each node
 */
struct chan_vec
{
	int nchans;
	struct dma_chan **chans;
};

static struct chan_vec *chan_vectors;
static int num_chan_vectors;

int dmachan_gather_func(void *data)
{
	int node_id;
	int cpu_id = smp_processor_id();
	struct dma_chan *chan;
	struct dmachan_gather_thread_data *thread_data = data;
	struct dma_chan **all_cpu_local_chans = thread_data->all_cpu_local_chans;

	chan= dma_find_channel(DMA_MEMCPY);
	if (chan) {
		node_id = chan->device->dev->numa_node;
		all_cpu_local_chans[cpu_id] = chan;
		pr_info("thread on cpu %d get channal %p on node %d\n", cpu_id,
				chan, node_id);
	}
	else {
		pr_info("can't find a channel on cpu %d\n", cpu_id);
	}
	atomic_inc(&thread_data->completed_nthreads);
	wake_up(&thread_data->wq);
	return 0;
}

/**
 * The current DMA engine breaks in a NUMA machine. The DMA channel
 * attached to the local CPU core may be from another NUMA node.
 * To avoid using a DMA engine on a wrong CPU, we gather all DMA channels
 * first.
 * The only way to gather all channels is to run threads on all CPUs
 * and read the channels on the CPUs.
 */
static void gather_all_dma_chans(void)
{
	int cpu;
	int ncpus = 0;
	int max_node_id = 0;
	int nchans = 0;
	int i;
	struct dma_chan **all_cpu_local_chans;
	struct task_struct **tasks;
	struct dmachan_gather_thread_data thread_data;

	all_cpu_local_chans = kmalloc(sizeof(struct dma_chan *) * NR_CPUS,
			GFP_KERNEL);
	tasks = kmalloc(sizeof(struct task_struct *) * NR_CPUS, GFP_KERNEL);
	if (all_cpu_local_chans == NULL || tasks == NULL) {
		pr_err("can't allocate memory\n");
		kfree(all_cpu_local_chans);
		kfree(tasks);
		return;
	}
	memset(all_cpu_local_chans, 0, sizeof(struct dma_chan *) * NR_CPUS);
	memset(tasks, 0, sizeof(struct task_struct *) * NR_CPUS);

	thread_data.all_cpu_local_chans = all_cpu_local_chans;
	atomic_set(&thread_data.completed_nthreads, 0);
	init_waitqueue_head(&thread_data.wq);

	get_online_cpus();
	for_each_online_cpu(cpu) {
		struct task_struct *task;

		task = kthread_create(&dmachan_gather_func, &thread_data,
				"dmachan_gather_%d", cpu);
		if (task) {
			tasks[ncpus] = task;
			get_task_struct(task);
			kthread_bind(task, cpu);
			wake_up_process(task);
			ncpus++;
		}
	}
	put_online_cpus();

	wait_event_interruptible(thread_data.wq,
			atomic_read(&thread_data.completed_nthreads) == ncpus);
	for (i = 0; i < ncpus; i++) {
		kthread_stop(tasks[i]);
		put_task_struct(tasks[i]);
	}
	kfree(tasks);

	if (all_cpu_local_chans[0])
		nchans++;
	/* Remove all redundant DMA channels. */
	for (cpu = 1; cpu < NR_CPUS; cpu++) {
		int i;

		if (all_cpu_local_chans[cpu] == NULL)
			continue;

		if (all_cpu_local_chans[cpu]->device->dev->numa_node > max_node_id)
			max_node_id = all_cpu_local_chans[cpu]->device->dev->numa_node;

		/* Check with all channels seen previously. */
		for (i = 0; i < cpu; i++) {
			/* If the channel has been seen before, we remove it. */
			if (all_cpu_local_chans[i] == all_cpu_local_chans[cpu]) {
				all_cpu_local_chans[cpu] = NULL;
				break;
			}
		}

		/* If the channel has been seen before. */
		if (i == cpu)
			nchans++;
	}

	/* Initialize the channel vectors. */
	num_chan_vectors = max_node_id + 1;
	chan_vectors = (struct chan_vec *) kmalloc(
			sizeof(struct chan_vec) * num_chan_vectors, GFP_KERNEL);
	if (chan_vectors == NULL) {
		pr_err("can't allocate memory for channel vectors\n");
		goto out;
	}
	memset(chan_vectors, 0, sizeof(struct chan_vec) * num_chan_vectors);
	for (i = 0; i < num_chan_vectors; i++) {
		/* Let's allocate more memory than we need. */
		chan_vectors[i].chans = (struct dma_chan **) kmalloc(
				sizeof(struct dma_chan *) * nchans, GFP_KERNEL);
		if (chan_vectors[i].chans == NULL) {
			pr_err("can't allocate memory for channel vector %d\n", i);
			goto cleanup;
		}
		memset(chan_vectors[i].chans, 0, sizeof(struct dma_chan *) * nchans);
	}

	/* Add channels to the channel vectors. */
	for (i = 0; i < NR_CPUS; i++) {
		int node_id;
		struct dma_chan *chan = all_cpu_local_chans[i];
		if (chan == NULL)
			continue;

		node_id = chan->device->dev->numa_node;
		chan_vectors[node_id].chans[chan_vectors[node_id].nchans++] = chan;
	}
	for (i = 0; i < num_chan_vectors; i++) {
		pr_info("node %d has %d channels\n", i, chan_vectors[i].nchans);
	}
out:
	kfree(all_cpu_local_chans);
	return;

cleanup:
	for (i = 0; i < num_chan_vectors; i++) {
		kfree(chan_vectors[i].chans);
	}
	kfree(chan_vectors);
	chan_vectors = NULL;
	kfree(all_cpu_local_chans);
}

#define NUM_PAGES (128 * 1024)
#define FROM_NODE 0
#define TO_NODE 1
#define NUM_THREADS 8

struct copy_worker_data
{
	void **from_addrs;
	void **to_addrs;
	int num_pages;

	wait_queue_head_t *wq;
	atomic_t *completed_nthreads;
};

int memcpy_worker(void *arg)
{
	int i;
	struct copy_worker_data *data = arg;

	for (i = 0; i < data->num_pages; i++) {
		memcpy(data->to_addrs[i], data->from_addrs[i], PAGE_SIZE);
	}

	atomic_inc(data->completed_nthreads);
	wake_up(data->wq);
	return 0;
}

int cpu_memcpy_test(void *arg)
{
	int i;
	struct page **from_pages = NULL, **to_pages = NULL;
	void **from_addrs = NULL, **to_addrs = NULL;
	struct copy_worker_data *worker_data_arr = NULL;
	struct task_struct **tasks = NULL;
	struct timeval start_time, end_time;
	int npages_per_thread;
	wait_queue_head_t wq;
	atomic_t completed_nthreads;
	int created_nthreads = 0;

	from_pages = kmalloc(sizeof(from_pages[0]) * NUM_PAGES, GFP_KERNEL);
	to_pages = kmalloc(sizeof(to_pages[0]) * NUM_PAGES, GFP_KERNEL);
	from_addrs = kmalloc(sizeof(from_addrs[0]) * NUM_PAGES, GFP_KERNEL);
	to_addrs = kmalloc(sizeof(to_addrs[0]) * NUM_PAGES, GFP_KERNEL);
	if (from_pages == NULL || to_pages == NULL
			|| from_addrs == NULL || to_addrs == NULL) {
		pr_err("can't allocate the page array\n");
		goto cleanup;
	}
	memset(from_pages, 0, sizeof(from_pages[0]) * NUM_PAGES);
	memset(to_pages, 0, sizeof(to_pages[0]) * NUM_PAGES);
	memset(from_addrs, 0, sizeof(from_addrs[0]) * NUM_PAGES);
	memset(to_addrs, 0, sizeof(to_addrs[0]) * NUM_PAGES);

	for (i = 0; i < NUM_PAGES; i++) {
		from_pages[i] = alloc_pages_exact_node(FROM_NODE, GFP_KERNEL, 0);
		if (from_pages[i] == NULL) {
			pr_err("can't allocate pages for from array\n");
			goto cleanup;
		}
		from_addrs[i] = page_address(from_pages[i]);
		to_pages[i] = alloc_pages_exact_node(TO_NODE, GFP_KERNEL, 0);
		if (to_pages[i] == NULL) {
			pr_err("can't allocate pages for to array\n");
			goto cleanup;
		}
		to_addrs[i] = page_address(to_pages[i]);
	}

	init_waitqueue_head(&wq);
	atomic_set(&completed_nthreads, 0);
	worker_data_arr = kmalloc(sizeof(*worker_data_arr) * NUM_THREADS,
			GFP_KERNEL);
	if (worker_data_arr == NULL) {
		pr_err("can't allocate worker data array\n");
		goto cleanup;
	}
	npages_per_thread = NUM_PAGES / NUM_THREADS;
	for (i = 0; i < NUM_THREADS; i++) {
		worker_data_arr[i].from_addrs = from_addrs + npages_per_thread * i;
		worker_data_arr[i].to_addrs = to_addrs + npages_per_thread * i;
		worker_data_arr[i].num_pages = npages_per_thread;
		worker_data_arr[i].wq = &wq;
		worker_data_arr[i].completed_nthreads = &completed_nthreads;
	}

	tasks = kmalloc(sizeof(*tasks) * NUM_THREADS, GFP_KERNEL);
	if (tasks == NULL) {
		pr_err("can't allocate task array\n");
		goto cleanup;
	}
	memset(tasks, 0, sizeof(*tasks) * NUM_THREADS);

	do_gettimeofday(&start_time);
	for (i = 0; i < NUM_THREADS; i++) {
		tasks[i] = kthread_create_on_node(&memcpy_worker, &worker_data_arr[i],
				TO_NODE, "memcpy_%d", i);
		if (tasks[i]) {
			get_task_struct(tasks[i]);
			wake_up_process(tasks[i]);
			created_nthreads++;
		}
		else
			pr_err("fail to create a memcpy worker thread\n");
	}
	wait_event_interruptible(wq,
			atomic_read(&completed_nthreads) == created_nthreads);
	
	do_gettimeofday(&end_time);
	pr_info("cpu memcopy takes %ldus\n", (end_time.tv_sec - start_time.tv_sec)
			* 1000000 + (end_time.tv_usec - start_time.tv_usec));

cleanup:
	if (tasks) {
		for (i = 0; i < NUM_THREADS; i++) {
			if (tasks[i]) {
				kthread_stop(tasks[i]);
				put_task_struct(tasks[i]);
			}
		}
	}
	if (from_pages) {
		for (i = 0; i < NUM_PAGES; i++) {
			if (from_pages[i])
				__free_page(from_pages[i]);
		}
	}
	if (to_pages) {
		for (i = 0; i < NUM_PAGES; i++) {
			if (to_pages[i])
				__free_page(to_pages[i]);
		}
	}
	kfree(from_pages);
	kfree(to_pages);
	kfree(from_addrs);
	kfree(to_addrs);
	kfree(worker_data_arr);
	kfree(tasks);

	return 0;
}

void dmacpy_wait_until(struct dma_chan *chan, dma_cookie_t cookie)
{
	enum dma_status status = DMA_IN_PROGRESS;

	while (status == DMA_IN_PROGRESS) {
		dma_cookie_t done, used;

		status = dma_async_memcpy_complete(chan, cookie, &done, &used);
	}
}

int dma_memcpy_test(void *arg)
{
	int i;
	struct page **from_pages = NULL, **to_pages = NULL;
	struct timeval start_time, end_time;
	struct dma_chan *chan;
	dma_cookie_t last_cookie = 0;
	int chan_node = * (int *) arg;
	int num_cleanup = 0;

	pr_info("dma_memcpy_test runs on cpu %d, use DMA channel on node %d\n",
			smp_processor_id(), chan_node);
	if (chan_node >= num_chan_vectors || chan_vectors[chan_node].nchans <= 0) {
		pr_err("can't get a dma channel\n");
		return -1;
	}
	chan = chan_vectors[chan_node].chans[0];

	from_pages = kmalloc(sizeof(from_pages[0]) * NUM_PAGES, GFP_KERNEL);
	to_pages = kmalloc(sizeof(to_pages[0]) * NUM_PAGES, GFP_KERNEL);
	if (from_pages == NULL || to_pages == NULL) {
		pr_err("can't allocate the page array\n");
		goto cleanup;
	}
	memset(from_pages, 0, sizeof(from_pages[0]) * NUM_PAGES);
	memset(to_pages, 0, sizeof(to_pages[0]) * NUM_PAGES);

	for (i = 0; i < NUM_PAGES; i++) {
		from_pages[i] = alloc_pages_exact_node(FROM_NODE, GFP_KERNEL, 0);
		if (from_pages[i] == NULL) {
			pr_err("can't allocate pages for from array\n");
			goto cleanup;
		}
		to_pages[i] = alloc_pages_exact_node(TO_NODE, GFP_KERNEL, 0);
		if (to_pages[i] == NULL) {
			pr_err("can't allocate pages for to array\n");
			goto cleanup;
		}
	}

	do_gettimeofday(&start_time);
	for (i = 0; i < NUM_PAGES; i++) {
		int err;
		
again:
		err = dma_async_memcpy_pg_to_pg(chan, to_pages[i], 0,
				from_pages[i], 0, PAGE_SIZE);
		if (err == -ENOMEM && last_cookie > 0) {
			num_cleanup++;
			dmacpy_wait_until(chan, last_cookie);
			goto again;
		}
		if (err < 0) {
			pr_info("%d pages have been copied\n", i);
			pr_err("dma err: %d\n", err);
			goto cleanup;
		}
		last_cookie = err;
	}
	dmacpy_wait_until(chan, last_cookie);
	last_cookie = 0;
	do_gettimeofday(&end_time);
	pr_info("dma memcopy takes %ldus, %d cleanups in the middle\n",
			(end_time.tv_sec - start_time.tv_sec)
			* 1000000 + (end_time.tv_usec - start_time.tv_usec), num_cleanup);

cleanup:
	if (last_cookie > 0)
		dmacpy_wait_until(chan, last_cookie);

	if (from_pages) {
		for (i = 0; i < NUM_PAGES; i++) {
			if (from_pages[i])
				__free_page(from_pages[i]);
		}
	}
	if (to_pages) {
		for (i = 0; i < NUM_PAGES; i++) {
			if (to_pages[i])
				__free_page(to_pages[i]);
		}
	}
	kfree(from_pages);
	kfree(to_pages);
	return 0;
}

static struct task_struct *test_thread;

/**************
 * Module core
 */

static int
dmacpy_init(void)
{
	int node_id = 0;
	int (*memcpy_test) (void *data) = dma_memcpy_test;

	printk(KERN_INFO "dmacpy version %d starting\n", version_id);
	dmaengine_get();
	gather_all_dma_chans();

	test_thread = kthread_create_on_node(memcpy_test, &node_id, node_id,
			"memcopy_test");
	if (test_thread) {
		const struct cpumask *cpumask = cpumask_of_node(node_id);

		// Bind the kernel thread to the specified NUMA node.
		get_task_struct(test_thread);
		if (!cpumask_empty(cpumask)) {
			pr_info("bind the testing thread on node %d\n", node_id);
			set_cpus_allowed_ptr(test_thread, cpumask);
		}
		wake_up_process(test_thread);
	}

	return 0;
}
module_init(dmacpy_init);

static void
dmacpy_exit(void)
{
	int i;

	kthread_stop(test_thread);
	put_task_struct(test_thread);

	dmaengine_put();
	if (chan_vectors) {
		for (i = 0; i < num_chan_vectors; i++) {
			kfree(chan_vectors[i].chans);
		}
		kfree(chan_vectors);
	}
	printk(KERN_INFO "dmacpy version %d terminating\n", version_id);
}
module_exit(dmacpy_exit);

MODULE_LICENSE("Dual BSD/GPL");
MODULE_AUTHOR("Da Zheng<zhengda1936@gmail.com>");
MODULE_DESCRIPTION("memory copy with DMA engine");
