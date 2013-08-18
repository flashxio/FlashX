#include <linux/kthread.h>
#include <linux/cpu.h>
#include <linux/module.h>	/* Needed by all modules */
#include <linux/kernel.h>	/* Needed for KERN_INFO */
#include <linux/dmaengine.h>
#include <linux/slab.h>

static int version_id = 0;

static atomic_t completed_nthreads;

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
	struct dma_chan **all_cpu_local_chans = data;

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
	atomic_inc(&completed_nthreads);
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

	get_online_cpus();
	for_each_online_cpu(cpu) {
		struct task_struct *task;

		task = kthread_create(&dmachan_gather_func, all_cpu_local_chans,
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
	/*
	 * This isn't a good way to wait, but it guarantees that all threads
	 * have exit and now we have all channels.
	 */
	while (atomic_read(&completed_nthreads) < ncpus);
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

void memcpy_test(void)
{

}

/**************
 * Module core
 */

static int
dmacpy_init(void)
{
	printk(KERN_INFO "dmacpy version %d starting\n", version_id);
	dmaengine_get();
	gather_all_dma_chans();

	return 0;
}
module_init(dmacpy_init);

static void
dmacpy_exit(void)
{
	int i;

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
