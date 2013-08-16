#include <linux/kthread.h>
#include <linux/cpu.h>
#include <linux/module.h>	/* Needed by all modules */
#include <linux/kernel.h>	/* Needed for KERN_INFO */
#include <linux/dmaengine.h>
#include <linux/slab.h>

static int version_id = 0;

static atomic_t completed_nthreads;
static struct dma_chan *all_cpu_local_chans[NR_CPUS];

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
	int cpu_id = smp_processor_id();
	all_cpu_local_chans[cpu_id] = dma_find_channel(DMA_MEMCPY);
	pr_info("thread on cpu %d get channal %p on node %d\n", cpu_id,
			all_cpu_local_chans[cpu_id],
			all_cpu_local_chans[cpu_id]->device->dev->numa_node);
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

	get_online_cpus();
	for_each_online_cpu(cpu) {
		struct task_struct *task = kthread_create(&dmachan_gather_func,
				NULL, "dmachan_gather_%d", cpu);
		if (task) {
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
	pr_info("gather all channels\n");

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
	pr_info("There are %d channels and max node id is %d\n", nchans, max_node_id);

	num_chan_vectors = max_node_id + 1;
	chan_vectors = (struct chan_vec *) kmalloc(
			sizeof(struct chan_vec) * num_chan_vectors, GFP_KERNEL);
	memset(chan_vectors, 0, sizeof(struct chan_vec) * num_chan_vectors);
	for (i = 0; i < num_chan_vectors; i++) {
		/* Let's allocate more memory than we need. */
		chan_vectors[i].chans = (struct dma_chan **) kmalloc(
				sizeof(struct dma_chan *) * nchans, GFP_KERNEL);
		memset(chan_vectors[i].chans, 0, sizeof(struct dma_chan *) * nchans);
	}

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
	dmaengine_put();
	printk(KERN_INFO "dmacpy version %d terminating\n", version_id);
}
module_exit(dmacpy_exit);

MODULE_LICENSE("Dual BSD/GPL");
MODULE_AUTHOR("Da Zheng<zhengda1936@gmail.com>");
MODULE_DESCRIPTION("memory copy with DMA engine");
