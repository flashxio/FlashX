#!/bin/bash

ctl_pattern="mpt.sas"
num_cpus=`grep processor /proc/cpuinfo | wc -l`
num_nodes=`grep "physical id" /proc/cpuinfo | awk 'BEGIN{max=0} {if (max < $4) max=$4} END{print max}'`
num_nodes=$((num_nodes + 1))
echo "There are ${num_cpus} CPUs in ${num_nodes} NUMA nodes"

set_affinity_node ()
{
	start_cpu=$1
	ctl_num=$2
	irqs=`grep ${ctl_pattern}${ctl_num}-msix /proc/interrupts | awk '{print $1}' | awk -F : '{print $1}'`
	i=0
	for irq in $irqs
	do
		cpu_num=$(($start_cpu + $i * $num_nodes))
		irq_num=$irq
		i=$(($i + 1))
		echo $cpu_num > /proc/irq/$irq_num/smp_affinity_list
		echo "set $irq_num to $cpu_num"
		cat /proc/irq/$irq_num/smp_affinity_list
	done
}

if [ -f /etc/init.d/irq_balancer ] ; then
    /etc/init.d/irq_balancer stop
elif [ -f /etc/init.d/irqbalance ]; then
    /etc/init.d/irqbalance stop
else
    exit 1;
fi;

ctl_num=0
node_id=0
while [ 0 ]
do
	irq=`grep ${ctl_pattern}${ctl_num}-msix0 /proc/interrupts`
	if [ $? == 1 ]
	then
		break
	fi

	echo "set for controller ${ctl_num}"
	set_affinity_node $node_id $ctl_num

	ctl_num=$((ctl_num+1))
	node_id=$ctl_num
done
