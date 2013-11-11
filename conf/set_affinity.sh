#!/bin/bash

ctl_name="mpt2sas"

set_affinity_node ()
{
	start_cpu=$1
	start_irq=$2
	for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
	do
		cpu_num=$(($start_cpu + $i * 4))
#		cpu_num=$(($cpu_num % 32))
		#cpu_num=$(($start_cpu))
		irq_num=$(($start_irq + $i))
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
	irq=`grep ${ctl_name}${ctl_num}-msix0 /proc/interrupts | awk '{print $1}' | awk -F : '{print $1}'`
	if [ "$irq" == "" ]
	then
		break
	fi
	echo "irq: $irq"

	echo "set for controller ${ctl_num}"
	set_affinity_node $node_id $irq

	ctl_num=$((ctl_num+1))
	node_id=$ctl_num
done
