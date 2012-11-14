#!/bin/sh

set_affinity_node ()
{
start_cpu=$1
start_irq=$2
for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15
do
cpu_num=$(($start_cpu + $i * 4))
#cpu_num=$(($start_cpu))
irq_num=$(($start_irq + $i))
echo $cpu_num > /proc/irq/$irq_num/smp_affinity_list
echo "set $irq_num to $cpu_num"
cat /proc/irq/$irq_num/smp_affinity_list
done
}

#echo 1 > /sys/block/sdb/queue/rq_affinity
#echo 1 > /sys/block/sdc/queue/rq_affinity
#echo 1 > /sys/block/sdd/queue/rq_affinity
#echo 1 > /sys/block/sde/queue/rq_affinity
#echo 1 > /sys/block/sdf/queue/rq_affinity
#echo 1 > /sys/block/sdg/queue/rq_affinity
#echo 1 > /sys/block/sdh/queue/rq_affinity
#echo 1 > /sys/block/sdi/queue/rq_affinity
#echo 1 > /sys/block/sdj/queue/rq_affinity
#echo 1 > /sys/block/sdk/queue/rq_affinity
#echo 1 > /sys/block/sdl/queue/rq_affinity
#echo 1 > /sys/block/sdm/queue/rq_affinity
#echo 1 > /sys/block/sdn/queue/rq_affinity

if [ -f /etc/init.d/irq_balancer ] ; then
    /etc/init.d/irq_balancer stop
elif [ -f /etc/init.d/irqbalance ]; then
    /etc/init.d/irqbalance stop
else
    exit 1;
fi;

echo "set for controller 0"
set_affinity_node 0 201
echo "set for controller 1"
set_affinity_node 1 217
echo "set for controller 2"
set_affinity_node 1 233
