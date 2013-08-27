#!/bin/bash
# Wrapper script around calibrator. Assumes calibrator and cpumhz are in the current
# directory

MHZ=`./cpumhz`
SIZE=$((13*1048576))
STRIDE=3932160
PREFIX=calibrator-results
TMPFILE=`mktemp`
TOLERANCE=2
MATCH_REQUIREMENT=3
MEASURED=0
FAILED_MEASURE=0

if [ "$TMPFILE" = "" ]; then
	echo ERROR: Failed to create tmpfile
	exit -1
fi
if [ "$MHZ" = "" ]; then
	echo ERROR: Failed to calculate CPU MHz
	exit -1
fi
if [ ! -x ./calibrator ]; then
	echo ERROR: Calibrator is not in current directory and executable
	exit -1
fi
trap "rm $TMPFILE*; exit" INT

MATCHED=0
LAST_LATENCY_CYCLES=-1
LAST_LATENCY_TIME=-1

# Keep increasing size until TLB latency is being measured consistently
while [ $MATCHED -lt $MATCH_REQUIREMENT ]; do
	echo -n "Running calibrator with size $SIZE: "
	./calibrator $MHZ $SIZE $PREFIX > $TMPFILE 2>&1
	if [ $? != 0 ]; then
		SIZE=$(($SIZE*2))
		continue
	fi

	LATENCY_CYCLES=`grep ^TLBs: -A 2 $TMPFILE | tail -1 | awk -F = '{print $2}'`
	LATENCY_CYCLES=`echo $LATENCY_CYCLES | awk '{print $1}'`
	LATENCY_TIME=`grep ^TLBs: -A 2 $TMPFILE | tail -1 | awk '{print $5}'`

	if [ "$LATENCY_CYCLES" = "" ]; then
		echo -n "No TLB Latency measured "
		FAILED_MEASURE=$(($FAILED_MEASURE+1))
		if [ $MEASURED -eq 0 ]; then
			SIZE=$(($SIZE*3/2))
			FAILED_MEASURE=0
		else
			if [ $FAILED_MEASURE -eq 3 ]; then
				SIZE=$(($SIZE+$STRIDE))
				FAILED_MEASURE=0
			else
				echo -n Retrying
			fi
		fi
		echo
		continue
	fi
	echo -n "$LATENCY_CYCLES cycles $LATENCY_TIME ns "
	LOW_TOLERANCE=$(($LATENCY_CYCLES-$TOLERANCE))
	HIGH_TOLERANCE=$(($LATENCY_CYCLES+$TOLERANCE))
	if [ $LAST_LATENCY_CYCLES -ge $LOW_TOLERANCE -a \
			$LAST_LATENCY_CYCLES -le $HIGH_TOLERANCE ]; then
		MATCHED=$(($MATCHED+1))
		echo -n "matched $MATCHED times"
	else
		MATCHED=0
	fi

	LAST_LATENCY_CYCLES=$LATENCY_CYCLES
	LAST_LATENCY_TIME=$LATENCY_TIME
	SIZE=$(($SIZE+$STRIDE))
	MEASURED=$(($MEASURED+1))
	FAILED_MEASURE=0
	echo
done
rm $TMPFILE*

echo
echo TLB_MISS_LATENCY_TIME=$LAST_LATENCY_TIME
echo TLB_MISS_LATENCY_CYCLES=$LAST_LATENCY_CYCLES
echo TLB_MISSES_COST_ONE_SECOND=$(($MHZ*1000000/$LAST_LATENCY_CYCLES))
