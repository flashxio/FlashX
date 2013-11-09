#include <stdio.h>
#include <unistd.h>

#include "timer.h"

class test_task: public timer_task
{
	long prev;
	int task_id;
public:
	test_task(int task_id): timer_task(10, true) {
		prev = get_curr_time_us();
		this->task_id = task_id;
	}

	void run() {
		long curr = get_curr_time_us();
		printf("task %d: timeout after %ldus\n", task_id, curr - prev);
		prev = curr;
	}
};

int main()
{
	timer t;

	for (int i = 1; i <= 10; i++)
		t.add_task(new test_task(i));
	sleep(100);
}
