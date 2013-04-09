#include "../slab_allocator.h"

int main()
{
	const int num_objs = 1000;
	slab_allocator::linked_obj_list list1;
	slab_allocator::linked_obj_list list2;
	slab_allocator::linked_obj objs[num_objs];
	for (int i = 0; i < num_objs / 2; i++)
		list1.add(&objs[i]);
	for (int i = num_objs / 2; i < num_objs; i++)
		list2.add(&objs[i]);
	printf("There are %d objs in list1\n", list1.get_size());
	printf("There are %d objs in list2\n", list2.get_size());

	slab_allocator::linked_obj_list list3;
	list3.add_list(&list1);
	printf("After moving objs from list 1 to list 3, there are %d objs in list1\n",
			list1.get_size());
	printf("There are %d objs in list3\n", list3.get_size());

	list3.add_list(&list2);
	printf("After moving objs from list 2 to list 3, there are %d objs in list2\n",
			list2.get_size());
	printf("There are %d objs in list3\n", list3.get_size());

	slab_allocator::linked_obj *obj = list3.pop(500);
	int num = 0;
	while (obj != NULL) {
		obj = obj->get_next();
		num++;
	}
	printf("pop 500 objects, there are %d objs in the retuend list\n", num);
	printf("There are %d objs in list 3\n", list3.get_size());

	obj = list3.pop(600);
	num = 0;
	while (obj != NULL) {
		obj = obj->get_next();
		num++;
	}
	printf("pop 600 objects, there are %d objs in the retuend list\n", num);
	printf("There are %d objs in list 3\n", list3.get_size());
}
