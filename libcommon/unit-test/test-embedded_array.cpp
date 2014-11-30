#include "container.h"

int main()
{
	embedded_array<int> arr1, arr2;
	arr2.resize(100);
	arr2 = arr1;

	embedded_array<int> arr3, arr4;
	arr3.resize(100);
	arr4.resize(100);
	arr3 = arr4;
}
