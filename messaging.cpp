#include "messaging.h"

//inline int min(int v1, int v2)
//{
//	return v1 > v2 ? v2 : v1;
//}
//
//template<class T>
//int bulk_queue<T>::fetch(T *entries, int num) {
//	pthread_spin_lock(&_lock);
//	int n = min(num, num_entries);
//	for (int i = 0; i < n; i++) {
//		entries[i] = buf[(start + i) % this->size];
//	}
//	start = (start + n) % this->size;
//	num_entries -= n;
//	pthread_spin_unlock(&_lock);
//	return n;
//}
//
//template<class T>
//int bulk_queue<T>::add(T *entries, int num) {
//	pthread_spin_lock(&_lock);
//	int n = min(num, this->size - num_entries);
//	int end = (start + num_entries) % this->size;
//	for (int i = 0; i < n; i++) {
//		buf[(end + i) % this->size] = entries[i];
//	}
//	num_entries += n;
//	pthread_spin_unlock(&_lock);
//	return n;
//}
