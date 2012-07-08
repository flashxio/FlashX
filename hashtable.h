#ifndef __HASHTABLE_H__
#define __HASHTABLE_H__

template<class KeyT, class ValueT>
class hashtable
{
public:
	virtual ValueT get(KeyT key) = 0;

	virtual bool remove(KeyT key, ValueT value) = 0;

	virtual ValueT putIfAbsent(KeyT key, ValueT value) = 0;

	virtual bool replace(KeyT key, ValueT expect, ValueT new_value) = 0;
};

#endif
