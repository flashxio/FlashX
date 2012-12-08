#ifndef __MY_PARAMETERS_H__
#define __MY_PARAMETERS_H__

#define USE_GCLOCK

const int AIO_DEPTH_PER_FILE = 32;

const int IO_QUEUE_SIZE = AIO_DEPTH_PER_FILE * 5;
const int MAX_FETCH_REQS = AIO_DEPTH_PER_FILE;
const int MSG_SEND_BUF_SIZE = 10;

/**
 * The number of requests issued by user applications
 * in one access.
 */
const int NUM_REQS_BY_USER = 100;

/**
 * The initial size of the queue for pending IO requests
 * in the global cache.
 */
const int INIT_GCACHE_PENDING_SIZE = NUM_REQS_BY_USER * 10;

#endif
