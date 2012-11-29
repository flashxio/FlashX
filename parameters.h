#ifndef __MY_PARAMETERS_H__
#define __MY_PARAMETERS_H__

#define USE_GCLOCK

const int AIO_DEPTH_PER_FILE = 32;

const int IO_QUEUE_SIZE = AIO_DEPTH_PER_FILE * 5;
const int MAX_FETCH_REQS = AIO_DEPTH_PER_FILE;
const int MSG_SEND_BUF_SIZE = 10;

#endif
