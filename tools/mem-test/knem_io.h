/*
 * Copyright © inria 2009-2010
 * Brice Goglin <Brice.Goglin@inria.fr>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. The name of the author may not be used to endorse or promote products
 *    derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
 * IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
 * NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
 * THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef KNEM_IO_H
#define KNEM_IO_H 1

#ifdef __KERNEL__
#include <linux/ioctl.h>
#include <linux/types.h>
#else
#include <sys/ioctl.h>
#include <stdint.h>
#endif



/* Pseudo-character device */
#define KNEM_DEVICE_NAME	"knem"
#define KNEM_DEVICE_FILENAME	("/dev/" KNEM_DEVICE_NAME)

/* Current driver binary interface version */
#define KNEM_ABI_VERSION	0x0000000d

/* ioctl base numbers */
#define KNEM_CMD_MAGIC	'K'
#define KNEM_CMD_INDEX	__IOC_NR(x)


/* Offset to give to mmap() when mapping the device file to get an array of status.
 *
 * This memory mapping consists in an array of knem_status_t where
 * asynchronous requests will report their status.
 *
 * The mapping may be as long as desired. If the application may need
 * N asynchronous requests to be pending at the same time, it will need
 * N status slots. N*sizeof(knem_status_t) should therefore be given
 * as a mapping length to mmap().
 *
 * You should also note that the kernel thread associated with this
 * file descriptor is also launched at mmap() since both are related
 * to asynchronous requests. Hence the binding of the thread cannot
 * be changed (with KNEM_CMD_BIND_OFFLOAD) after invoking mmap().
 */
#define KNEM_STATUS_ARRAY_FILE_OFFSET	0



/* ioctl KNEM_CMD_GET_INFO retrieves information about the currently running driver.
 * Takes a struct knem_cmd_info parameter.
 * Nothing is needed in input.
 * The whole struct knem_cmd_info is filled in return.
 *
 * Returns 0 on success. Otherwise returns -1 with the following value in errno:
 * - EFAULT: an invalid pointer was passed.
 */
struct knem_cmd_info {
	uint32_t abi;			/* Driver binary interface version, to be compared with KNEM_ABI_VERSION */
	uint32_t features;		/* Bitmask of features supported by the driver */
	uint32_t ignored_flags;		/* ioctl flags that will be ignored by the driver */
	uint32_t forced_flags;		/* ioctl flags that will be forced by the driver */
};
#define KNEM_CMD_GET_INFO	_IOW(KNEM_CMD_MAGIC, 0x10, struct knem_cmd_info)

/* Feature bits in the features field of struct knem_cmd_info after the KNEM_CMD_GET_INFO ioctl */
#define KNEM_FEATURE_DMA	(1<<0)	/* Offload on DMA Engine is supported, KNEM_FLAG_DMA may be used. */



/* ioctl KNEM_CMD_BIND_OFFLOAD to bind the kernel thread attached to file descriptor.
 * Takes a struct knem_cmd_bind_offload parameter.
 * All fields of the structure are used as input parameters, nothing is modified in return.
 *
 * This ioctl takes precedence over the default kernel module binding policy which is
 * specified by the "binding" module parameter:
 * If set to 1 (default), the kernel thread would be bound to the processor that mmap'ed
 * the status array. If -1, the thread would be bound to anything but this processor.
 * If set to 0, the thread would not be bound.
 *
 * Returns 0 on success. Otherwise returns -1 with the following value in errno:
 * - EBUSY: binding was not possible because the kernel thread was already launched
 *   (the device was mmap'ed before this ioctl).
 * - EFAULT: an invalid pointer was passed.
 */
struct knem_cmd_bind_offload {
	uint32_t flags;		/* Bitmask of binding flags */
	uint32_t mask_len;	/* Binding mask length, in bytes */
	uint64_t mask_ptr;	/* Pointer to the binding mask */
};
#define KNEM_CMD_BIND_OFFLOAD	_IOR(KNEM_CMD_MAGIC, 0x11, struct knem_cmd_bind_offload)

/* Binding flag bits that may be given to the driver in the KNEM_CMD_BIND_OFFLOAD ioctl */
#define KNEM_BIND_FLAG_CUSTOM (1<<0)		/* Bind using given custom mask */
#define KNEM_BIND_FLAG_CURRENT (1<<1)		/* Bind using current process mask */
#define KNEM_BIND_FLAG_CURRENT_REVERSED (1<<2)	/* Bind using current process reversed mask */
typedef uint32_t knem_bind_flags;



/* Flag bits that may be given to the driver in the flags fields.
 *
 * Flags are always honoured. If an unsupported flag is given, the command
 * will fail with EINVAL. The caller should consult the feature mask returned
 * by the KNEM_CMD_GET_INFO command to find out which flags are supported.
 *
 * KNEM_FLAG_DMA is the only flag whose support might be missing for now.
 */
#define KNEM_FLAG_DMA			(1<<0)	/* Offload copies onto DMA engine.
						 * Should only be given when KNEM_FEATURE_DMA is set in the feature mask,
						 * otherwise the command will fail with EINVAL.
						 * Only useful when submitting a copy request.
						 */
#define KNEM_FLAG_ASYNCDMACOMPLETE	(1<<1)	/* Report DMA completion asynchronously.
						 * Do not poll for DMA completion, let it update the status in the background to maximize overlap.
						 * Only useful when submitting a copy request with KNEM_FLAG_DMA also given.
						 */
#define KNEM_FLAG_DMATHREAD		(1<<2)	/* Offload DMA management to a thread.
						 * Have the ioctl return to user-space immediately to allow overlap
						 * while a kernel thread performs the DMA copy management in the background.
						 * Only useful when submitting a copy request with KNEM_FLAG_DMA also given.
						 */
#define KNEM_FLAG_MEMCPYTHREAD		(1<<3)	/* Offload memcpy processing to a thread.
						 * Have the ioctl return to user-space immediately to allow overlap
						 * while a kernel thread performs the copy in the background.
						 * Only useful when submitting a copy request without KNEM_FLAG_DMA.
						 */
#define KNEM_FLAG_PINLOCAL		(1<<4)	/* Always pin the pages of the local process, even in synchronous mode.
						 * For performance debugging purpose only.
						 * Only useful when submitting a copy request without KNEM_FLAG_DMA.
						 */
#define KNEM_FLAG_SINGLEUSE		(1<<5)	/* When creating a region, destroy it after its first use instead of waiting for a destroy ioctl.
						 * Only useful when creating a region.
						 */

#define KNEM_FLAG_ANY_THREAD_MASK	(KNEM_FLAG_DMATHREAD \
					 |KNEM_FLAG_MEMCPYTHREAD)	/* Any kind of offloading to a thread.
									 */
#define KNEM_FLAG_ANY_ASYNC_MASK	(KNEM_FLAG_ASYNCDMACOMPLETE \
					 |KNEM_FLAG_ANY_THREAD_MASK)	/* Any kind of asynchronism, either with a thread
									 * or with asynchronous DMA completion.
									 */
#define KNEM_FLAG_ANY_DMA_MASK		(KNEM_FLAG_DMA \
					 |KNEM_FLAG_ASYNCDMACOMPLETE \
					 |KNEM_FLAG_DMATHREAD)	/* Any DMA-related flag.
								 */

#define KNEM_FLAG_ANY_CREATE_MASK	(KNEM_FLAG_SINGLEUSE)	/* Any flag that may be passed to the create_region ioctl.
								 */
#define KNEM_FLAG_ANY_COPY_MASK		(KNEM_FLAG_DMA \
					 |KNEM_FLAG_ASYNCDMACOMPLETE \
					 |KNEM_FLAG_DMATHREAD \
					 |KNEM_FLAG_MEMCPYTHREAD \
					 |KNEM_FLAG_PINLOCAL)	/* Any flag that may be passed to a copy ioctl.
								 */

typedef uint32_t knem_flags;



/* Cookie identifier returned when creating a region,
 * to be given back to the driver when submitting a copy from the same or another process.
 * Also used to explicitly destroy a region (only from the process that created the region).
 */
typedef uint64_t knem_cookie_t;



/* A contigous segment a data */
struct knem_cmd_param_iovec {
	uint64_t base;	/* Base pointer */
	uint64_t len;	/* Segment size */
};



/* ioctl to declare a persistent memory region.
 * Takes a struct knem_create_region parameter.
 * The cookie is modified and returned to the application,
 * all other fields of the structure are input parameters
 *
 * Returns 0 on success. Otherwise returns -1 with the following value in errno:
 * - ENOMEM: the driver failed to allocate the required memory.
 * - EFAULT: an invalid pointer or memory range was passed.
 */
struct knem_cmd_create_region {
	uint64_t iovec_array;	/* Pointer to the array of source segments (input) */
	uint32_t iovec_nr;	/* Number of source segments (input) */
	uint32_t flags;		/* Region flags bitmask within KNEM_FLAG_ANY_CREATE_MASK (input) */
	uint32_t protection;	/* bitwise-union of PROT_READ and PROT_WRITE (input) */
	uint32_t pad1;
	uint64_t cookie;	/* Region cookie identifier (output) */
};
#define KNEM_CMD_CREATE_REGION	_IOWR(KNEM_CMD_MAGIC, 0x21, struct knem_cmd_create_region)



/* ioctl to destroy a persistent memory region.
 * Takes a knem_cookie_t parameter.
 *
 * Returns 0 on success. Otherwise returns -1 with the following value in errno:
 * - EINVAL: the given cookie does not match an existing region of the current process.
 */
#define KNEM_CMD_DESTROY_REGION	_IOR(KNEM_CMD_MAGIC, 0x22, knem_cookie_t)



/* ioctl to initiate a data transfer between two regions.
 * Takes a struct knem_cmd_copy_bounded parameter. All fields of the
 * structure are used as input parameters, except status.current_status.
 *
 * Bytes are copied until length is reached or until the end of any region
 * is reached.
 *
 * If any asynchronous flag is given (within KNEM_FLAG_ANY_ASYNC_MASK),
 * status.current_status will be set to KNEM_STATUS_PENDING on return.
 * Later, the actual status will be updated in the mmap'ed status array
 * at index given by status.async_status_index. It will change from
 * PENDING to SUCCESS or FAILED once the request will complete in the
 * background.
 *
 * Otherwise, if no asynchronous flag is given, a synchronous request is
 * performed and its completion status is immediately returned in
 * status.current_status. The async_status_index does not need to be valid
 * in this case.
 *
 * This command may be used even for copying between two regions that
 * belong to other processes. Regions may be used multiple times, enabling
 * possible regcache-like user-space optimizations, and also enabling a
 * single region to be used multiple times by different parts of a same
 * collective.
 *
 * Returns 0 on success. Otherwise returns -1 with the following value in errno:
 * - EINVAL: a cookie is invalid.
 * - EINVAL: an asynchronous copy is requested while the status array
 *   has not been mmap'ed yet.
 * - EINVAL: the async_status_index is outside of the array.
 * - EINVAL: an unsupported flag is given, for instance a DMA copy without
 *   any hardware DMA channel.
 * - EACCES: a cookie matches a region whose memory protection does not
 *   allow the requested copy.
 *
 * Once the ioctl returns 0, the copy initialization went fine. Any
 * problem during the actual copy will be reported in the request status,
 * either immediately in current_status or later in the mmap'ed array.
 */
struct knem_cmd_copy_bounded {
	uint64_t src_cookie;		/* Cookie identifier obtained when creating the source region (input) */
	uint64_t src_offset;		/* Offset within the source cookie region (input) */
	uint64_t dst_cookie;		/* Cookie identifier obtained when creating the destination region (input) */
	uint64_t dst_offset;		/* Offset within the destination cookie region (input) */
	uint64_t length;		/* Bytes to copy */
	uint32_t flags;			/* Copy flags bitmask within KNEM_FLAG_ANY_COPY_MASK (input) */
	uint32_t current_status;	/* Current status of the request, PENDING if asynchronous, SUCCESS or FAILED if synchronous (output) */
	uint32_t async_status_index;	/* Index of the box in the mmap'ed status array where the status of the request will be set if asynchronous (input) */
	uint32_t pad;
};
#define KNEM_CMD_COPY_BOUNDED	_IOR(KNEM_CMD_MAGIC, 0x34, struct knem_cmd_copy_bounded)

/* ioctl to copy without explicit boundary.
 * Contrary to knem_copy_bounded which copies up to a given length,
 * this ioctl copies until it reaches the end of one of the regions.
 */
struct knem_cmd_copy {
	uint64_t src_cookie;
	uint64_t src_offset;
	uint64_t dst_cookie;
	uint64_t dst_offset;
	uint32_t flags;
	uint32_t current_status;
	uint32_t async_status_index;
	uint32_t pad;
};
#define KNEM_CMD_COPY	_IOR(KNEM_CMD_MAGIC, 0x32, struct knem_cmd_copy)



/* ioctl to initiate a data transfer between a region and the given local memory segments.
 * Takes a struct knem_cmd_inline_copy_bounded parameter. All fields of the
 * structure are used as input parameters, except status.current_status.
 *
 * To be used when the local segments do not need to be declared as a persistent region.
 * It enables some possible optimizations since they cannot be accessed by other processes.
 *
 * Bytes are copied until length, or until the end of the region is reached, or
 * until the end of local iovecs is reached.
 *
 * If any asynchronous flag is given (within KNEM_FLAG_ANY_ASYNC_MASK),
 * status.current_status will be set to KNEM_STATUS_PENDING on return.
 * Later, the actual status will be updated in the mmap'ed status array
 * at index given by status.async_status_index. It will change from
 * PENDING to SUCCESS or FAILED once the request will complete in the
 * background.
 *
 * Otherwise, if no asynchronous flag is given, a synchronous request is
 * performed and its completion status is immediately returned in
 * status.current_status. The async_status_index does not need to be valid
 * in this case.
 *
 * Return values and behaviors in case of errors during initialization
 * or during the actual copy are similar to those of KNEM_CMD_COPY above.
 */
struct knem_cmd_inline_copy_bounded {
	uint64_t local_iovec_array;	/* Pointer to the array of destination segments (input) */
	uint32_t local_iovec_nr;	/* Number of destination segments (input) */
	uint32_t write;			/* Set if writing to the remote cookie region, unset if reading (input) */
	uint64_t remote_cookie;		/* Cookie identifier obtained when creating the remote region (input) */
	uint64_t remote_offset;		/* Offset within the remote cookie region (input) */
	uint64_t length;		/* Bytes to copy */
	uint32_t flags;			/* Copy flags bitmask within KNEM_FLAG_ANY_COPY_MASK (input) */
	uint32_t current_status;	/* Current status of the request, PENDING if asynchronous, SUCCESS or FAILED if synchronous (output) */
	uint32_t async_status_index;	/* Index of the box in the mmap'ed status array where the status of the request will be set if asynchronous (input) */
	uint32_t pad;
};
#define KNEM_CMD_INLINE_COPY_BOUNDED	_IOR(KNEM_CMD_MAGIC, 0x35, struct knem_cmd_inline_copy_bounded)

/* ioctl to copy without explicit boundary.
 * Contrary to knem_inline_copy_bounded which copies up to a given length,
 * this ioctl copies until it reaches the end of the region or the end
 * of the local iovecs.
 */
struct knem_cmd_inline_copy {
	uint64_t local_iovec_array;
	uint32_t local_iovec_nr;
	uint32_t write;
	uint64_t remote_cookie;
	uint64_t remote_offset;
	uint32_t flags;
	uint32_t current_status;
	uint32_t async_status_index;
	uint32_t pad;
};
#define KNEM_CMD_INLINE_COPY	_IOR(KNEM_CMD_MAGIC, 0x33, struct knem_cmd_inline_copy)



/* Status values */
enum knem_status_e {
	KNEM_STATUS_PENDING,	/* Request is still pending */
	KNEM_STATUS_SUCCESS,	/* Request completed successfully (even if no bytes were copied because of region boundaries) */
	KNEM_STATUS_FAILED	/* Request failed */
};
typedef uint8_t knem_status_t;



/******************************************************
 * Old deprecated API to be dropped in the near future
 */

/*
 * Converting from the old API to the new API:
 * - Replace init_send with create_region (see below).
 * - Replace both init_async_recv and init_sync_recv with inline_copy
 *   (see below).
 * - After every inline_copy ioctl, check the value of current_status.
 *   If it is KNEM_STATUS_SUCCESS or KNEM_STATUS_FAILED, the request is
 *   already completed. If it is KNEM_STATUS_PENDING, the request is still
 *   being processed by the driver in the background. You should poll the
 *   status array until the slot at async_status_index changes from
 *   KNEM_STATUS_PENDING into KNEM_STATUS_SUCCESS or KNEM_STATUS_FAILED.
 */

/* Similar to create_region with protection = PROT_READ (so that the region
 * is only accessible for reading) and flags = KNEM_FLAG_SINGLEUSE (so that
 * the region is only accessible once).
 */
struct knem_cmd_init_send_param {
	uint64_t send_iovec_array;	/* Pointer to the array of source segments (input) */
	uint32_t send_iovec_nr;		/* Number of source segments (input) */
	uint32_t flags;			/* Send bitmask flags, unused yet (input) */
	uint64_t send_cookie;		/* Send cookie identifier (output) */
};
#define KNEM_CMD_INIT_SEND	_IOWR(KNEM_CMD_MAGIC, 0x20, struct knem_cmd_init_send_param)

/* Similar to inline_copy with write = 0, remote_offset = 0,
 * and some asynchronous flags.
 * There was no current_status, the final status was always reported at status_index,
 * even if the request completed synchronously.
 */
struct knem_cmd_init_async_recv_param {
	uint64_t recv_iovec_array;	/* Pointer to the array of destination segments (input) */
	uint32_t recv_iovec_nr;		/* Number of destination segments (input) */
	uint32_t status_index;		/* Index of the box in the mmap'ed status array where the status of the request will be set (input) */
	uint64_t send_cookie;		/* Cookie identifier obtained by the corresponding send ioctl (input) */
	uint32_t flags;			/* Receive bitmask flags (input) */
	uint32_t pad;
};
#define KNEM_CMD_INIT_ASYNC_RECV	_IOR(KNEM_CMD_MAGIC, 0x30, struct knem_cmd_init_async_recv_param)

/* Similar to inline_copy with write = 0, remote_offset = 0,
 * and no asynchronous flag in flags.
 * The final status was always returned in status since all requests
 * were processed synchronously.
 */
struct knem_cmd_sync_recv_param {
	uint64_t recv_iovec_array;	/* Pointer to the array of destination segments (input) */
	uint32_t recv_iovec_nr;		/* Number of destination segments (input) */
	uint32_t status;		/* Status of the request that was just processed (output) */
	uint64_t send_cookie;		/* Cookie identifier obtained by the corresponding send ioctl (input) */
	uint32_t flags;			/* Receive bitmask flags (input) */
	uint32_t pad;
};
#define KNEM_CMD_SYNC_RECV	_IOR(KNEM_CMD_MAGIC, 0x31, struct knem_cmd_sync_recv_param)



#endif /* KNEM_IO_H */

/*
 * Local variables:
 *  tab-width: 8
 *  c-basic-offset: 8
 *  c-indent-level: 8
 * End:
 */
