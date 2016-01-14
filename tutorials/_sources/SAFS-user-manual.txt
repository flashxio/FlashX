SAFS is user space filesystem designed to access a large SSD array. The
goal of SAFS is to maximize the I/O performance of the SSD array on a
NUMA machine while still providing a filesystem interface to users. It
provides basic operations on files: create, delete, read and write. A
file exposed by SAFS is partitioned and each partition is stored as a
physical file on an SSD. SAFS currently does not support directory
operations. SAFSlib is a user-space C++ library that implements SAFS.

Programming interface
=====================

File metadata operations
------------------------

The class ``safs_file`` represents an SAFS file and provides a few
methods for some metadata operations such as creating a file, deleting a
file and renaming a file.

.. raw:: html

   <pre>
   class safs_file
   {
   public:
       /* The constructor method. The file doesn't need to exist. */
       safs_file(const RAID_config &conf, const std::string &file_name);
       /* Test whether the SAFS file exists. */
       bool exist() const;
       /* Get the size of the SAFS file. */
       ssize_t get_size() const;
       /* Create the SAFS file with the specified size. */
       bool create_file(size_t file_size);
       /* Delete the SAFS file. */
       bool delete_file();
       /* Rename the SAFS file to a new name. */
       bool rename(const std::string &new_name);
   };
   </pre>

SAFS does not support directories. The function ``get_all_safs_files``
returns all files in SAFS.

.. raw:: html

   <pre>
   size_t get_all_safs_files(std::set<std::string> &files);
   </pre>

File access
-----------

Two classes (``file_io_factory`` and ``io_interface``) are used for
accessing data in a file. The class ``file_io_factory`` creates and
destroys ``io_interface`` objects, which provides methods to read and
write an SAFS file. An ``io_interface`` instance can only access a
single file and can only be used in a single thread. We intentionally
make the implementations of ``io_interface`` **not thread-safe** for the
sake of performance.

File open and close
~~~~~~~~~~~~~~~~~~~

The function ``create_io_factory`` creates a ``file_io_factory``
instance for a file. It allows a user to specify an access option, which
decides what type of the ``file_io_factory`` instance is created. Right
now, SAFS supports two access options:

-  REMOTE\_ACCESS: this corresponds to direct I/O in Linux. The
   ``io_interface`` instance created by such a ``file_io_factory``
   doesn't use page cache in SAFS.
-  GLOBAL\_CACHE\_ACCESS: this corresponds to buffered I/O in Linux. The
   ``io_interface`` instance uses page cache in SAFS.

Opening a file involves in two steps: invoking ``create_io_factory`` to
create a ``file_io_factory`` object; invoking the ``create_io`` method
of ``file_io_factory`` to create an ``io_interface`` object. Files are
closed implicitly when the ``file_io_factory`` object is destroyed.

Synchronous read and write
~~~~~~~~~~~~~~~~~~~~~~~~~~

| A user can use the following method of ``io_interface`` to issue
  synchronous I/O requests. ``access_method`` determines whether it is a
  read or write request: 0 indicates read and 1 indicates write.
| 
| class io\_interface
| {
| public:
| io\_status access(char \*buf, off\_t off, ssize\_t size, int
  access\_method);
| ...
| };
| 

Asynchronous read and write
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users can use the following set of methods to use asynchronous I/O.
First, users need to implement the ``callback`` interface and register
it to an ``io_interface`` object to get notification of completion of
I/O requests before issuing any I/O requests. Then, they use the
asynchronous version of the ``access`` method to issue I/O requests.
When a request completes, the ``callback`` is invoked. It is guaranteed
that the ``callback`` object will be invoked in the same thread where
the object is associated to. An ``io_interface`` instance does not limit
the number of parallel I/O requests that can be issued to it. Users need
to monitor and limit the number of incomplete I/O requests themselves.

.. raw:: html

   <pre>
   class io_interface
   {
   public:
       ...
       /* Issue asynchronous I/O requests. */
       void access(io_request *, int, io_status *);
       /* Flash I/O requests buffered by the io_interface instance. */
       void flush_requests();
       /* Wait for at least the specified number of I/O requests to complete. */
       int wait4complete(int);
       /* Get the number of pending I/O requests. */
       int num_pending_ios() const;
       /* set the callback function. */
       bool set_callback(callback *);
   };
   </pre>

.. raw:: html

   <pre>
   class callback
   {
   public:
       virtual int invoke(io_request *reqs[], int num) = 0;
   };
   </pre>

A simple example of using the library
-------------------------------------

The following pseudocode illustrates a simple use case of SAFS, which
uses its synchronous I/O interface to read data from a file.

.. raw:: html

   <pre>
   config_map::ptr configs = config_map::create(conf_file);
   init_io_system(configs);
   file_io_factory::shared_ptr factory = create_io_factory(graph_file,
           REMOTE_ACCESS);
   io_interface::ptr io = factory->create_io(
           thread::get_curr_thread());

   char *buf = NULL;
   size_t buf_capacity = 0;
   BOOST_FOREACH(task t, tasks) {
       // This is directed I/O. Memory buffer, I/O offset and I/O size
       // all need to be aligned to the I/O block size.
       size_t io_size = ROUNDUP_PAGE(t.get_size());
       data_loc_t loc(factory->get_file_id(), t.get_offset());
       if (io_size > buf_capacity) {
           free(buf);
           buf_capacity = io_size;
           buf = (char *) valloc(buf_capacity);
       }
       assert(buf_capacity >= io_size);
       io_request req(buf, loc, io_size, READ);
       io->access(&req, 1);
       io->wait4complete(1);
       run_computation(buf, t.get_size());
   }
   free(buf);
   </pre>

The following pseudocode illustrates a use case of SAFS' asynchronous
I/O interface to read data from a file. It is slightly more complex than
the synchronous I/O interface. It requires to define a callback class
and the computation is performed in the callback class.

.. raw:: html

   <pre>
   class multiply_callback: public callback
   {
       const std::vector<task> &tasks;
   public:
       multiply_callback(const std::vector<task> &_tasks): tasks(_tasks) {
       }

       virtual int invoke(io_request *reqs[], int num);
   };

   struct comp_task
   {
       bool operator()(const task &t1, const task &t2) const {
           return t1.get_offset() < t2.get_offset();
       }
   };

   int multiply_callback::invoke(io_request *reqs[], int num)
   {
       for (int i = 0; i < num; i++) {
           off_t off = reqs[i]->get_offset();
           ext_mem_vertex_info info(0, off, 0);
           std::vector<task>::const_iterator it = std::lower_bound(
                   tasks.begin(), tasks.end(),
                   task(info), comp_task());
           assert(it != tasks.end());
           assert(it->get_offset() == off);
           char *buf = reqs[i]->get_buf();
           run_multiply(buf, it->get_size());
           free(buf);
       }
       return 0;
   }

   config_map::ptr configs = config_map::create(conf_file);
   init_io_system(configs);
   file_io_factory::shared_ptr factory = create_io_factory(graph_file,
           REMOTE_ACCESS);
   io_interface::ptr io = factory->create_io(thread::get_curr_thread());
   io->set_callback(new multiply_callback(tasks));

   int max_ios = 20;
   BOOST_FOREACH(task t, tasks) {
       while (io->num_pending_ios() >= max_ios)
           io->wait4complete(1);

       // This is directed I/O. Memory buffer, I/O offset and I/O size
       // all need to be aligned to the I/O block size.
       size_t io_size = ROUNDUP_PAGE(t.get_size());
       data_loc_t loc(factory->get_file_id(), t.get_offset());
       char *buf = (char *) valloc(io_size);
       io_request req(buf, loc, io_size, READ);
       io->access(&req, 1);
   }
   io->wait4complete(io->num_pending_ios());
   </pre>

Utility tool in SAFS
--------------------

SAFS-util is a tool that helps to manage SAFS. It provides a few
commands to operate SAFS:

-  create: create a file in SAFS.
-  delete file\_name: delete a file in SAFS.
-  list: list all existing files in SAFS.
-  load: load a file from an external filesystem to a file in SAFS.
-  verify: verify the data of a file in SAFS. It’s mainly used for
   testing.

Configurations
==============

System configurations for SAFS:
-------------------------------

SAFS requires proper system configurations to get the maximal
performance from an SSD array. It includes evenly distributing
interrupts to all CPU cores and using the noop I/O scheduler;

The following two scripts automate the process.

-  ``conf/set_affinity.sh``: The script disables irqbalance and
   distributes IRQs (Interrupt Requests) to CPU cores evenly. It is only
   required when we run SAFS on a NUMA machine, and it is written for an
   LSI host bus adapter (HBA). For other HBAs, users need to adapt the
   script for their specific hardware.
-  ``conf/set_ssds.pl``: The script takes an input file that contains
   device files to run SAFS on. The input file has a device file on each
   line. The script identifies the SSDs that connect to the same HBA
   controller, and assigns all SSDs in the same HBA to the same NUMA
   node. It sets I/O request affinity for each device file to ``2`` to
   force the request completion on the requesting CPU core. It
   configures each device file to use the noop I/O scheduler. Then it
   mounts SSDs to the system, and creates conf/data\_files.txt, which is
   used as the configuration file of the root directories by the
   library. In conf/data\_files.txt, each line has the format:
   ``node_id:path_to_mountpoint``.

The configurations in SAFS
--------------------------

SAFS defines the following parameters for users to customize SAFS. When
SAFS is initialized, users have the opportunity to set them.

-  ``root_conf``: a config file that specifies the directories on SSDs
   where SAFS files are stored. The config file has a line for each
   directory and the format of each line is
   ``node_id:path_of_directory``.
-  ``RAID_block_size``: defines the size of a data block on an SSD.
   Users can specify the size with the format x(k, K, m, M, g, G). e.g.,
   4k = 4096 bytes. The default block size is 512KB.
-  ``RAID_mapping``: defines how data blocks of a file are mapped to
   SSDs. Currently, the library provides three mapping functions: RAID0,
   RAID5, HASH. The default mapping is RAID5.
-  ``cache_size``: define the size of the page cache. It uses the format
   of x(k, K, m, M, g, G). The default cache size is 512MB.
-  ``num_nodes``: defines the number of NUMA nodes where the page cache
   is allocated. The default number is 1.
-  ``virt_aio``: enable virtual AIO for debugging and performance
   evaluation. It is mainly used for testing when SSDs are not
   available. By default, it’s not enabled.
