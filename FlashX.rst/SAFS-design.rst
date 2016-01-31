File layout.
============

SAFS splits a file into multiple partitions and each SSD gets a
partition of the file. SAFS deploys a Linux filesystem on the SSDs and
each partition is stored as a native file in the Linux filesystem. To
map data in a SAFS file to SSDs, SAFS first splits the file into blocks
of the user-defined size and map each block to a particular SSD and a
particular location in the native file on the SSD with a map function.
Therefore, SAFS does not need to maintain metadata to get the physical
location of a data block in an SAFS file. There are three map functions
available right now:

-  Striping: Data are divided into fixed-size small blocks placed on
   successive disks in increasing order.
-  Rotated Striping: Data are divided into stripes but
   the start disk for each stripe is rotated, much like distributed
   parity in RAID5.
-  Hash mapping: The placement of each block is randomized among all
   disks.

| The figure below shows the data layout of an SAFS file. In this
  example, the SAFS file has 9 blocks and each SSD gets 3 blocks. The 3
  blocks on an SSD are stored in a native file. This figure also
  illustrates the Striping map function. The first data block in the
  SAFS file is mapped to the first data block in the native file on the
  first SSD, and the second data block in the SAFS file is mapped to the
  first data block in the native file on the second SSD, and so on.
| |data organization|.

The details of SAFS design is described in the
`paper <http://dl.acm.org/citation.cfm?id=2503225>`__ published in
SuperComputing'13.

.. |data organization| image:: http://www.cs.jhu.edu/~zhengda/libsafs/SAFS_file_layout.png
