Run FlashX in the Amazon cloud
==============================

Launch an instance with an FlashX AMI
-------------------------------------

When launching an Amazon instance, the first step is to choose an Amazon
Machine Image (AMI). Instead of choosing an AMI in "Quick Start" tab,
users can choose one from "Community AMIs". Users can search for
"FlashX". Currently, there is an AMI named "FlashX-v1" in the search
result. After choosing the wanted AMI, users can launch an Amazon
instance. Users can select an instance type based on their needs. To run
FlashX on SSDs, users can choose ``i2`` instances. The largest ``i2``
instance is ``i2.8xlarge``, which has 8 SSDs.

Configure FlashX in an Amazon instance
--------------------------------------

The FlashX AMI has preinstalled all of the FlashX components in
``/home/ubuntu/FlashX``. FlashX is compiled in
``/home/ubuntu/FlashX/build`` and FlashR has been preinstalled.

| If users want to use FlashX with SSDs, users need to configure SSDs
  for FlashX. To simplify the configuration, we provide a perl script
  ``conf/set_amazon_ssds.pl``. It takes a device file, where each line
  is an SSD device. An example of the device file is shown in
  ``conf/ssds.txt``. Users only need to provide device names such as
  ``xvd?``. When users run
  ``perl conf/set_amazon_ssds.pl conf/ssds.txt`` in the top directory of
  FlashX, it configures the system parameters and outputs a file named
  ``conf/data_files.txt``. **NOTE: ``set_amazon_ssds.pl`` formats the
  SSD devices by default**. Please use it with caution.
| After ``data_files.txt`` is generated, users need to set the parameter
  ``root_conf`` in the configuration file to specify the location where
  SSDs are mounted. In the default configuration file
  ``conf/run_amazon.txt``, ``root_conf`` points to
  ``/home/ubuntu/FlashX/conf/data_files.txt``, which is the location
  where the ``conf/set_amazon_ssds.pl`` script generates
  ``data_files.txt``.
