# Copyright 2014 Open Connectome Project (http;//openconnecto.me)
# Written by Da Zheng (zhengda1936@gmail.com)
#
# This file is part of SAFSlib.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

include Makefile.common

all: build_lib unit_test tools flash-graph test utils

build_lib:
	$(MAKE) -C libcommon
	$(MAKE) -C libsafs

unit_test: build_lib
ifndef MEMCHECK
	$(MAKE) -C unit-test
endif

test: build_lib
	$(MAKE) -C test

tools: build_lib
	$(MAKE) -C tools

utils: build_lib
	$(MAKE) -C utils

flash-graph: build_lib
	$(MAKE) -C flash-graph

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f include/*~
	find -name core -delete
	make --ignore-errors -C unit-test clean
	make --ignore-errors -C test clean
	make --ignore-errors -C libsafs clean
	make --ignore-errors -C libcommon clean
	make --ignore-errors -C tools clean
	make --ignore-errors -C utils clean
	make --ignore-errors -C flash-graph clean

-include $(DEPS) 
