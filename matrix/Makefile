# Copyright 2014 Open Connectome Project (http://openconnecto.me)
# Written by Da Zheng (zhengda1936@gmail.com)
#
# This file is part of FlashMatrix.
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

include ../Makefile.common

OMP_FLAG = -fopenmp
CXXFLAGS += -I../include -I../libcommon -I../flash-graph $(OMP_FLAG)
LDFLAGS += $(OMP_FLAG)

all: libmatrix test unit-test eigensolver

libmatrix: $(OBJS)
	rm -f libFMatrix.a
	ar -cvq libFMatrix.a $(OBJS)

test: libmatrix eigensolver
	$(MAKE) -C test

unit-test: libmatrix
	$(MAKE) -C unit-test

eigensolver: libmatrix
	$(MAKE) -C eigensolver

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f libFMatrix.a
	find -name core -delete
	$(MAKE) -C test clean
	$(MAKE) -C unit-test clean
	$(MAKE) -C eigensolver clean

-include $(DEPS) 
