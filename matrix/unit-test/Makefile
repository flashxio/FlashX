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

include ../../Makefile.common

OMP_FLAG = -fopenmp
CXXFLAGS += -I../../include -I../../libcommon -I../../flash-graph -I.. $(OMP_FLAG)
LDFLAGS += $(OMP_FLAG) -L../ -lmatrix -L../../flash-graph -lgraph

all: test_io_gen

test_io_gen: test_io_gen.o
	$(CXX) -o test_io_gen test_io_gen.o $(LDFLAGS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm test_io_gen

-include $(DEPS) 
