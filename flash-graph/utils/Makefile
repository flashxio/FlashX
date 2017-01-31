# Copyright 2017 Open Connectome Project (http://openconnecto.me)
# Written by Disa Mhembere (disa@jhu.edu)
#
# This file is part of FlashGraph.
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

LDFLAGS := $(LDFLAGS)
CXXFLAGS += -I../../../../libcommon -I../../libsafs -I.. -I.

FILES = convert_matrix create_matrix

all: dirs $(FILES) test

mkdir = mkdir
utils_dir = ../../../../utils/

.PHONY: dirs

dirs: ${utils_dir}

${utils_dir}:
	${mkdir} ${utils_dir}

convert_matrix: convert.o convert_matrix.o
	$(CXX) -o ${utils_dir}convert_matrix convert.o convert_matrix.o $(LDFLAGS)

create_matrix: create_matrix.o
	$(CXX) -o ${utils_dir}create_matrix create_matrix.o $(LDFLAGS)

test: convert_matrix
	$(MAKE) -C test

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f $(FILES)
	$(MAKE) clean -C test

-include $(DEPS)
