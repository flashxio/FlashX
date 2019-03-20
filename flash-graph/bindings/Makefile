# Copyright 2014 Open Connectome Project (http://openconnecto.me)
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

LDFLAGS := -L../libgraph-algs -lgraph-algs -L..\
	-lgraph -L../../libsafs -lsafs \
	-lrt -lz $(LDFLAGS)
CXXFLAGS += -I../../libsafs -I.. -I.

all: testit

testit: testit.o ../libgraph.a ../libgraph-algs/libgraph-algs.a
	$(CXX) -o testit testit.o $(LDFLAGS)

clean:
	rm -f *.d
	rm -f *.o
	rm -f *~
	rm -f testit

-include $(DEPS)
