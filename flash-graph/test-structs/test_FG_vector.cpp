/**
 * Copyright 2014 Open Connectome Project (http://openconnecto.me)
 * Written by Disa Mhembere (disa@jhu.edu)
 *
 * This file is part of FlashGraph.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#include <stdio.h>
#ifdef PROFILER
#include <google/profiler.h>
#endif

#define VERBOSE 1

#include "../FG_vector.h"

void test_shallow_copy() {

  FG_vector<float>::ptr fgv1 = FG_vector<float>::create(5);
  fgv1->set(1, 1); fgv1->set(2,2);

  FG_vector<float>::ptr fgv2 = FG_vector<float>::create(5);
  fgv2->shallow_copy(fgv1);
  
  fgv2->set(3, 3.6);
  
#if VERBOSE
  printf("Printing fgv1:");
  fgv1->print();

  printf("Printing fgv2:");
  fgv2->print();
#endif

  assert(fgv1->get(3) != fgv2->get(3));
  printf("Shallow copy test successful!\n");
}

int main(int argc, char* argv[]) {
  
  test_shallow_copy();
  return 0;
}
