/**
 * Copyright 2019
 *
 * This file is part of SAFSlib.
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

#include "FileManager.h"

int main() {
    std::string conf =
        "/home/disa/Research/graphyti/graphyti/src/conf/run_graph.txt";
    std::string fn = "tbrmd.txt";
    std::string renamefn = "tbrmd2.txt";

    // Create file manager
    printf("Create FileManager\n\n");
    fg::FileManager fm(conf);

    // Load up a file
    fm.to_ex_mem(fn, fn); // Arg1 is SAFS, Arg2 is local

    // Check size
    auto filesize = fm.file_size(fn);
    printf("\n\nINFO:\n");
    for (auto& kv : fm.info(fn))
        std::cout << kv.first.c_str() << " : " <<   kv.second.c_str() << "\n";

    printf("\n");

    assert(filesize > 0);
    assert(fm.file_exists(fn));

    // Check size
    printf("Rename %s --> %s\n", fn.c_str(), renamefn.c_str());
    fm.rename(fn, renamefn);

    assert(!fm.file_exists(fn));
    assert(fm.file_exists(renamefn));

    assert(filesize == fm.file_size(renamefn));


    printf("Deleting renamed file\n");
    fm.delete_file(renamefn);

    printf("Successful! FileManager test complete!\n");
}
