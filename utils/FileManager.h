#ifndef __FILE_MANAGER_H__
#define __FILE_MANAGER_H__
/**
 * Copyright 2019
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

#include <string>
#include <stdexcept>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unordered_map>

#include "io_interface.h"
#include "native_file.h"
#include "safs_file.h"
#include "file_mapper.h"
#include "RAID_config.h"

using namespace safs;

namespace fg {

class FileManager {
    private:
    std::string configuration_file;
    config_map::ptr configs;
    const int BUF_SIZE = 1024 * 64 * 4096;

    public:
    FileManager(const std::string& configuration_file):
        configuration_file(configuration_file) {
        configs = config_map::create(configuration_file);
    }

    void delete_file(const std::string& file_name) {
        init_io_system(configs, false);
        configs->add_options("writable=1");

        safs_file file(get_sys_RAID_conf(), file_name);

        if (!file.exist())
            throw std::runtime_error(
                    file_name + std::string(" doesn't exist\n"));

        file.delete_file();
        destroy_io_system();
    }

    void from_ex_mem(const std::string& file_name, const std::string& ext_file) {
        init_io_system(configs, false);
        FILE *f = fopen(ext_file.c_str(), "w");
        if (f == NULL) {
            fprintf(stderr, "can't open %s: %s\n", ext_file.c_str(),
                    strerror(errno));

            throw std::runtime_error(std::string("can't open ") + ext_file +
                    std::string("\n"));
        }

        file_io_factory::shared_ptr io_factory = create_io_factory(file_name,
                REMOTE_ACCESS);
        io_interface::ptr io = create_io(io_factory, thread::get_curr_thread());
        size_t phy_file_size = io_factory->get_file_size();
        size_t file_size = io_factory->get_header().get_size();
        assert(file_size <= phy_file_size);
        if (file_size == 0)
            file_size = phy_file_size;

        size_t buf_size = 16 * 1024 * 1024;
        char *buf = (char *) valloc(buf_size);
        assert(buf);
        size_t off = 0;
        while (off < file_size) {
            data_loc_t loc(io->get_file_id(), off);
            size_t read_size = buf_size;
            size_t write_size = buf_size;
            if (off + buf_size > file_size) {
                // The physical storage size of SAFS is always rounded to
                // the page size, and we access the SAFS file with direct I/O,
                // so we need to round up the read size.
                read_size = ROUNDUP(file_size - off, PAGE_SIZE);
                write_size = file_size - off;
            }
            io_request req(buf, loc, read_size, READ);
            io->access(&req, 1);
            io->wait4complete(1);

            size_t ret = fwrite(buf, write_size, 1, f);
            assert(ret == 1);
            off += read_size;
        }

        fclose(f);
        destroy_io_system();
    }

    void to_ex_mem(const std::string& int_file_name,
            const std::string& ext_file) {

        init_io_system(configs, false);
        safs_file file(get_sys_RAID_conf(), int_file_name);
        bool ret = file.load_data(ext_file);

        if (!ret)
            throw std::runtime_error(std::string("SAFS write fail with code") +
                    std::to_string(ret) + std::string("\n"));
        destroy_io_system();
    }

    void rename(const std::string& file_name,
            const std::string& new_name) {

        init_io_system(configs, false);
        safs_file file(get_sys_RAID_conf(), file_name);
        if (!file.exist()) {
            throw std::runtime_error(
                    file_name + std::string(" doesn't exist\n"));
        }

        bool ret = file.rename(new_name);
        if (!ret)
            throw std::runtime_error(std::string("Can't rename ") + file_name +
                    std::string(" to ") + new_name);
        destroy_io_system();
    }

    std::vector<std::pair<std::string, size_t> > list_files() {
        init_io_system(configs, false);
        std::set<std::string> files;
        std::vector<std::pair<std::string, size_t> > files_set;

        const RAID_config &conf = get_sys_RAID_conf();

        // First find all individual file names in the root directories.
        for (int i = 0; i < conf.get_num_disks(); i++) {
            std::string dir_name = conf.get_disk(i).get_file_name();
            native_dir dir(dir_name);
            std::vector<std::string> file_names;
            dir.read_all_files(file_names);
            files.insert(file_names.begin(), file_names.end());
        }

        for (std::set<std::string>::const_iterator it = files.begin();
                it != files.end(); it++) {
            safs_file file(conf, *it);
            if (file.exist()) {
                files_set.push_back(std::pair<std::string, size_t>(
                            file.get_name(), file.get_size()));
            } else {
                // Corrupted
                files_set.push_back(std::pair<std::string, size_t>(
                            file.get_name(), 0));
            }
        }
        destroy_io_system();
        return files_set;
    }

    const bool file_exists(const std::string& file_name) {
        init_io_system(configs, false);
        safs_file file(get_sys_RAID_conf(), file_name);
        auto ret = false;

        if (file.exist())
            ret = true;

        destroy_io_system();
        return ret;
    }

    size_t file_size(const std::string& file_name) {
        init_io_system(configs, false);
        file_io_factory::shared_ptr io_factory = create_io_factory(file_name,
                REMOTE_ACCESS);
        safs_header header = io_factory->get_header();

        auto ret = header.get_size();
        destroy_io_system();
        return ret;
    }

    std::unordered_map<std::string, std::string> info(
            const std::string& file_name) {
        init_io_system(configs, false);
        file_io_factory::shared_ptr io_factory = create_io_factory(file_name,
                REMOTE_ACCESS);
        safs_header header = io_factory->get_header();

        std::unordered_map<std::string, std::string> ret = {
            {"File", file_name},
            {"RAID block size",
                std::to_string((header.get_block_size() * PAGE_SIZE))},
            {"RAID mapping option",
                std::to_string(header.get_mapping_option())},
            {"File size", std::to_string(header.get_size())}
        };

        destroy_io_system();
        return ret;
    }

    std::string to_str() {
        return std::string("FileManager\nConfig file: ") + configuration_file
            + std::string("\n");
    }
};
} // End namespace fg
#endif
