#include <iostream>
#include "FG_vector.h"
#include <map>

using namespace fg;

int main(int argc, char* argv[]) {
    std::string fn = "TEST_OUT.bin";
    FG_vector<double>::ptr test_out = FG_vector<double>::create(10);
    test_out->init_rand();
    
    std::cout << "Printing test vector:\n";
    test_out->print();
    std::cout << "Writing out ...\n";
    test_out->to_file(fn);

    try {
        std::cout << "Reading in ...\n";
        FG_vector<double>::ptr test_in = FG_vector<double>::from_file(fn);
        std::cout << "Printing test vector read in:\n";
        test_in->print();

        BOOST_VERIFY(*test_out == *test_in);
        std::cout << "Success Comparing ==\n";

        *test_out += *test_out;
        std::cout << "Printing Addition:\n";
        test_out->print();

        *test_out -= *test_in;
        std::cout << "Printing subtraction:\n";
        test_out->print();

        BOOST_VERIFY(*test_out == *test_in);
        std::cout << "Success Comparing\n";

        // *************************************** //
        std::map<size_t, double> change_map;
        std::vector<size_t> change_idx;

        change_map[0] = 2.0; change_map[3] = 23.3; change_map[9]= .3;

        for (std::map<size_t, double>::iterator it = change_map.begin();
                it != change_map.end(); ++it) {
            test_out->set(it->first, it->second);
            change_idx.push_back(it->first);
        }

        BOOST_VERIFY(!(*test_out == *test_in));

        std::cout << "nequal test ...";
        std::unique_ptr<std::vector<size_t> > neq = test_out->where_nequal(test_in);
        BOOST_VERIFY(std::equal(change_idx.begin(), change_idx.end(), neq->begin()));
        std::cout << " SUCCESS! ...\n";

        std::vector<size_t>::iterator it = neq->begin();
        for (; it != neq->end(); ++it) {
            std::cout << *it << " ";
        }
        std::cout << "\n";

    } catch (int e) {
        fprintf(stderr, "[ERROR]: Exception caught!\n");
    }

    if(remove(fn.c_str()))
        perror("Error deleting file");
    else
        puts("File successfully deleted");
    return EXIT_SUCCESS;
}
