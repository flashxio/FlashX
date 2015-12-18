#include <iostream>
#include "FG_vector.h"

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
        std::cout << "Success Comparing ==!\n";

        *test_out += *test_out;
        std::cout << "Printing Addition:\n";
        test_out->print();

        *test_out -= *test_in;
        std::cout << "Printing subtraction:\n";
        test_out->print();

        BOOST_VERIFY(*test_out == *test_in);
        std::cout << "Success Comparing!\n";
    } catch (int e) {
        fprintf(stderr, "[ERROR]: Exception caught!\n");
    }

    if(remove(fn.c_str()))
        perror("Error deleting file");
    else
        puts("File successfully deleted");
    return EXIT_SUCCESS;
}
