#ifndef EXTENSION_TEST_H
#define EXTENSION_TEST_H

#include "test.hh"
#include "../src/extension.hh"

class Extension_test : Test
{
private:
    int n;

    void test_add_inverse();
    void test_associativity();
    void test_mul();
    void test_even_tau();
    void test_is_even();

    void run()
    {
        test_add_inverse();
        test_associativity();
        test_mul();
        test_even_tau();
        test_is_even();
    }

public:
    Extension_test(int deg);

};

#endif
