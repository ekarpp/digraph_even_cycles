/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef TEST_H
#define TEST_H

#include <iostream>

class Test
{
private:


public:
    int tests = 10000;

    Test(const int t = 0)
    {
        if (t)
            this->tests = t;
    }

    bool end_test(int err)
    {
        if (err)
        {
            std::cout << "\033[31m";
            std::cout << err << " error(s)";
        }
        else
        {
            std::cout << "\033[32m";
            std::cout << "CLEAR";
        }
        std::cout << "\033[0m" << std::endl;

        /* did tests fail? */
        return (err != 0);
    }

    void start_tests(const std::string &name)
    {
        std::string whole = "testing ";
        whole.append(name);
        for (auto & c: whole)
            c = std::toupper(c);

        std::string dash(whole.length(), '-');

        std::cout << "\033[36m" << dash << std::endl;
        std::cout << whole <<  std::endl;
        std::cout << dash << "\033[0m" << std::endl;
    }
};

#endif
