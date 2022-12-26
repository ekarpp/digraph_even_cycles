/* Copyright 2022 Eetu Karppinen. Subject to the MIT license. */
#ifndef TEST_H
#define TEST_H

#include <iostream>

class Test
{
private:
    void run();
public:
    int tests = 10000;
    bool end_test(int err)
    {
        if (err)
            std::cout << err << " error(s)";
        else
            std::cout << "CLEAR";

        std::cout << std::endl;

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

        std::cout << dash << std::endl;
        std::cout << whole << std::endl;
        std::cout << dash << std::endl;
    }
};

#endif
