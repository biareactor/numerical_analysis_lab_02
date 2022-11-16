#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>

namespace functions
{
    double k_1_test(double x)
    {
        return 1. / 1.4;
    }

    double k_2_test(double x)
    {
        return 1. / 0.4;
    }

    double q_1_test(double x)
    {
        return 0.4;
    }

    double q_2_test(double x)
    {
        return 0.16;
    }

    double f_1_test(double x)
    {
        return 0.4;
    }

    double f_2_test(double x)
    {
        return exp(-0.4);
    }

    double k_1_main(double x)
    {
        return 1. / (x + 1);
    }

    double k_2_main(double x)
    {
        return 1. / x;
    }

    double q_1_main(double x)
    {
        return x;
    }

    double q_2_main(double x)
    {
        return x*x;
    }

    double f_1_main(double x)
    {
        return x;
    }

    double f_2_main(double x)
    {
        return exp(-x);
    }

    double u_1_test(double x)
    {
        return 0.060557222866651 * exp(sqrt(2./7.)*x) -1.060557222866651 * exp(-sqrt(2./7.)*x) + 1;
    }

    double u_2_test(double x)
    {
        return -0.472024550734437 * exp(sqrt(2./5.)*x) -4.331084823580059 * exp(-sqrt(2./5.)*x) + 25./(4*exp(2./5.));
    }
}


#endif // FUNCTIONS_H
