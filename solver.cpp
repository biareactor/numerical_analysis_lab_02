#include "solver.h"

double Solver::integrate(func function, double from, double to)
{
    return function( (from + to) / 2. ) * (to - from);
}

vecvec Solver::calc_coefs()
{
    const double UNUSED_VALUE = 0;
    vecvec res = { {UNUSED_VALUE}, {UNUSED_VALUE}, {UNUSED_VALUE} };

    double x_curr = step;
    for (size_t i = 1; i < n + 1; i++)
    {
        if (x_curr < break_point)
        {
            double a = 1. / (1. / step * integrate(k_1, x_curr - step, x_curr));
            res[0].push_back(a);
        }
        else if (x_curr - break_point < step)
        {
            double a = 1. / (1. / step * integrate(k_1, x_curr - step, break_point)
                           + 1. / step * integrate(k_2, break_point, x_curr));
            res[0].push_back(a);
        }
        else
        {
            double a = 1. / (1. / step * integrate(k_2, x_curr - step, x_curr));
            res[0].push_back(a);
        }

        x_curr += step;
    }

    x_curr = 3 * step / 2.;

    for (size_t i = 1; i < n; i++)
    {
        if (x_curr < break_point)
        {
            double fi = 1. / step * integrate(f_1, x_curr - step, x_curr);
            double d = 1. / step * integrate(q_1, x_curr - step, x_curr);

            res[1].push_back(d);
            res[2].push_back(fi);
        }
        else if (x_curr - break_point < step)
        {
            double fi = 1. / step * integrate(f_1, x_curr - step, break_point)
                     + 1. / step * integrate(f_2, break_point, x_curr);
            double d = 1. / step * integrate(q_1, x_curr - step, break_point)
                     + 1. / step * integrate(q_2, break_point, x_curr);

            res[1].push_back(d);
            res[2].push_back(fi);
        }
        else
        {
            double fi = 1. / step * integrate(f_2, x_curr - step, x_curr);
            double d = 1. / step * integrate(q_2, x_curr - step, x_curr);

            res[1].push_back(d);
            res[2].push_back(fi);
        }

        x_curr += step;
    }

    return res;
}

vecvec Solver::coefs_to_system(vecvec coefs)
{
    vecvec res = { {0}, {0}, {0} , coefs[2]};

    for (size_t i = 1; i < n; i++)
    {
        double a_i  = coefs[0][i];
        double a_i1 = coefs[0][i + 1];
        double d_i  = coefs[1][i];

        double a = a_i  / (step * step);
        double b = a_i1 / (step * step);
        double c = a + b + d_i;

        res[0].push_back(a);
        res[1].push_back(c);
        res[2].push_back(b);
    }

    return res;
}

vec Solver::solve_matrix(vecvec matrix)
{
    vec a  = matrix[0];
    vec c  = matrix[1];
    vec b  = matrix[2];
    vec fi = matrix[3];

    vec y(n+1);
    vec alpha(n+1);
    vec beta(n+1);

    alpha[1] = 0;
    beta[1]  = mu_1;

    for (size_t i = 2; i < n+1; i++)
    {
        double denom = c[i-1] - alpha[i-1] * a[i-1];
        alpha[i] = b[i-1] / denom;
        beta[i]  = (fi[i-1] + beta[i-1] * a[i-1]) / denom;
    }

    y[n] = mu_2;

    for(int i = n-1; i >= 0; i--)
    {
        y[i] = alpha[i+1] * y[i+1] + beta[i+1];
    }

    return y;
}

Solver::Solver(
    int n_,
    double break_point_,
    func k_1_, func q_1_, func f_1_,
    func k_2_, func q_2_, func f_2_,
    double mu_1_, double mu_2_
)
{
    n = n_;
    step = 1. / n_;
    break_point = break_point_;
    k_1 = k_1_; q_1 = q_1_; f_1 = f_1_;
    k_2 = k_2_; q_2 = q_2_; f_2 = f_2_;
    mu_1 = mu_1_; mu_2 = mu_2_;

}

vec Solver::solve()
{
    return(solve_matrix(coefs_to_system(calc_coefs())));
}

