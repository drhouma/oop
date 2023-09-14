#pragma once

#include <cmath>
#include <stdexcept>
#include "statistics.h"
// плотность, дисперсия, к-ты эксцесса и ассиметрии
struct distribution {
    double v;
    double mu;
    double lambda;
};

//Плотность смеси
// 0 <= p <= 1
std::function<double(double)> mixture_density(distribution f1, distribution f2, double p);

// вычисляет матожидание смеси распределений
// пока вместо мат ожидания соответствующих распределений будет их сдвиг
double mixture_mathematical_expectation(distribution f1, distribution f2, double p);

// вычисляет дисперсию смеси распределений 0 <= p <= 1
double mixture_variance(distribution f1, distribution f2, double p);

// вычисляет эксцесс смеси распределений
double mixture_mathematical_excess(distribution f1, distribution f2, double p);

double mixture_mathematical_asymmetry(distribution f1, distribution f2, double p);

double mixture_generate_value(distribution f1, distribution f2, double p);