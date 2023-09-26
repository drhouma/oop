#pragma once

#include "statistics.h"

struct distribution {
    double v;
    double mu;
    double lambda;
};

// возвращает функцию плотности смеси распределений f1 = {v, mu, lambda}, f2 = {v, mu, lambda} - параметры v, смещения и растяжения (0 <= p <= 1)
std::function<double(double)> mixture_density(distribution f1, distribution f2, double p);

// вычисляет матожидание смеси распределений
double mixture_mathematical_expectation(distribution f1, distribution f2, double p);

// вычисляет дисперсию смеси распределений 0 <= p <= 1
double mixture_variance(distribution f1, distribution f2, double p);

// вычисляет эксцесс смеси распределения 0 <= p <= 1
double mixture_mathematical_excess(distribution f1, distribution f2, double p);

// вычисляет коэффициент ассиметрии распределения смаси распределений 0 <= p <= 1
double mixture_mathematical_asymmetry(distribution f1, distribution f2, double p);

// генерирует случайную величину распределенную из смеси
double mixture_generate_value(distribution f1, distribution f2, double p);