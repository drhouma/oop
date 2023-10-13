#include "smesi.h"
#include <fstream>


// плотность, дисперсия, к-ты эксцесса и ассиметрии,

// 0 <= p <= 1

std::function<double(double)> mixture_density(distribution f1, distribution f2,
                                              double p) {
    if (abs(p) > 1)
        throw std::invalid_argument("invalid fraction parameter p (0 <= p <= 1)");

    return [f1, f2, p](double x) {
        std::function<double(double)> func1 =
                shift_scaled_density(cosine_power_density(f1.v), f1.mu, f1.lambda);
        std::function<double(double)> func2 = shift_scaled_density(cosine_power_density(f2.v), f2.mu, f2.lambda);
        return (1 - p) * func1(x) + p * func2(x);
    };
}

double mixture_mathematical_expectation(distribution f1, distribution f2, double p) {
    double me_1 = cosine_power_mathematical_expectation(f1.v, f1.mu, f1.lambda);
    double me_2 = cosine_power_mathematical_expectation(f2.v, f2.mu, f2.lambda);
    return (1 - p) * me_1 + p * me_2;
}

double mixture_variance(distribution f1, distribution f2, double p) {
    if (abs(p) > 1)
        throw std::invalid_argument("invalid fraction parameter p (0 <= p <= 1)");
    double variance_1 =
            cosine_power_mathematical_variance(f1.v, f1.mu, f1.lambda);
    double variance_2 =
            cosine_power_mathematical_variance(f2.v, f2.mu, f2.lambda);
    double me_1 = cosine_power_mathematical_expectation(f1.v, f1.mu, f1.lambda);
    double me_2 = cosine_power_mathematical_expectation(f2.v, f2.mu, f2.lambda);
    double me_mixture = mixture_mathematical_expectation(f1, f2, p);

    return (1 - p) * (me_1 * me_1 + variance_1) + p * (me_2 * me_2 + variance_2) -
           me_mixture * me_mixture;
}

// gamma 1
double mixture_mathematical_asymmetry(distribution f1, distribution f2, double p) {
    if (abs(p) > 1)
        throw std::invalid_argument("invalid distribution parameter p (0 <= p <= 1)");

    double variance_1 =
            cosine_power_mathematical_variance(f1.v, f1.mu, f1.lambda);
    double variance_2 =
            cosine_power_mathematical_variance(f2.v, f2.mu, f2.lambda);
    double me_1 = cosine_power_mathematical_expectation(f1.v, f1.mu, f1.lambda);
    double me_2 = cosine_power_mathematical_expectation(f2.v, f2.mu, f2.lambda);
    double me_mixture = mixture_mathematical_expectation(f1, f2, p);
    double mixture_var = mixture_variance(f1, f2, p);
    return 1 / pow(mixture_var, 3.0/2) * ((1 - p)* (pow(me_1 - me_mixture, 3) + 3 * (me_1 - me_mixture)* variance_1 )+ p * (pow(me_2 - me_mixture, 3) + 3 * (me_2 - me_mixture)* variance_2));
    return 1 / pow(mixture_var, 1.5) * ((1 - p) * (pow(me_1 - me_mixture, 3) + 3 * (me_1 - me_mixture)*variance_1) + p * (pow(me_2 - me_mixture, 3) + 3 * (me_2 - me_mixture)*variance_2));
}

// gamma 2
double mixture_mathematical_excess(distribution f1, distribution f2, double p) {
    if (abs(p) > 1)
        throw std::invalid_argument("invalid fraction parameter p (0 <= p <= 1)");

    double variance_1 =
        cosine_power_mathematical_variance(f1.v, f1.mu, f1.lambda);
    double variance_2 =
        cosine_power_mathematical_variance(f2.v, f2.mu, f2.lambda);
    double me_1 = cosine_power_mathematical_expectation(f1.v, f1.mu, f1.lambda);
    double me_2 = cosine_power_mathematical_expectation(f2.v, f2.mu, f2.lambda);
    double me_mixture = mixture_mathematical_expectation(f1, f2, p);
    double mixture_var = mixture_variance(f1, f2, p);
    double mat_excess_1 = cosine_power_mathematical_excess(f1.mu, f1.mu, f1.lambda);
    double mat_excess_2 = cosine_power_mathematical_excess(f2.mu, f2.mu, f2.lambda);
    double m_as_1 = cosine_power_mathematical_asymmetry(f1.v, f1.mu, f1.lambda), m_as_2 = cosine_power_mathematical_asymmetry(f1.v, f1.mu, f1.lambda);

    return (1 / (mixture_var * mixture_var)) * (
        (1 - p) * (pow(me_1 - me_mixture, 4) + 6 * pow(me_1 - me_mixture, 2) * variance_1 + 4 * (me_1 - me_mixture) * pow(variance_1, 3 / 2) * m_as_1 + pow(variance_1, 2) * (mat_excess_1 + 3))
        +
        p * (pow(me_2 - me_mixture, 4) + 6 * pow(me_2 - me_mixture, 2) * variance_2 + 4 * (me_2 - me_mixture) * pow(variance_2, 3 / 2) * m_as_2 + pow(variance_2, 2) * (mat_excess_2 + 3))
        ) - 3;

}




double mixture_generate_value(distribution f1, distribution f2, double p) {
    if (p < 0 || p > 1) {
        throw std::invalid_argument("Parameter p must be in interval [0, 1]");
    }
    double r = random();
    if (r >= p) {
        return generate_shift_scaled_cosine_power_value(f1.v, f1.mu, f1.lambda);
    }
    else {
        return generate_shift_scaled_cosine_power_value(f2.v, f2.mu, f2.lambda);
    }
}

void mixture_eval_theor_and_emperical_chars(int n, distribution f1, distribution f2, double p) {
    // Записываем теоретическую плотность
    auto func = mixture_density(f1, f2, p);
    // Генерируем случайные x
    std::vector<double> data(n);
    for (int i = 0; i < n; ++i)
        data[i] = mixture_generate_value(f1, f2, p);
   std::ofstream file1("data/data.txt");
    // generate mathematical expectation (theorethical and emperical)
   file1 << n << std::endl;
    double me_t = mixture_mathematical_expectation(f1, f2, p), me_e = emperical_mathematical_expectation(data);
    file1 << me_t << ' ' << me_e << std::endl;

    // generate variance (theorethical and emperical)
    double var_t = mixture_variance(f1, f2, p), var_e = emperical_variance(data);
    file1<< var_t << ' ' << var_e << std::endl;

    // generate asymmetry coef (theorethical and emperical)
    double as_t = mixture_mathematical_asymmetry(f1, f2, p), as_e = emperical_asymmetry(data);
    file1 << as_t << ' ' << as_e << std::endl;

    // generate excess coef (theorethical and emperical)
    double excess_t = mixture_mathematical_excess(f1, f2, p), excess_e = emperical_excess(data);
    file1 << excess_t << ' ' << excess_e << std::endl;

    //file1.close();
}