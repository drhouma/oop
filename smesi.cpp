#include "smesi.h"

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

    return (1 / (pow(mixture_var, 3/2))) *
           (
                   (1 - p)* (pow(me_1 - me_mixture, 3) +3 * (me_1 - me_mixture)* variance_1 + pow(variance_1, 3/2) * cosine_power_mathematical_excess(f1.v, f1.mu, f1.lambda))
                   +    p * (pow(me_2 - me_mixture, 3) +3 * (me_2 - me_mixture)* variance_2 + pow(variance_2, 3/2) * cosine_power_mathematical_excess(f2.v, f2.mu, f2.lambda))
           );
}


double mixture_mathematical_asymmetry(distribution f1, distribution f2, double p) {
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
    double mixture_mat_excess = mixture_mathematical_excess(f1, f2, p);
    double m_ass_1 = 0, m_ass_2 = 0;

    return (1 / (mixture_var * mixture_var)) * (
            (1 - p) * (pow(me_1 - me_mixture, 4) + 6 * pow(me_1 - me_mixture, 2) * variance_1 + 4 * (me_1 - me_mixture) * pow(variance_1, 3/2) * mat_excess_1 + pow(variance_1, 2) * m_ass_1 )
            +
            p * (pow(me_2 - me_mixture, 4) + 6 * pow(me_2 - me_mixture, 2) * variance_2 + 4 * (me_2 - me_mixture) * pow(variance_2, 3/2) * mat_excess_2 + pow(variance_2, 2) * m_ass_2 )
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