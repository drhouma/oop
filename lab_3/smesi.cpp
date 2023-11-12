#include "smesi.h"



function<double(double)> MixtureDistribution::Density() {
    return [&, this](double x) {
        std::function<double(double)> func1 = this->GetFirstFunction().Density();
        std::function<double(double)> func2 = this->GetSecondFunction().Density();
        return (1 - _p) * (func1(x)) + _p * func2(x);
  };
}


double MixtureDistribution::Expectation() {
    return (1-_p) * _d1.Expectation() + _p *_d2.Expectation();
}

double MixtureDistribution::Variance() {
    double variance_1 =_d1.Variance();
    double variance_2 = _d2.Variance();
    double me_1 = _d1.Expectation();
    double me_2 = _d2.Expectation();
    double me_mixture = Expectation();

    return (1 - _p) * (me_1 * me_1 + variance_1) + _p * (me_2 * me_2 + variance_2) -
           me_mixture * me_mixture;
}

double MixtureDistribution::Asymmetry() {
    double variance_1 =_d1.Variance();
    double variance_2 = _d2.Variance();
    double me_1 = _d1.Expectation();
    double me_2 = _d2.Expectation();
    double me_mixture = Expectation();
    double mixture_var = Variance();
    return 1 / pow(mixture_var, 3.0/2) * ((1 - _p)* (pow(me_1 - me_mixture, 3) + 3 * (me_1 - me_mixture)* variance_1 )+ 
                                        _p * (pow(me_2 - me_mixture, 3) + 3 * (me_2 - me_mixture)* variance_2));
}

double MixtureDistribution::Excess() {
     double variance_1 =_d1.Variance();
    double variance_2 = _d2.Variance();
    double me_1 = _d1.Expectation();
    double me_2 = _d2.Expectation();
    double me_mixture = Expectation();
    double mixture_var = Variance();
    double mat_excess_1 = _d1.Excess();
    double mat_excess_2 = _d2.Excess();
    double m_as_1 = _d1.Asymmetry(), m_as_2 = _d2.Asymmetry();

    return (1 / (mixture_var * mixture_var)) * (
        (1 - _p) * (pow(me_1 - me_mixture, 4) + 6 * pow(me_1 - me_mixture, 2) * variance_1 + 4 * (me_1 - me_mixture) * pow(variance_1, 3 / 2) * m_as_1 + pow(variance_1, 2) * (mat_excess_1 + 3))
        +
        _p * (pow(me_2 - me_mixture, 4) + 6 * pow(me_2 - me_mixture, 2) * variance_2 + 4 * (me_2 - me_mixture) * pow(variance_2, 3 / 2) * m_as_2 + pow(variance_2, 2) * (mat_excess_2 + 3))
        ) - 3;

}