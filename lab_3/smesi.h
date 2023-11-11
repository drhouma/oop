#pragma once
#include "standart_distribution.h"


class MixtureDistribution {
    private:
    CosinePower _d1;
    CosinePower _d2;
    double _p{0};

    public:
    // must be only one constuctor with 2 standart distributions as arguments
    MixtureDistribution() = delete;
    MixtureDistribution(CosinePower &d1, CosinePower &d2, double p) : _d1(d1), _d2(d2), _p(p) {
        if (p < 0 || p > 1) throw std::invalid_argument("p must be in range [0:1]");
    };


    private:
    CosinePower& GetFirstFunction() {return _d1;}
    CosinePower& GetSecondFunction() {return _d2;}

    public:
    function<double(double)> Density();

    // Returns math expectation of cosine-power distribution
    double Expectation();

    // Returns variance of cosine-power distribution with
    double Variance();

    // Returns variance of cosine-power distribution
    double Asymmetry();

    // Returns excess of cosine-power distribution
    double Excess();
    
    ~MixtureDistribution() {}
};