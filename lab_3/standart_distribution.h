#include <fstream>
#include <iostream>
#include<cmath>
#include<functional>
#include<vector>
#include<stdexcept>
#include <string>

struct DistributionParameters {
    double v;
    double mu;
    double lambda;
};



using namespace std;

class CosinePower {
private:

    double v;//form coefficient
    double mu;//shift coefficient
    double lambda;//scale coefficient

public:
    //For standart distribution default: shift=0, scale=1. V - any.
    //For shift-scaled use other parameters
    CosinePower(double form = 1, double shift = 0, double scale = 1):
    v(form), mu(shift), lambda(scale) {
        if (lambda == 0) {
            throw std::invalid_argument("Parameter lambda can't be equal zero");
        }
    }

    CosinePower(DistributionParameters p):
    v(p.v), mu(p.mu), lambda(p.lambda) {
        if (lambda == 0) {
            throw std::invalid_argument("Parameter lambda can't be equal zero");
        }
    }
    CosinePower(ifstream& file);

    void SetForm(double form) {v = form;}
    void SetShift(double shift) {mu = shift;}
    void SetScale(double scale);

    double GetForm() {return v;}
    double GetShift() {return mu;}
    double GetScale() {return lambda;}

    void Save(ofstream& file);
    void Load(ifstream& file, vector<double> &options);

    // Returns function taking x and retrieving value of cosine-power density
    function<double(double)> Density();

    // Returns math expectation of cosine-power distribution
    double Expectation();

    // Returns variance of cosine-power distribution with
    double Variance();

    // Returns variance of cosine-power distribution
    double Asymmetry();

    // Returns excess of cosine-power distribution
    double Excess();

    // Generate random value from shift-scaled cosine_power distribution
    double Generate_cosine_power_value();

};