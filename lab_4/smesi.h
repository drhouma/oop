#pragma once
#include "standart_distribution.h"
#include "interface.h"


template<class Distribution1, class Distribution2>
class MixtureDistribution : public DInterface, public PInterface {
    private:
    Distribution1 _d1;
    Distribution2 _d2;
    double _p{0};

    public:
    // must be only one constuctor with 2 standart distributions as arguments
    MixtureDistribution() = delete;
    
    MixtureDistribution(Distribution1 &d1, Distribution2 &d2, double p) : _d1(d1), _d2(d2), _p(p) {
        if (p < 0 || p > 1) throw std::invalid_argument("p must be in range [0:1]");
    };

    
    Distribution1& GetFirstFunction() {return _d1;}
    
    Distribution2& GetSecondFunction() {return _d2;}
    double GetP() {return _p;}


   
    function<double(double)> Density() override {
        std::function<double(double)> func1 = _d1.Density();
        std::function<double(double)> func2 = _d2.Density();
        return [this, func1, func2](double x) {
            return (1 - _p) * (func1(x)) + _p * func2(x); 
      };
    }

    
    double Expectation() override {
        return (1-_p) * _d1.Expectation() + _p *_d2.Expectation();
    }

    
    double Variance() override {
        double variance_1 =_d1.Variance();
        double variance_2 = _d2.Variance();
        double me_1 = _d1.Expectation();
        double me_2 = _d2.Expectation();
        double me_mixture = Expectation();

        return (1 - _p) * (me_1 * me_1 + variance_1) + _p * (me_2 * me_2 + variance_2) -
               me_mixture * me_mixture;
    }

    
    double Asymmetry() override {
        double variance_1 =_d1.Variance();
        double variance_2 = _d2.Variance();
        double me_1 = _d1.Expectation();
        double me_2 = _d2.Expectation();
        double me_mixture = Expectation();
        double mixture_var = Variance();
        return 1 / pow(mixture_var, 3.0/2) * ((1 - _p)* (pow(me_1 - me_mixture, 3) + 3 * (me_1 - me_mixture)* variance_1 )+ 
                                            _p * (pow(me_2 - me_mixture, 3) + 3 * (me_2 - me_mixture)* variance_2));
    }

    
    
    double Excess() override {
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
    
    double GenerateValue() override {
        double r = nstu::random();
        if (r >= _p) {
            return _d1.GenerateValue();
        }
        else {
            return _d2.GenerateValue();
        }
    }

    ~MixtureDistribution() {}

    void Save(ofstream& file) override {
        if(!file.is_open()){
            throw invalid_argument("File didn't open.");
        }
        file << this->Density()(0) << endl;
        file << this->Expectation() << endl;
        file << this->Variance() << endl;
        file << this->Asymmetry() << endl;
        file << this->Excess() << endl;

        cout << "Data for mixture distribution was saved." << endl;
    }
    void Load(ifstream& file, vector<double> &options) override {
        if (options.size() < 5){
            throw invalid_argument("Vector for options must have a minimum size of 5 (Density, Expectation, Variance, Asymmetry, Excess).");
        }
        if(!file.is_open()){
            throw invalid_argument("File didn't open.");
        }

        file >> options[0];
        file >> options[1];
        file >> options[2];
        file >> options[3];
        file >> options[4];
        cout << "Data for mixture distribution was loaded." << endl;
    }
};

