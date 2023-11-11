#include "empiric.h"


Empiric::Empiric(int n0, Primary& prim, int k0=1) :
    n(n0>1 ? n0 : throw std::invalid_argument("numbers quantity must be > 1")), k(k0>1 ? k0 : int(log(double(n))/log(double(2.)))+1.),
    data(new double[n]), fr(new double[k]){
    if (k < 1) throw std::invalid_argument("k must be >= 1");
    for (int i = 0; i < n; ++i)
        data[i] = prim.Generate_cosine_power_value();
    
}

Empiric::Empiric(int n0, MixtureDistribution& mixt, int k0=1) :
    n(n0>1 ? n0 : throw std::invalid_argument("numbers quantity must be > 1")), k(k0>1 ? k0 : int(log(double(n))/log(double(2.)))+1.),
    data(new double[n]), fr(new double[k]) {
        if (k < 1) throw std::invalid_argument("k must be >= 1");
}

Empiric::Empiric(int n0, Empiric& prim, int k0=1) :
    n(n0>1 ? n0 : throw std::invalid_argument("numbers quantity must be > 1")), k(k0>1 ? k0 : int(log(double(n))/log(double(2.)))+1.),
    data(new double[n]), fr(new double[k]){
        if (k < 1) throw std::invalid_argument("k must be >= 1");
}

Empiric::Empiric& operator=(const Empiric& emp){
    // проверка на присваивание самому себе
    if(this==&emp) return *this;
    // если размеры изменяются, то перераспределяем память
    if(n!=emp.n){
        delete [] data;
        n=emp.n;
        data=new double[n];
    }

    if(k!=emp.k){
        delete [] fr;
        k=emp.k;
        fr=new double[k];
    }
    for (int i = 0; i < n; i++) {
        data[i] = emp.data[i];
    }
    for (int i = 0; i < k; i++) {
        fr[i] = emp.fr[i];
    }
    return *this;
}

// no need checks because emp is already checked
Empiric::Empiric(const Empiric& emp) :
    n(emp.n), k(emp.k), data(new double[n]), fr(new double[k]) {
        for (int i = 0; i < n; i++) {
            data[i] = emp.data[i];
        }
        for (int i = 0; i < k; i++) {
            fr[i] = emp.fr[i];
        }
}


double Empiric::emperical_moment(int p, bool central) {
    double mean = (!central) ? 0 : emperical_moment(data, 1, false);
    double output = 0;
    for (auto elem : data) {
        output += pow(elem - mean, p);
    }
    return output / size(data);
}



double Empiric::Expectation() {
    return emperical_moment(data, 1, false);
}

double Empiric::Variance() {
    return emperical_moment(data, 2, true);
}

double Empiric::Asymmetry() {
    return emperical_moment(data, 3, true) / pow(emperical_moment(data, 2, true), 1.5);
}

double Empiric::Excess() {
    return emperical_moment(data, 4, true) / pow(emperical_moment(data, 2, true), 2) - 3;
}