#include "empiric.h"

double Empiric::GenerateValue() {
    double val = nstu::random(), proba = 0;
    for (auto& segment : fr) {
        proba += segment.second*(segment.first.second - segment.first.first);
        if (val <= proba) {
            double output = nstu::random() * (segment.first.second - segment.first.first) + segment.first.first;
            return output;
        }
    }
}

Empiric::Empiric(int n0, CosinePower& prim, int k0) :
    n(n0>1 ? n0 : throw std::invalid_argument("numbers quantity must be > 1")), k(k0>1 ? k0 : int(log(double(n))/log(double(2.)))+1.),
    data(new double[n]) {
    if (k < 1) throw std::invalid_argument("k must be >= 1");
    for (int i = 0; i < n; ++i) data[i] = prim.GenerateValue();
    fr = this->GetEmpericalDensity();
}

Empiric::Empiric(int n0, MixtureDistribution& mixt, int k0) :
    n(n0>1 ? n0 : throw std::invalid_argument("numbers quantity must be > 1")), k(k0>1 ? k0 : int(log(double(n))/log(double(2.)))+1.),
    data(new double[n]) {
        if (k < 1) throw std::invalid_argument("k must be >= 1");
    for (int i = 0; i < n; ++i)
        data[i] = mixt.GenerateValue();
    fr = GetEmpericalDensity();
}

/// @brief 
/// @param n0 - кол-во чисел 
/// @param prim - Эмпирическое распределение, основанное на стандартном. 
/// @param k0 - кол-во интервалов
Empiric::Empiric(int n0, Empiric& prim, int k0) :
    n(n0>1 ? n0 : throw std::invalid_argument("numbers quantity must be > 1")), k(k0>1 ? k0 : int(log(double(n))/log(double(2.)))+1.),
    data(new double[n]){
        if (k < 1) throw std::invalid_argument("k must be >= 1");
    for (int i = 0; i < n; i++) {
        data[i] = prim.GenerateValue();
    }
    fr = GetEmpericalDensity();
}

Empiric& Empiric::operator=(const Empiric& other){
    // проверка на присваивание самому себе
    if(this==&other) return *this;
    // если размеры изменяются, то перераспределяем память
    if(n!=other.n){
        delete [] data;
        n=other.n;
        data=new double[n];
    }

    for (int i = 0; i < n; i++) {
        data[i] = other.data[i];
    }
    fr = other.fr;

    return *this;
}

// no need checks because emp is already checked
Empiric::Empiric(const Empiric& other) :
    n(other.n), k(other.k), data(new double[n]) {
        for (int i = 0; i < n; i++) {
            data[i] = other.data[i];
        }
        fr = other.fr;
}


double Empiric::emperical_moment(int p, bool central) {
    double mean = (!central) ? 0 : emperical_moment(1, false);
    double output = 0;
    for (int i = 0; i < n; i++) {
        output += pow(data[i] - mean, p);
    }
    return output / n;
}

function<double(double)> Empiric::Density() {

    return [&, this](double x) {
        auto func = this->GetFr();
        for (auto& segment : func) {
        if ((segment.first.first <= x) && (x <= segment.first.second))
            return segment.second;
        }
        return 0.0;
  };
}

double Empiric::Expectation() {
    return emperical_moment(1, false);
}

double Empiric::Variance() {
    return emperical_moment(2, true);
}

double Empiric::Asymmetry() {
    return emperical_moment(3, true) / pow(emperical_moment(2, true), 1.5);
}

double Empiric::Excess() {
    return emperical_moment(4, true) / pow(emperical_moment(2, true), 2) - 3;
}

std::map<std::pair<double, double>, double> Empiric::GetEmpericalDensity() {
    double min_elem = std::numeric_limits<double>::max();
    double max_elem = std::numeric_limits<double>::min();
    for(int i=0;i<n;i++){
        if (data[i]<min_elem) min_elem = data[i];
        if (data[i] > max_elem) max_elem = data[i];
    }
   

    std::map<std::pair<double, double>, double> density;
    double h = (max_elem-min_elem) / k;
    for (int i = 1; i <= k; ++i) {
        std::pair<double, double> key = std::make_pair(min_elem + (i - 1) * h, min_elem + i * h);
        density[key] = 0;
    }
    for (int j = 0;j < n; j++) {
        for (int i = 1; i <= k; ++i) {
            double left = min_elem + (i - 1) * h, right = min_elem + i * h;
            if ((i != k) && (left <= data[j]) && (data[j] < right)) {
                ++density[std::make_pair(left, right)];
                break;
            } else if ((i == k) && (left <= data[j]) && (data[j] <= right)) {
                ++density[std::make_pair(left, right)];
                break;
            }
        }
    }
    for (auto& key : density) {
        density[key.first] /= double(h * n);
    }
    return density;
}

