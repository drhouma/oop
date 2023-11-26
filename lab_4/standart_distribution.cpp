#include "standart_distribution.h"
#include "math/smath.h"

int sturges_rule(int n) { return 1 + ceil(log2(n)); }

// don't work on linux (function with same name exists) so add namespace

function<double(double)> cosine_power_density(double v) {
    return [v](double x) {
        if ((-1 <= x) && (x <= 1)) {
            return sqrt(acos(-1.0)) * tgamma(v / 2.0 + 1.0) *
                   pow(cos(acos(-1.0) * x / 2.0), v) / 2.0 / tgamma((v + 1.0) / 2.0);
        } else {
            return 0.0;
        }
    };
}

double RandomCosinePowerValue(double v) {
    double r1 = nstu::random(), r2 = nstu::random();
    return 2.0 * asin(sqrt(1 - pow(r1, 2.0 / v)) * cos(2.0 * acos(-1.0) * r2)) /
           acos(-1.0);
}

double nstu::random() {
    return rand() / double(RAND_MAX); 
}


function<double(double)> CosinePower::Density() {
    if (mu == 0 && lambda == 1) {
        return cosine_power_density(v);
    } else {
        auto func = cosine_power_density(v);
        double shift = mu;
        double scale = lambda;
        return [func, shift, scale](double x)  { return func((x - shift) / scale) / scale; };
    }
}

void CosinePower::SetScale(double scale) {
    if (scale == 0) {
        throw std::invalid_argument("Parameter scale(lambda) can't be equal zero");
    }
    lambda = scale;
}

double CosinePower::Expectation() {
    return 0.0 + mu;
}

double CosinePower::Variance() {
    return 2.0 * polygamma((v + 2) / 2.0, 1) / pow(acos(-1.0), 2) * lambda * lambda;
}

double CosinePower::Asymmetry() {
    return 0.0;
}

double CosinePower::Excess() {
    return -polygamma((v + 2) / 2.0, 3) / 2.0 /
           pow(polygamma((v + 2) / 2.0, 1), 2);
}

double CosinePower::GenerateValue() {
    return mu+lambda*RandomCosinePowerValue(v);
}

void CosinePower::Save(std::ofstream &file) {
    if(!file.is_open()){
        throw invalid_argument("File didn't open.");
    }

    file << v << " " <<  mu << " " << lambda << endl;
    file << this->Density()(0) << endl;
    file << this->Expectation() << endl;
    file << this->Variance() << endl;
    file << this->Asymmetry() << endl;
    file << this->Excess() << endl;

    cout << "Data was saved." << endl;
}

void CosinePower::Load(std::ifstream &file, vector<double> &options) {
    if (options.size() < 5){
        throw invalid_argument("Vector for options must have a minimum size of 5 (Density, Expectation, Variance, Asymmetry, Excess).");
    }
    if(!file.is_open()){
        throw invalid_argument("File didn't open.");
    }

    file >> v  >> mu >> lambda;
    if (lambda == 0) {
        throw std::invalid_argument("Parameter scale(lambda) can't be equal zero");
    }

    file >> options[0];
    file >> options[1];
    file >> options[2];
    file >> options[3];
    file >> options[4];
    cout << "Data was loaded." << endl;

}

CosinePower::CosinePower(std::ifstream &file) {
    if(!file.is_open()){
        throw invalid_argument("File didn't open.");
    }

    file >> v  >> mu >> lambda;
    if (lambda == 0) {
        throw std::invalid_argument("Parameter scale(lambda) can't be equal zero");
    }
}