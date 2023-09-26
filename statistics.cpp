#include "statistics.h"

int sturges_rule(int n) { return 1 + ceil(log2(n)); }

double random() { return rand() / double(RAND_MAX); }

// Theoretical distribution

std::function<double(double)> cosine_power_density(double v) {
    return [v](double x) {
        if ((-1 <= x) && (x <= 1)) {
            return sqrt(acos(-1.0)) * tgamma(v / 2.0 + 1.0) *
                   pow(cos(acos(-1.0) * x / 2.0), v) / 2.0 / tgamma((v + 1.0) / 2.0);
        } else {
            return 0.0;
        }
    };
}

std::function<double(double)> shift_scaled_density(std::function<double(double)> func, double mu, double lambda) {
    if (lambda == 0) {
        throw std::invalid_argument("Parameter lambda can't be equal zero");
    }
    return
            [func, mu, lambda](double x) { return func((x - mu) / lambda) / lambda; };
}

double cosine_power_mathematical_expectation(double v, double mu, double lambda) {
    return 0.0 + mu;
}

double cosine_power_mathematical_variance(double v, double mu, double lambda) {
    if (lambda == 0) {
        throw std::invalid_argument("Parameter lambda can't be equal zero");
    }
    return 2.0 * polygamma((v + 2) / 2.0, 1) / pow(acos(-1.0), 2) * lambda * lambda;
}

double cosine_power_mathematical_asymmetry(double v, double mu, double lambda) {
    return 0.0;
}

double cosine_power_mathematical_excess(double v, double mu, double lambda) {
    return -polygamma((v + 2) / 2.0, 3) / 2.0 /
           pow(polygamma((v + 2) / 2.0, 1), 2);
}

// Emperical distribution

double generate_cosine_power_value(double v) {
    double r1 = random(), r2 = random();
    return 2.0 * asin(sqrt(1 - pow(r1, 2.0 / v)) * cos(2.0 * acos(-1.0) * r2)) /
           acos(-1.0);
}

double generate_shift_scaled_cosine_power_value(double v, double mu, double lambda) {
    return mu+lambda*generate_cosine_power_value(v);
}

std::map<std::pair<double, double>, double> get_emperical_density(std::vector<double>& data) {
    if (data.size() == 0) {
        throw std::invalid_argument("Data can't be empty");
    }
    int k = sturges_rule(data.size());
    double min_elem = *(std::min_element(data.begin(), data.end())), max_elem = *(std::max_element(data.begin(), data.end()));
    std::map<std::pair<double, double>, double> density;
    double h = (max_elem-min_elem) / k;
    for (int i = 1; i <= k; ++i) {
        std::pair<double, double> key = std::make_pair(min_elem + (i - 1) * h, min_elem + i * h);
        density[key] = 0;
    }
    for (auto elem : data) {
        for (int i = 1; i <= k; ++i) {
            double left = min_elem + (i - 1) * h, right = min_elem + i * h;
            if ((i != k) && (left <= elem) && (elem < right)) {
                ++density[std::make_pair(left, right)];
                break;
            } else if ((i == k) && (left <= elem) && (elem <= right)) {
                ++density[std::make_pair(left, right)];
                break;
            }
        }
    }
    for (auto& key : density) {
        density[key.first] /= double(h * data.size());
    }
    return density;
}

double emperical_moment(std::vector<double>& data, int p, bool central) {
    if (size(data) == 0) {
        throw std::length_error("Vector size must be higher then 0");
    }
    double mean = (!central) ? 0 : emperical_moment(data, 1, false);
    double output = 0;
    for (auto elem : data) {
        output += pow(elem - mean, p);
    }
    return output / size(data);
}

double emperical_mathematical_expectation(std::vector<double>& data) {
    return emperical_moment(data, 1, false);
}

double emperical_variance(std::vector<double>& data) {
    return emperical_moment(data, 2, true);
}

double emperical_asymmetry(std::vector<double>& data) {
    return emperical_moment(data, 3, true) / pow(emperical_moment(data, 2, true), 1.5);
}

double emperical_excess(std::vector<double>& data) {
    return emperical_moment(data, 4, true) / pow(emperical_moment(data, 2, true), 2) - 3;
}

double generate_from_emperical_distribution(std::map<std::pair<double, double>, double> func) {
    double val = random(), proba = 0;
    for (auto& segment : func) {
        proba += segment.second*(segment.first.second - segment.first.first);
        if (val <= proba) {
            double output = random() * (segment.first.second - segment.first.first) + segment.first.first;
            return output;
        }
    }
}

double get_value_of_emperical_density(std::map<std::pair<double, double>, double>& func, double x) {
    for (auto& segment : func) {
        if ((segment.first.first <= x) && (x <= segment.first.second))
            return segment.second;
    }
    return 0.0;
}