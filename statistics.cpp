#include "statistics.h"

int sturges_rule(int n) { return 1 + ceil(log2(n)); }

double random() { return rand() / double(RAND_MAX); }

// Theoretical distribution

std::function<double(double)> cosine_power_density(double v) {
  return [v](double x) {
    return sqrt(acos(-1.0)) * tgamma(v / 2.0 + 1.0) *
           pow(cos(acos(-1.0) * x / 2.0), v) / 2.0 / tgamma((v + 1.0) / 2.0);
  };
}

std::function<double(double)> shift_scaled_density(
    std::function<double(double)> func, double mu, double lambda) {
  if (lambda == 0) {
    throw std::invalid_argument("Parameter lambda can't be equal zero");
  }
  return
      [func, mu, lambda](double x) { return func((x - mu) / lambda) / lambda; };
}

double cosine_power_mathematical_expectation(double v, double mu,
                                             double lambda) {
  return 0.0 + mu;
}

double cosine_power_mathematical_variance(double v, double mu, double lambda) {
  if (lambda == 0) {
    throw std::invalid_argument("Parameter lambda can't be equal zero");
  }
  return 2.0 * polygamma((v + 2) / 2.0, 1) / pow(acos(-1.0), 2) * lambda *
         lambda;
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

std::map<std::pair<double, double>, double> get_emperic_density(
    std::vector<double>& data) {
  if (data.size() == 0) {
    throw std::invalid_argument("Data can't be empty");
  }
  int k = sturges_rule(data.size());
  std::map<std::pair<double, double>, double> density;
  double h = 2.0 / k;
  for (int i = 1; i <= k; ++i) {
    std::pair<double, double> key =
        std::make_pair(-1.0 + (i - 1) * h, -1.0 + i * h);
    density[key] = 0;
  }
  for (auto elem : data) {
    for (int i = 1; i <= k; ++i) {
      double left = -1.0 + (i - 1) * h, right = -1.0 + i * h;
      if ((i != k) && (left <= elem) && (elem < right)) {
        ++density[std::make_pair(left, right)];
        break;
      } else if ((i == k) && (left <= elem) && (elem <= right)) {
        ++density[std::make_pair(left, right)];
        break;
      }
    }
  }
  for (auto key : density) {
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