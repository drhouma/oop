#pragma once

#include<cmath>
#include<functional>
#include<vector>
#include<stdexcept>
#include<map>
#include<algorithm>

#include "smath.h"

// Calculates K using sturges rule
int sturges_rule(int n);

// Generates random value in range [0, 1]
double random(void);

// Returns function taking x and retrieving value of cosine-power density with parameter equal v in that point
std::function<double(double)> cosine_power_density(double v);

// Returns shifted and scaled specified density (lambda != 0)
std::function<double(double)> shift_scaled_density(std::function<double(double)> func, double mu, double lambda);

// Returns math expectation of cosine-power distribution with paramter v, shift coefficient mu, scale coefficient lambda
double cosine_power_mathematical_expectation(double v, double mu, double lambda);

// Returns variance of cosine-power distribution with paramter v, shift coefficient mu, scale coefficient lambda
double cosine_power_mathematical_variance(double v, double mu, double lambda);

// Returns variance of cosine-power distribution with paramter v, shift coefficient mu, scale coefficient lambda
double cosine_power_mathematical_asymmetry(double v, double mu, double lambda);

// Returns excess of cosine-power distribution with paramter v, shift coefficient mu, scale coefficient lambda
double cosine_power_mathematical_excess(double v, double mu, double lambda);

// Generates random value from cosine_power distribution with parameter equal v
double generate_cosine_power_value(double v);

// Generate random value from shift-scaled cosine_power distribution with parameter v, mu, lambda
double generate_shift_scaled_cosine_power_value(double v, double mu, double lambda);




// Emperical functions




// Returns emperical central/not-central mometh with p-s order
double emperical_moment(std::vector<double>& data, int p, bool central);

// Returns emperical mathematical expectation
double emperical_mathematical_expectation(std::vector<double>& data);

// Returns emperical vatiance
double emperical_variance(std::vector<double>& data);

// Returns emperical asymmetry coefficient
double emperical_asymmetry(std::vector<double>& data);

// Returns emperical excess
double emperical_excess(std::vector<double>& data);

// Returns emperic density (For each segment exclude last segment left border included right border excluded. Last segment includes boths borders)
std::map<std::pair<double, double>, double> get_emperical_density(std::vector<double>& data);

// Generates random value from emperical distribution
double generate_from_emperical_distribution(std::map<std::pair<double, double>, double> func);