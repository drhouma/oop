#define CATCH_CONFIG_RUNNER
#include<iostream>
#include<fstream>

#include "statistics.h"
#include "smesi.h"
#include "catch.h"

void export_theoretical_data(double v) {
    // Exports data for plotting theoretical density
    auto func = cosine_power_density(v);
    std::ofstream file("data/1_1.txt");
    for (double val = -0.99; abs(val - 1) > 0.0001; val += 0.01)
        file << func(val) << std::endl;
    file.close();
}

void export_emperical_data(int n, double v) {
    // Exports data for plotting emperic distribution
    std::ofstream file("data/2_1.txt");
    for (int i = 0; i <= 100000; ++i)
        file << generate_cosine_power_value(v) << std::endl;
    file.close();
}

void export_emperical_density(int n, double v) {
    // Exports emperical density for plotting
    std::vector<double> data(n);
    for (int i = 0; i < n; ++i) {
        data[i] = generate_cosine_power_value(v);
    }
    auto func = get_emperical_density(data);
    std::ofstream file("data/3_1.txt");
    for (auto& segment : func) {
        file << segment.first.first << '\t' << segment.first.second << '\t' << segment.second << std::endl;
    }
    file.close();
}

void export_mixture_emperical_density(int n, distribution f1, distribution f2, double p) {

    std::ofstream file1("data/4_1.txt");
    std::vector<double> data(n);
    for (int i = 1; i <= n; ++i) {
        double val = mixture_generate_value(f1, f2, p);
        data[i - 1] = val;
        file1 << val << std::endl;
    }
    file1.close();

    auto func = get_emperical_density(data);
    std::ofstream file2("data/4_2.txt");
    for (auto& segment : func) {
        file2 << segment.first.first << '\t' << segment.first.second << '\t' << segment.second << std::endl;
    }
    file2.close();
}

void export_from_emperical_distribution_data(int n, distribution f) {
    std::ofstream file("data/5_1.txt");
    std::vector<double> data(n);
    for (int i = 0; i < n; ++i) {
        data[i] = generate_shift_scaled_cosine_power_value(f.v, f.mu, f.lambda);
        file << data[i] << std::endl;
    }
    file.close();

    auto emperical_func = get_emperical_density(data);
    std::ofstream file1("data/5_2.txt");
    for (int i = 0; i < n; ++i) {
        file1 << generate_from_emperical_distribution(emperical_func) << std::endl;
    }
    file1.close();

    std::ofstream file2("data/5_3.txt");
    for (auto& segment : emperical_func) {
        file2 << segment.first.first << '\t' << segment.first.second << '\t' << segment.second << std::endl;
    }
    file2.close();

    auto theor = cosine_power_density(f.v);
    auto func = shift_scaled_density(theor, f.mu, f.lambda);
    std::ofstream file3("data/5_4.txt");
    double step = (2 * f.lambda + f.mu) / n;
    for (double val = f.mu-f.lambda; abs(val - (f.mu+f.lambda)) > 0.0001; val += step)
        file3 << val << ' ' << func(val) << std::endl;
    file3.close();
}

int main(int argc, char* argv[]) {
    srand(time(0)); // random state for random

    double v = 0.5; // parameter v for distribution density
    int n = 10000;
    distribution f1 = {0.7, -1, 2.0}, f2 = {0.3, 1, 0.5};
    double p = 0.7;

    export_emperical_data(n, v);		// Generating emperical data
    export_theoretical_data(v);	// Generating theoretical data
    export_emperical_density(n, v); // Generating emperical density
    export_mixture_emperical_density(n, f1, f2, p); // Generating mixture data
    export_from_emperical_distribution_data(n, { 1.5, 2, 0.5 }); // Generating data from emperical distribution

    int result = Catch::Session().run(argc, argv);
    return result;
}