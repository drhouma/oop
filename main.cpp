#define CATCH_CONFIG_RUNNER
#include<iostream>
#include<fstream>
#include "catch.h"
#include "statistics.h"

void generate_theoretical_data(double v) {
	// Exports data for plotting theoretical density
	auto func = cosine_power_density(v);
	std::ofstream file("theoretical_values.txt");
	for (double val = -0.99; abs(val - 1) > 0.0001; val += 0.01)
		file << func(val) << std::endl;
	file.close();
}

void generate_emperic_data(int n, double v) {
	// Exports data for plotting emperic distribution
	std::ofstream file("emperical_values.txt");
	for (int i = 0; i <= 100000; ++i)
		file << generate_cosine_power_value(v) << std::endl;
	file.close();
}

void generate_emperical_density(int n, double v) {
	// Exports emperical density for plotting
	std::vector<double> data(n);
	for (int i = 0; i < n; ++i) {
		data[i] = generate_cosine_power_value(v);
	}
	auto func = get_emperic_density(data);
	std::ofstream file("emperical_density.txt");
	for (auto segment : func) {
		file << segment.first.first << '\t' << segment.first.second << '\t' << segment.second << std::endl;
	}
	file.close();
}

int main(int argc, char* argv[]) {
    srand(time(0)); // random state for random

    double v = 1.5; // parameter v for distribution density
    int n = 1000;

    generate_emperic_data(n, v);		// Generating emperic data
    generate_theoretical_data(v);	// Generating theoretical data
    generate_emperical_density(n, v);

    int result = Catch::Session().run(argc, argv);
    return result;
}