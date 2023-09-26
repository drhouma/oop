#define CATCH_CONFIG_RUNNER
#include<iostream>
#include<fstream>

#include "statistics.h"
#include "smesi.h"
#include "catch.h"

void export_main_distribution_graph_data(int n, double v) {
    // Записываем теоретическую плотность
    auto func = cosine_power_density(v);
    std::ofstream file1("data/1_theor_dens.txt");
    file1 << v << std::endl;
    for (double val = -0.99; abs(val - 1) > 0.001 && val < 1; val += 0.01)
        file1 << val << ' ' <<  func(val) << std::endl;
    file1.close();
    // Генерируем случайные x
    std::vector<double> data(n);

    for (int i = 0; i < n; ++i)
        data[i] =  generate_cosine_power_value(v);
    auto emperic_dens = get_emperical_density(data);
    // Записываем граничные значения эмпирической плотности
    std::ofstream file2("data/1_emperical_dens.txt");
    for (auto& segment : emperic_dens) {
        file2 << segment.first.first << ' ' << segment.second << std::endl;
        file2 << segment.first.second << ' ' << segment.second << std::endl;
    }
    file2.close();
    std::ofstream file3("data/1_generated_dots.txt");
    // Записываем значения сгенерированных x и соответствующие значения теор
    // плотности и эмпирической плотности в данных x
    for (auto x : data) {
        file3 << x << ' ' << func(x) << ' ' << get_value_of_emperical_density(emperic_dens, x) << std::endl;
    }
    file3.close();
}

void export_shift_scaled_distribution_draph_data(int n, distribution params) {
    // Записываем теоретическую плотность
    auto main_func = cosine_power_density(params.v);
    auto func = shift_scaled_density(main_func, params.mu, params.lambda);
    std::ofstream file1("data/2_theor_dens.txt");
    file1 << params.v << ' ' << params.mu << ' ' << params.lambda << std::endl;
    std::vector<double> elements = { params.mu - params.lambda, params.mu+params.lambda };
    double left = *(std::min_element(elements.begin(), elements.end()));
    double right = *(std::max_element(elements.begin(), elements.end()));
    for (double val = left; abs(val - right) > 0.001 && val < right; val += (right - left) / double(1000)) {
        file1 << val << ' ' << func(val) << std::endl;
    }
    file1.close();
    // Генерируем случайные x
    std::vector<double> data(n);

    for (int i = 0; i < n; ++i)
        data[i] = generate_shift_scaled_cosine_power_value(params.v, params.mu, params.lambda);
    auto emperic_dens = get_emperical_density(data);
    // Записываем граничные значения эмпирической плотности
    std::ofstream file2("data/2_emperical_dens.txt");
    for (auto& segment : emperic_dens) {
        file2 << segment.first.first << ' ' << segment.second << std::endl;
        file2 << segment.first.second << ' ' << segment.second << std::endl;
    }
    file2.close();
    std::ofstream file3("data/2_generated_dots.txt");
    // Записываем значения сгенерированных x и соответствующие значения теор
    // плотности и эмпирической плотности в данных x
    for (auto x : data) {
        file3 << x << ' ' << func(x) << ' ' << get_value_of_emperical_density(emperic_dens, x) << std::endl;
    }
    file3.close();
}

void export_mixture_distribution_graph_data(int n, distribution f1, distribution f2, double p) {
    // Записываем теоретическую плотность
    auto func = mixture_density(f1, f2, p);
    std::ofstream file1("data/3_theor_dens.txt");
    file1 << f1.v << ' ' << f1.mu << ' ' << f1.lambda << ' ' << f2.v << ' ' << f2.mu << ' ' << f2.lambda << ' ' << p << std::endl;
    std::vector<double> elements = { f1.mu - f1.lambda, f1.mu + f1.lambda, f2.mu - f2.lambda, f2.mu + f2.lambda };
    double left = *(std::min_element(elements.begin(), elements.end()));
    double right = *(std::max_element(elements.begin(), elements.end()));
    for (double val = left+0.001; abs(val - right) > 0.001 && val < right; val += (right-left) / double(1000)) {
        file1 << val << ' ' << func(val) << std::endl;
    }
    file1.close();
    // Генерируем случайные x
    std::vector<double> data(n);

    for (int i = 0; i < n; ++i)
        data[i] = mixture_generate_value(f1, f2, p);
    auto emperic_dens = get_emperical_density(data);
    // Записываем граничные значения эмпирической плотности
    std::ofstream file2("data/3_emperical_dens.txt");
    for (auto& segment : emperic_dens) {
        file2 << segment.first.first << ' ' << segment.second << std::endl;
        file2 << segment.first.second << ' ' << segment.second << std::endl;
    }
    file2.close();
    std::ofstream file3("data/3_generated_dots.txt");
    // Записываем значения сгенерированных x и соответствующие значения теор
    // плотности и эмпирической плотности в данных x
    for (auto x : data) {
        file3 << x << ' ' << func(x) << ' ' << get_value_of_emperical_density(emperic_dens, x) << std::endl;
    }
    file3.close();
}

void export_emperical_from_emperical_graph_data(int n, distribution f1) {
    // Записываем теоретическую плотность
    auto main_func = cosine_power_density(f1.v);
    auto func = shift_scaled_density(main_func, f1.mu, f1.lambda);
    std::ofstream file1("data/4_theor_dens.txt");
    file1 << f1.v << ' ' << f1.mu << ' ' << f1.lambda << std::endl;
    std::vector<double> elements = { f1.mu - f1.lambda, f1.mu + f1.lambda };
    double left = *(std::min_element(elements.begin(), elements.end()));
    double right = *(std::max_element(elements.begin(), elements.end()));
    for (double val = left + 0.001; abs(val - right) > 0.001 && val < right; val += (right - left) / double(1000)) {
        file1 << val << ' ' << func(val) << std::endl;
    }
    file1.close();
    // Генерируем случайные x
    std::vector<double> data(n);

    for (int i = 0; i < n; ++i)
        data[i] = generate_shift_scaled_cosine_power_value(f1.v, f1.mu, f1.lambda);
    auto emperic_dens = get_emperical_density(data);
    // Записываем граничные значения эмпирической плотности первой
    std::ofstream file2("data/4_emperical_dens.txt");
    for (auto& segment : emperic_dens) {
        file2 << segment.first.first << ' ' << segment.second << std::endl;
        file2 << segment.first.second << ' ' << segment.second << std::endl;
    }
    file2.close();
    for (int i = 0; i < n; ++i) {
        data[i] = generate_from_emperical_distribution(emperic_dens);
    }
    auto emperic_dens2 = get_emperical_density(data);
    // Записываем вторую эмперическую плотность
    std::ofstream file3("data/4_emperical_dens2.txt");
    for (auto& segment : emperic_dens2) {
        file3 << segment.first.first << ' ' << segment.second << std::endl;
        file3 << segment.first.second << ' ' << segment.second << std::endl;
    }
    file3.close();
    std::ofstream file4("data/4_generated_dots.txt");
    // Записываем значения сгенерированных x и соответствующие значения теор
    // плотности и первой и второй эмпирических плотностей в данных x
    for (auto x : data) {
        file4 << x << ' ' << func(x) << ' ' << get_value_of_emperical_density(emperic_dens, x)
              << ' ' << get_value_of_emperical_density(emperic_dens2, x) << std::endl;
    }
    file4.close();
}

int main(int argc, char* argv[]) {
    srand(time(0)); // random state for random

    export_main_distribution_graph_data(250, 3);
    export_shift_scaled_distribution_draph_data(250, { 2, 5, 0.3 });
    export_mixture_distribution_graph_data(1000, { 2, -3, 1.5 }, { 0.5, 2, 0.5 }, 0.3);
    export_emperical_from_emperical_graph_data(1000, { 4, 2, 2.2 });

    int result = Catch::Session().run(argc, argv);
    return result;
    return 0;
}