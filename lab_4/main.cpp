#define CATCH_CONFIG_RUNNER
#include "standart_distribution.h"
#include "empiric.h"
#include "smesi.h"
#include <iostream>
#include "catch.h"


// double get_value_of_emperical_density(std::map<std::pair<double, double>, double>& func, double x) {
//     for (auto& segment : func) {
//         if ((segment.first.first <= x) && (x <= segment.first.second))
//             return segment.second;
//     }
//     return 0.0;
// }


// //Graphics
// void export_main_distribution_graph_data(int n, double v) {
//     // Записываем теоретическую плотность
//     CosinePower obj{v};
//     auto func = obj.Density();
//     std::ofstream file1("data/1_theor_dens.txt");
//     file1 << v << std::endl;
//     for (double val = -0.99; abs(val - 1) > 0.001 && val < 1; val += 0.01)
//         file1 << val << ' ' <<  func(val) << std::endl;
//     file1.close();
//     // Генерируем случайные x
//     Empiric emp(n, obj);
//     auto emperic_dens = emp.GetEmpericalDensity();
//     // Записываем граничные значения эмпирической плотности
//     std::ofstream file2("data/1_emperical_dens.txt");
//     for (auto& segment : emperic_dens) {
//         file2 << segment.first.first << ' ' << segment.second << std::endl;
//         file2 << segment.first.second << ' ' << segment.second << std::endl;
//     }
//     file2.close();
//     std::ofstream file3("data/1_generated_dots.txt");
//     // Записываем значения сгенерированных x и соответствующие значения теор
//     // плотности и эмпирической плотности в данных x
//     auto empFunc = emp.Density();
//     auto data = emp.GetData();
//     for (int i = 0; i < n; i++) {
//         file3 << data[i] << ' ' << data[i] << ' ' << empFunc(data[i]) << std::endl;
//     }
//     file3.close();
// }


// void export_emperical_from_emperical_graph_data(int n, CosinePower obj) {
//     // Записываем теоретическую плотность
//     double v = obj.GetForm(), mu = obj.GetShift(), lambda = obj.GetScale();
//     auto func = obj.Density();
//     //auto func = shift_scaled_density(main_func, mu, lambda);
//     std::ofstream file1("data/4_theor_dens.txt");
//     file1 << v << ' ' << mu << ' ' << lambda << std::endl;

//     std::vector<double> elements = { mu - lambda, mu + lambda };
//     double left = *(std::min_element(elements.begin(), elements.end()));
//     double right = *(std::max_element(elements.begin(), elements.end()));
//     for (double val = left + 0.001; abs(val - right) > 0.001 && val < right; val += (right - left) / double(1000)) {
//         file1 << val << ' ' << func(val) << std::endl;
//     }
//     file1.close();
//     // Генерируем случайные x
//     Empiric emp(n, obj, sturges_rule(n));
//     auto emperic_dens = emp.GetEmpericalDensity();
//     // Записываем граничные значения эмпирической плотности первой
//     std::ofstream file2("data/4_emperical_dens.txt");
//     for (auto& segment : emperic_dens) {
//         file2 << segment.first.first << ' ' << segment.second << std::endl;
//         file2 << segment.first.second << ' ' << segment.second << std::endl;
//     }
//     file2.close();
// }

// void export_mixture_distribution_graph_data(int n, CosinePower f1, CosinePower f2, double p) {
//     // Записываем теоретическую плотность
//     MixtureDistribution md(f1, f2, p);
//     auto func = md.Density();
//     std::ofstream file1("data/3_theor_dens.txt");
//     file1 << f1.GetForm() << ' ' << f1.GetShift() << ' ' << f1.GetScale() << ' ' << f2.GetForm() << ' ' << f2.GetShift() << ' ' << f2.GetScale() << ' ' << p << std::endl;
//     std::vector<double> elements = { f1.GetShift() - f1.GetScale(), f1.GetShift() + f1.GetScale(), f2.GetShift() - f2.GetScale(), f2.GetShift() + f2.GetScale() };
//     double left = *(std::min_element(elements.begin(), elements.end()));
//     double right = *(std::max_element(elements.begin(), elements.end()));
//     for (double val = left+0.001; abs(val - right) > 0.001 && val < right; val += (right-left) / double(1000)) {
//         file1 << val << ' ' << func(val) << std::endl;
//     }
//     file1.close();
//     // Генерируем случайные x

//     Empiric tmp(n, md, sturges_rule(n));
//     auto emperic_dens = tmp.GetEmpericalDensity();
//     // Записываем граничные значения эмпирической плотности
//     std::ofstream file2("data/3_emperical_dens.txt");
//     for (auto& segment : emperic_dens) {
//         file2 << segment.first.first << ' ' << segment.second << std::endl;
//         file2 << segment.first.second << ' ' << segment.second << std::endl;
//     }
//     file2.close();
//     std::ofstream file3("data/3_generated_dots.txt");
//     // Записываем значения сгенерированных x и соответствующие значения теор
//     // плотности и эмпирической плотности в данных x
//     auto data = tmp.GetData();
//     auto empFunc = tmp.Density();
//     for (int i = 0; i < n; i++) {
//         file3 << data[i] << ' ' << func(data[i]) << ' ' << empFunc(data[i]) << std::endl;
//     }
//     file3.close();
// }

// void export_shift_scaled_distribution_draph_data(int n, CosinePower d1) {
//     // Записываем теоретическую плотность
//     auto func = d1.Density();
//     std::ofstream file1("data/2_theor_dens.txt");
//     file1 << d1.GetForm() << ' ' << d1.GetShift() << ' ' << d1.GetScale() << std::endl;
//     std::vector<double> elements = { d1.GetShift() - d1.GetScale(), d1.GetShift()+d1.GetScale() };
//     double left = *(std::min_element(elements.begin(), elements.end()));
//     double right = *(std::max_element(elements.begin(), elements.end()));
//     for (double val = left; abs(val - right) > 0.001 && val < right; val += (right - left) / double(1000)) {
//         file1 << val << ' ' << func(val) << std::endl;
//     }
//     file1.close();
//     // Генерируем случайные x
//     Empiric emp(n, d1, sturges_rule(n));
//     auto emperic_dens = emp.GetEmpericalDensity();
//     // Записываем граничные значения эмпирической плотности
//     std::ofstream file2("data/2_emperical_dens.txt");
//     for (auto& segment : emperic_dens) {
//         file2 << segment.first.first << ' ' << segment.second << std::endl;
//         file2 << segment.first.second << ' ' << segment.second << std::endl;
//     }
//     file2.close();
//     std::ofstream file3("data/2_generated_dots.txt");
//     // Записываем значения сгенерированных x и соответствующие значения теор
//     // плотности и эмпирической плотности в данных x
//     auto data = emp.GetData();
//     auto empFunc = emp.Density();
//     for (int i = 0; i < n; i++) {
//         file3 << data[i] << ' ' << func(data[i]) << ' ' << empFunc(data[i]) << std::endl;
//     }
//     file3.close();
// }



int main(int argc, char* argv[]) {
    // try {
    //     CosinePower obj{1,1,2};
    //     ofstream file_out("options.txt");
    //     obj.Save(file_out);

    //     //Stream constructor
    //     ifstream  file_in("options.txt");
    //     CosinePower fromfile(file_in);
    //     cout << fromfile.GetForm() << " " << fromfile.GetShift() << " " << fromfile.GetScale() << endl;
    //     file_in.close();

    //     //Load
    //     vector<double> options(5);
    //     CosinePower fromload{};
    //     file_in.open("options.txt");
    //     fromload.Load(file_in, options);
    //     cout << fromload.GetForm() << " " << fromload.GetShift() << " " << fromload.GetScale() << endl;
    //     for(auto elem: options){
    //         cout << elem << endl;
    //     }

    //     //Generate
    //     cout << "Generate data:\n";
    //     cout << fromload.GenerateValue();

    //     file_out.close();
    //     file_in.close();

    //     CosinePower unit{3};
    //     export_emperical_from_emperical_graph_data(100000, unit);
    //     export_main_distribution_graph_data(100000, 3);
    //     export_mixture_distribution_graph_data(1000, CosinePower(2, -3, 1.5), CosinePower(0.5, 2, 0.5), 0.3);
    //     export_shift_scaled_distribution_draph_data(250, CosinePower(2, 5, 0.3));

    // }

    // catch (const invalid_argument& e)
    // {
    //     cerr << e.what() << endl;
    //     return -1;
    // }

    // int result = Catch::Session().run(argc, argv);
    // return result;


    // return 0;

    DInterface *interface;
    CosinePower p(1, 1, 2), p2(2,1,2);
    interface = new CosinePower(1, 1, 2);
    // or 
    // interface = &p;
    auto func = interface->Density();
    std::cout << "Virtual: "<< func(0.5) << " Obj: " << p.Density()(0.5) << std::endl;

    MixtureDistribution<CosinePower, CosinePower> md(p2, p, 0.5);
    md.Density();
    md.Asymmetry();

    delete interface;
}
