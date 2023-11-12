#define CATCH_CONFIG_RUNNER
#include "standart_distribution.h"
#include "empiric.h"
#include "smesi.h"
#include <iostream>
#include "catch.h"


double get_value_of_emperical_density(std::map<std::pair<double, double>, double>& func, double x) {
    for (auto& segment : func) {
        if ((segment.first.first <= x) && (x <= segment.first.second))
            return segment.second;
    }
    return 0.0;
}


//Graphics
void export_main_distribution_graph_data(int n, double v) {
    // Записываем теоретическую плотность
    CosinePower obj{v};
    auto func = obj.Density();
    std::ofstream file1("data/1_theor_dens.txt");
    file1 << v << std::endl;
    for (double val = -0.99; abs(val - 1) > 0.001 && val < 1; val += 0.01)
        file1 << val << ' ' <<  func(val) << std::endl;
    file1.close();
    // Генерируем случайные x
    Empiric emp(n, obj);
    auto emperic_dens = emp.GetEmpericalDensity();
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
    for (int i = 0; i < n; i++) {
        file3 << emp.GetData()[i] << ' ' << func(emp.GetData()[i]) << ' ' << get_value_of_emperical_density(emperic_dens, emp.GetData()[i]) << std::endl;
    }
    file3.close();
}


void export_emperical_from_emperical_graph_data(int n, CosinePower obj) {
    // Записываем теоретическую плотность
    double v = obj.GetForm(), mu = obj.GetShift(), lambda = obj.GetScale();
    auto func = obj.Density();
    //auto func = shift_scaled_density(main_func, mu, lambda);
    std::ofstream file1("data/4_theor_dens.txt");
    file1 << v << ' ' << mu << ' ' << lambda << std::endl;

    std::vector<double> elements = { mu - lambda, mu + lambda };
    double left = *(std::min_element(elements.begin(), elements.end()));
    double right = *(std::max_element(elements.begin(), elements.end()));
    for (double val = left + 0.001; abs(val - right) > 0.001 && val < right; val += (right - left) / double(1000)) {
        file1 << val << ' ' << func(val) << std::endl;
    }
    file1.close();
    // Генерируем случайные x
    Empiric emp(n, obj);
    auto emperic_dens = emp.GetEmpericalDensity();
    // Записываем граничные значения эмпирической плотности первой
    std::ofstream file2("data/4_emperical_dens.txt");
    for (auto& segment : emperic_dens) {
        file2 << segment.first.first << ' ' << segment.second << std::endl;
        file2 << segment.first.second << ' ' << segment.second << std::endl;
    }
    file2.close();
}


int main(int argc, char* argv[]) {
    try {
        //Example
        //Save
        CosinePower obj{1,1,2};
        ofstream file_out("options.txt");
        obj.Save(file_out);

        //Stream constructor
        ifstream  file_in("options.txt");
        CosinePower fromfile(file_in);
        cout << fromfile.GetForm() << " " << fromfile.GetShift() << " " << fromfile.GetScale() << endl;
        file_in.close();

        //Load
        vector<double> options(5);
        CosinePower fromload{};
        file_in.open("options.txt");
        fromload.Load(file_in, options);
        cout << fromload.GetForm() << " " << fromload.GetShift() << " " << fromload.GetScale() << endl;
        for(auto elem: options){
            cout << elem << endl;
        }

        //Generate
        cout << "Generate data:\n";
        cout << fromload.Generate_cosine_power_value();

        file_out.close();
        file_in.close();

        CosinePower unit{3};
        export_emperical_from_emperical_graph_data(100000, unit);
        export_main_distribution_graph_data(100000, 3);

    }

    catch (const invalid_argument& e)
    {
        cerr << e.what() << endl;
        return -1;
    }

    int result = Catch::Session().run(argc, argv);
    return result;


    return 0;
}
