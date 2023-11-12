#pragma once
#include "smesi.h"
#include "standart_distribution.h"
#include <map>
#include <algorithm>

class Empiric{
    int n; // объем выборки
    int k; // количество интервалов
    double* data; // массив данных
    std::map<std::pair<double, double>, double> fr;

    public:
    // для k = 1 - формула Стерджесса k >= 1
    // n > 2
    Empiric(int n0, CosinePower& prim, int k0=1); // три конструктора
    Empiric(int n0, MixtureDistribution& mixt, int k0=1); // с моделированием
    Empiric(int n0, Empiric& emp, int k0=1); // случайных величин
    // конструктор копирования и оператор присваивания
    // для глубокого копирования
    Empiric(const Empiric& emp);
    Empiric& operator=(const Empiric & emp);


    double* GetData() {return data;}
    std::map<std::pair<double, double>, double> GetFr() {return fr;};


    // Returns math expectation of cosine-power distribution
    double Expectation();
    // Returns variance of cosine-power distribution with
    double Variance();
    // Returns variance of cosine-power distribution
    double Asymmetry();
    // Returns excess of cosine-power distribution
    double Excess();

    function<double(double)> Density();

   

    double GenerateValue();

    //Ключ - интервал, значение - относительная частота
    std::map<std::pair<double, double>, double> GetEmpericalDensity();

    double emperical_moment(int p, bool central);

    ~Empiric(){
        delete [] data;
    }
};

