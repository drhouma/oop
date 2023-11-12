#pragma once
#include "smesi.h"
#include "standart_distribution.h"

class Empiric{
    int n; // объем выборки
    int k; // количество интервалов
    double* data; // массив данных
    double* fr;// массив относительных частот

    public:
    // для k = 1 - формула Стерджесса k >= 1
    // n > 2
    Empiric(int n0, CosinePower& prim, int k0); // три конструктора
    Empiric(int n0, MixtureDistribution& mixt, int k0); // с моделированием
    Empiric(int n0, Empiric& emp, int k0); // случайных величин
    // конструктор копирования и оператор присваивания
    // для глубокого копирования
    Empiric(const Empiric& emp);
    Empiric& operator=(const Empiric & emp);


    // Returns math expectation of cosine-power distribution
    double Expectation();
    // Returns variance of cosine-power distribution with
    double Variance();
    // Returns variance of cosine-power distribution
    double Asymmetry();
    // Returns excess of cosine-power distribution
    double Excess();
    

    double emperical_moment(int p, bool central);

    ~Empiric(){
        delete [] data;
        delete [] fr;
    }
};

