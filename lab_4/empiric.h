 #pragma once
 #include "smesi.h"
 #include "standart_distribution.h"
 #include "interface.h"
 #include <map>
 #include <algorithm>


 class Empiric: public DInterface, public PInterface{
     int n; // объем выборки
     int k; // количество интервалов
     double* data; // массив данных
     std::map<std::pair<double, double>, double> fr;

     public:
     Empiric(int n0, DInterface& d, int k0=1);
     // для глубокого копирования
     Empiric(const Empiric& emp);
     Empiric& operator=(const Empiric & emp);


     double* GetData() {return data;}
     std::map<std::pair<double, double>, double> GetFr() {return fr;};


     // Returns math expectation of cosine-power distribution
     double Expectation() override;
     // Returns variance of cosine-power distribution with
     double Variance() override;
     // Returns variance of cosine-power distribution
     double Asymmetry() override;
     // Returns excess of cosine-power distribution
     double Excess() override;

     function<double(double)> Density() override;



     double GenerateValue() override;

     void Save(ofstream& file) override;
     void Load(ifstream& file, vector<double> &options) override;

     //Ключ - интервал, значение - относительная частота
     std::map<std::pair<double, double>, double> GetEmpericalDensity();

     double emperical_moment(int p, bool central);

     ~Empiric(){
         delete [] data;
     }
 };

