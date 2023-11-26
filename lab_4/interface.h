#pragma once
#include <fstream>
#include <iostream>
#include<cmath>
#include<functional>
#include<vector>
#include<stdexcept>
#include <string>
#include <map>
#include <algorithm>


class DInterface
{
private:
    /* data */
public:
    DInterface(/* args */) {}
    ~DInterface() {}

    virtual std::function<double(double)> Density() = 0;
    // Returns math expectation of cosine-power distribution
    virtual double Expectation() = 0;

    // Returns variance of cosine-power distribution with
    virtual double Variance() = 0;

    // Returns variance of cosine-power distribution
    virtual double Asymmetry() = 0;

    // Returns excess of cosine-power distribution
    virtual double Excess() = 0;

    virtual double GenerateValue() = 0;

};

class PInterface {
    public:
    PInterface(/* args */) {}
    ~PInterface() {}

    virtual void Save(std::ofstream& file) = 0;
    virtual void Load(std::ifstream& file, std::vector<double> &options) = 0;

};

