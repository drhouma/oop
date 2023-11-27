#include "standart_distribution.h"
#include "smesi.h"
#include "empiric.h"
#include "catch.h"
#include "interface.h"

using namespace std;

//-------------------------------------------------------------------------
// getters
TEST_CASE("[Get_form]", "[getters]") {
    CosinePower unit{ 2, 0, 1 };
    double expected = 2;
    auto res = unit.GetForm();
    CHECK(res == expected);
}

TEST_CASE("[Get_shift]", "[getters]") {
    CosinePower unit{ 0, 3, 1 };
    double expected = 3;
    auto res = unit.GetShift();
    CHECK(res == expected);

}
TEST_CASE("[Get_scale]", "[getters]") {
    CosinePower unit{ 0, 3, 10 };
    double expected = 10;
    auto res = unit.GetScale();
    CHECK(res == expected);

}
//Using interfaces
//-------------------------------------------------------------------------
//Density

TEST_CASE("density: [Cosine-power density: mu=0, lambda=1, v=3]", "[density]") {
    DInterface* interface;
    double expected = 1.178;
    interface = new CosinePower{ 3, 0, 1 };
    auto func = interface->Density();
    CHECK(round(func(0) * 1000) / 1000 == expected);

}
//--------------------------------------------------------------------------------
// Math Excpectation
TEST_CASE("[Cosine-power ME: mu=0, lambda=1, v=0]", "[ME]") {
    DInterface* interface = new CosinePower{ 0, 0, 1 };
    double expected = 0;
    auto res = interface->Expectation();
    CHECK(res == expected);

}

TEST_CASE("[Cosine-power ME: mu=1, lambda=1, v=0]", "[ME]") {
    double expected = 1;
    DInterface* interface = new CosinePower{ 0, 1, 1 };
    auto res = interface->Expectation();
    CHECK(res == expected);

}

//--------------------------------------------------------------------------------
//Variance


TEST_CASE("[Cosine-power variance: mu=0, lambda=1, v=3]", "[variance]") {
    double expected = 0.099;
    DInterface* interface = new CosinePower{ 3, 0, 1 };
    auto variance = interface->Variance();
    CHECK(round(variance * 1000) / 1000 == expected);

}



//---------------------------------------------------------------------------------------------
//Excess


TEST_CASE("[Cosine-power excess: mu=0, lambda=1, v=3]", "[excess]") {
    double expected = -0.466;
    DInterface* interface = new CosinePower{3, 0, 1 };
    auto variance = interface->Excess();
    CHECK(round(variance * 1000) / 1000 == expected);

}

//----------------------------------------------------------
//Mixtures

//getters

TEST_CASE("[Get 1 component for mixtures]", "[getters]") {
    CosinePower d1{4.0, 3.0, 1.5};
    CosinePower d2{3.0, 5.0, 2.5};
    MixtureDistribution<CosinePower, CosinePower> mix{d1, d2, 0.5};
    CHECK(mix.GetFirstFunction() == d1);
    CHECK(mix.GetSecondFunction() == d2);

}

//Density
TEST_CASE("[Mixtures density: mu1=mu2=1, l1=l2=2, p=0.3, v=6]", "[smesi]"){
    CosinePower d1{6.0, 1.0, 2.0};
    CosinePower d2{6.0, 1.0, 2.0};
    MixtureDistribution<CosinePower, CosinePower> mix{d1, d2, 0.3};
    double expected = 0.1;
    auto func = mix.Density();
    CHECK(round(func(0)*1000)/1000 == expected);
}

//--------------------------------------------------------------------------------
//Math.Expectation
TEST_CASE("[Math.Expectation: mu1=3, l1=1.5, v1=4; mu2=5, l2=2.5, v2=3; p=0.5]", "[smesi]"){
    CosinePower d1{4.0, 3.0, 1.5};
    CosinePower d2{3.0, 5.0, 2.5};
    MixtureDistribution<CosinePower, CosinePower> mix{d1, d2, 0.5};
    double expected = 4.0;
    double res = mix.Expectation();
    double exc = mix.Excess();
    CHECK(round(res*1000)/1000 == expected);
    CHECK(exc != 0);
}

//--------------------------------------------------------------------------------
//Variance

TEST_CASE("[Variance: mu1=mu2=0, l1=1, l2=3, v1=v2=3, p=0.5]", "[smesi]"){
    CosinePower d1{3.0, 0.0, 1};
    CosinePower d2{3.0, 0.0, 3};
    MixtureDistribution<CosinePower, CosinePower> mix{d1, d2, 0.5};
    double expected = 0.497;
    double res = mix.Variance();
    CHECK(round(res*1000)/1000 == expected);
}

//Asymmetry

TEST_CASE("[Asymmetry: mu1=3, mu2=5, l1=1.5, l2=2.5, v1=4, v2=3, p=0.5]", "[smesi]"){
    CosinePower d1{4.0, 3.0, 1.5};
    CosinePower d2{3.0, 5.0, 2.5};
    MixtureDistribution<CosinePower, CosinePower> mix{d1, d2, 0.5};
    double expected = 0.399;
    double res = mix.Asymmetry();
    CHECK(round(res*1000)/1000 == expected);
}

//Empiric
TEST_CASE("[Empiric Expectation for cosine power; n = 100, v =4, mu=3, l=1.5]", "[smesi]"){
    CosinePower cs{4.0, 3.0, 1.5};

    Empiric emp(100, cs);
    double expected = 2.957;
    auto res = emp.Expectation();
    CHECK(round(res*1000)/1000 == expected);
}

TEST_CASE("[Empiric Asymmetry for mixture distribution (using CosinePower); n = 100,  mu1=3, mu2=5, l1=1.5, l2=2.5, v1=4, v2=3, p=0.5]", "[smesi]"){
    CosinePower d1{4.0, 3.0, 1.5};
    CosinePower d2{3.0, 5.0, 2.5};
    MixtureDistribution<CosinePower, CosinePower> mix{d1, d2, 0.5};
    Empiric emp(100, mix);
    double expected = 0.317;
    auto res = emp.Asymmetry();
    CHECK(round(res*1000)/1000 == expected);
}

TEST_CASE("[Empiric Excess for empiric distribution for Mixtures (using CosinePower); n = 100,  mu1=3, mu2=5, l1=1.5, l2=2.5, v1=4, v2=3, p=0.5]", "[smesi]"){
    CosinePower d1{4.0, 3.0, 1.5};
    CosinePower d2{3.0, 5.0, 2.5};
    MixtureDistribution<CosinePower, CosinePower> mix{d1, d2, 0.5};
    Empiric emp(100, mix);
    Empiric emp2(100, emp);
    double expected = -0.7;
    auto res = emp.Excess();
    CHECK(round(res*1000)/1000 == expected);
}
