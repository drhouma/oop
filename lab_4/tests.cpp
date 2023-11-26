#include "standart_distribution.h"
#include "smesi.h"
#include "catch.h"

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
//-------------------------------------------------------------------------
//Density

TEST_CASE("[Cosine-power density: mu=0, lambda=1, v=0]", "[density]") {
    CosinePower unit{0, 0, 1};
    double expected = 0.500;

    auto func = unit.Density();
    CHECK(round(func(0) * 1000) / 1000 == expected);

}



TEST_CASE("density: [Cosine-power density: mu=0, lambda=1, v=1]", "[density]") {
    CosinePower unit{ 1, 0, 1 };
    double expected = 0.785;
    auto func = unit.Density();
    CHECK(round(func(0) * 1000) / 1000 == expected);

}
TEST_CASE("density: [Cosine-power density: mu=0, lambda=1, v=2]", "[density]") {
    double expected = 1.000;
    CosinePower unit{2, 0, 1};
    auto func = unit.Density();
    CHECK(round(func(0) * 1000) / 1000 == expected);

}
TEST_CASE("density: [Cosine-power density: mu=0, lambda=1, v=3]", "[density]") {
    double expected = 1.178;
    CosinePower unit{ 3, 0, 1 };
    auto func = unit.Density();
    CHECK(round(func(0) * 1000) / 1000 == expected);

}
//--------------------------------------------------------------------------------
// Math Excpectation
TEST_CASE("[Cosine-power ME: mu=0, lambda=1, v=0]", "[ME]") {
    double expected = 0;
    CosinePower unit{ 0, 0, 1 };
    auto res = unit.Expectation();
    CHECK(res == expected);

}

TEST_CASE("[Cosine-power ME: mu=1, lambda=1, v=0]", "[ME]") {
    double expected = 1;
    CosinePower unit{ 0, 1, 1 };
    auto res = unit.Expectation();
    CHECK(res == expected);

}

//--------------------------------------------------------------------------------
//Variance

TEST_CASE("[Cosine-power variance: mu=0, lambda=1, v=0]", "[variance]") {
    double expected = 0.333;
    CosinePower unit{ 0, 0, 1 };
    auto variance = unit.Variance();
    CHECK(round(variance * 1000) / 1000 == expected);

}

TEST_CASE("[Cosine-power variance: mu=0, lambda=1, v=1]", "[variance]") {
    double expected = 0.189;
    CosinePower unit{ 1, 0, 1 };
    auto variance = unit.Variance();
    CHECK(round(variance * 1000) / 1000 == expected);

}
TEST_CASE("[Cosine-power variance: mu=0, lambda=1, v=2]", "[variance]") {
    double expected = 0.131;
    CosinePower unit{ 2, 0, 1 };
    auto variance = unit.Variance();
    CHECK(round(variance * 1000) / 1000 == expected);

}

TEST_CASE("[Cosine-power variance: mu=0, lambda=1, v=3]", "[variance]") {
    double expected = 0.099;
    CosinePower unit{ 3, 0, 1 };
    auto variance = unit.Variance();
    CHECK(round(variance * 1000) / 1000 == expected);

}



//---------------------------------------------------------------------------------------------
//Excess

TEST_CASE("[Cosine-power excess: mu=0, lambda=1, v=0]", "[excess]") {
    double expected = -1.200;
    CosinePower unit{ 0, 0, 1 };
    auto variance = unit.Excess();
    CHECK(round(variance * 1000) / 1000 == expected);

}

TEST_CASE("[Cosine-power excess: mu=0, lambda=1, v=1]", "[excess]") {
    double expected = -0.806;
    CosinePower unit{ 1, 0, 1 };
    auto variance = unit.Excess();
    CHECK(round(variance * 1000) / 1000 == expected);

}


TEST_CASE("[Cosine-power excess: mu=0, lambda=1, v=2]", "[excess]") {
    double expected = -0.594;
    CosinePower unit{ 2, 0, 1 };
    auto variance = unit.Excess();
    CHECK(round(variance * 1000) / 1000 == expected);

}


TEST_CASE("[Cosine-power excess: mu=0, lambda=1, v=3]", "[excess]") {
    double expected = -0.466;
    CosinePower unit{3, 0, 1 };
    auto variance = unit.Excess();
    CHECK(round(variance * 1000) / 1000 == expected);

}

//----------------------------------------------------------
//Mixtures

//getters

TEST_CASE("[Get 1 component for mixtures]", "[getters]") {
    CosinePower d1{4.0, 3.0, 1.5};
    CosinePower d2{3.0, 5.0, 2.5};
    MixtureDistribution mix{d1, d2, 0.5};
    CHECK(mix.GetFirstFunction() == d1);
    CHECK(mix.GetSecondFunction() == d2);

}

//Density
TEST_CASE("[Mixtures density: mu1=mu2=1, l1=l2=2, p=0.3, v=6]", "[smesi]"){
    CosinePower d1{6.0, 1.0, 2.0};
    CosinePower d2{6.0, 1.0, 2.0};
    MixtureDistribution mix{d1, d2, 0.3};
    double expected = 0.1;
    auto func = mix.Density();
    CHECK(round(func(0)*1000)/1000 == expected);
}

//--------------------------------------------------------------------------------
//Math.Expectation
TEST_CASE("[Math.Expectation: mu1=3, l1=1.5, v1=4; mu2=5, l2=2.5, v2=3; p=0.5]", "[smesi]"){
    CosinePower d1{4.0, 3.0, 1.5};
    CosinePower d2{3.0, 5.0, 2.5};
    MixtureDistribution mix{d1, d2, 0.5};
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
    MixtureDistribution mix{d1, d2, 0.5};
    double expected = 0.497;
    double res = mix.Variance();
    CHECK(round(res*1000)/1000 == expected);
}

//Asymmetry

TEST_CASE("[Asymmetry: mu1=mu2=0, l1=1, l2=3, v1=v2=3, p=0.5]", "[smesi]"){
    CosinePower d1{4.0, 3.0, 1.5};
    CosinePower d2{3.0, 5.0, 2.5};
    MixtureDistribution mix{d1, d2, 0.5};
    double expected = 0.399;
    double res = mix.Asymmetry();
    CHECK(round(res*1000)/1000 == expected);
}
