#include "../smesi.h"
#include "../catch.h"

using namespace std;

//---------------------------------------------------------------------------
//Density

TEST_CASE("[Mixtures density: mu1=mu2=1, l1=l2=2, p=0.3, v=0]", "[smesi]"){
    distribution d1{0.0, 1.0, 2.0};
    distribution d2{0.0, 1.0, 2.0};
    double expected = 0.250;
    auto func = mixture_density(d1, d2, 0.3);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Mixtures density: mu1=mu2=5, l1=l2=2, p=0.3, v=0]", "[smesi]"){
    distribution d1{0.0, 5.0, 2.0};
    distribution d2{0.0, 5.0, 2.0};
    double expected = 0.0;
    auto func = mixture_density(d1, d2, 0.3);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Mixtures density: mu1=mu2=1, l1=l2=2, p=0.3, v=1]", "[smesi]"){
    distribution d1{1.0, 1.0, 2.0};
    distribution d2{1.0, 1.0, 2.0};
    double expected = 0.278;
    auto func = mixture_density(d1, d2, 0.3);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Mixtures density: mu1=mu2=1, l1=l2=2, p=0.3, v=2]", "[smesi]"){
    distribution d1{2.0, 1.0, 2.0};
    distribution d2{2.0, 1.0, 2.0};
    double expected = 0.25;
    auto func = mixture_density(d1, d2, 0.3);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Mixtures density: mu1=mu2=1, l1=l2=2, p=0.3, v=3]", "[smesi]"){
    distribution d1{3.0, 1.0, 2.0};
    distribution d2{3.0, 1.0, 2.0};
    double expected = 0.208;
    auto func = mixture_density(d1, d2, 0.3);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Mixtures density: mu1=mu2=1, l1=l2=2, p=0.3, v=4]", "[smesi]"){
    distribution d1{4.0, 1.0, 2.0};
    distribution d2{4.0, 1.0, 2.0};
    double expected = 0.167;
    auto func = mixture_density(d1, d2, 0.3);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Mixtures density: mu1=mu2=1, l1=l2=2, p=0.3, v=5]", "[smesi]"){
    distribution d1{5.0, 1.0, 2.0};
    distribution d2{5.0, 1.0, 2.0};
    double expected = 0.130;
    auto func = mixture_density(d1, d2, 0.3);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Mixtures density: mu1=mu2=1, l1=l2=2, p=0.3, v=6]", "[smesi]"){
    distribution d1{6.0, 1.0, 2.0};
    distribution d2{6.0, 1.0, 2.0};
    double expected = 0.1;
    auto func = mixture_density(d1, d2, 0.3);
    CHECK(round(func(0)*1000)/1000 == expected);
}

//--------------------------------------------------------------------------------
//Math.Expectation
TEST_CASE("[Math.Expectation: mu1=3, l1=1.5, v1=4; mu2=5, l2=2.5, v2=3; p=0.5]", "[smesi]"){
    distribution d1{4.0, 3.0, 1.5};
    distribution d2{3.0, 5.0, 2.5};
    double expected = 4.0;
    double res = mixture_mathematical_expectation(d1, d2, 0.5);
    double exc = mixture_mathematical_excess(d1,d2,0.5);
    CHECK(round(res*1000)/1000 == expected);
    CHECK(exc != 0);
}

//--------------------------------------------------------------------------------
//Variance

TEST_CASE("[Variance: mu1=mu2=0, l1=1, l2=3, v1=v2=3, p=0.5]", "[smesi]"){
    distribution d1{3.0, 0.0, 1};
    distribution d2{3.0, 0.0, 3};
    double expected = 0.497;
    double res = mixture_variance(d1, d2, 0.5);
    CHECK(round(res*1000)/1000 == expected);
}