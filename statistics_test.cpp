#include "statistics.h"
#include "catch.h"

using namespace std;

//-------------------------------------------------------------------------

TEST_CASE("[Cosine-power density: mu=0, lambda=1, v=0]", "[statistics]") {
    double expected = 0.500;
    double v = 0;
    auto func = cosine_power_density(v);
    CHECK(round(func(0)*1000)/1000 == expected);

}

TEST_CASE("[Cosine-power density: mu=0, lambda=1, v=1]", "[statistics]") {
    double expected = 0.785;
    double v = 1;
    auto func = cosine_power_density(v);
    CHECK(round(func(0)*1000)/1000 == expected);

}
TEST_CASE("[Cosine-power density: mu=0, lambda=1, v=2]", "[statistics]") {
    double expected = 1.000;
    double v = 2;
    auto func = cosine_power_density(v);
    CHECK(round(func(0)*1000)/1000 == expected);

}
TEST_CASE("[Cosine-power density: mu=0, lambda=1, v=3]", "[statistics]") {
    double expected = 1.178;
    double v = 3;
    auto func = cosine_power_density(v);
    CHECK(round(func(0)*1000)/1000 == expected);

}
TEST_CASE("[Cosine-power density: mu=0, lambda=1, v=4]", "[statistics]") {
    double expected = 1.333;
    double v = 4;
    auto func = cosine_power_density(v);
    CHECK(round(func(0)*1000)/1000 == expected);

}
TEST_CASE("[Cosine-power density: mu=0, lambda=1, v=5]", "[statistics]") {
    double expected = 1.473;
    double v = 5;
    auto func = cosine_power_density(v);
    CHECK(round(func(0)*1000)/1000 == expected);

}
TEST_CASE("[Cosine-power density: mu=0, lambda=1, v=6]", "[statistics]") {
    double expected = 1.600;
    double v = 6;
    auto func = cosine_power_density(v);
    CHECK(round(func(0)*1000)/1000 == expected);

}

//--------------------------------------------------------------------------------

TEST_CASE("[Cosine-power variance: mu=0, lambda=1, v=0]", "[statistics]") {
    double expected = 0.333;
    double v = 0;
    auto variance = cosine_power_mathematical_variance(v, 0, 1);
    CHECK(round(variance*1000)/1000 == expected);

}
