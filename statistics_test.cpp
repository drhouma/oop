#include "statistics.h"
#include "catch.h"

using namespace std;

//-------------------------------------------------------------------------
//Density

TEST_CASE("[Cosine-power density: mu=0, lambda=1, v=0]", "[statistics]") {
    double expected = 0.500;
    double v = 0.0;
    auto func = cosine_power_density(v);
    CHECK(round(func(0)*1000)/1000 == expected);

}

TEST_CASE("[Cosine-power density: mu=0, lambda=1, v=1]", "[statistics]") {
    double expected = 0.785;
    double v = 1.0;
    auto func = cosine_power_density(v);
    CHECK(round(func(0)*1000)/1000 == expected);

}
TEST_CASE("[Cosine-power density: mu=0, lambda=1, v=2]", "[statistics]") {
    double expected = 1.000;
    double v = 2.0;
    auto func = cosine_power_density(v);
    CHECK(round(func(0)*1000)/1000 == expected);

}
TEST_CASE("[Cosine-power density: mu=0, lambda=1, v=3]", "[statistics]") {
    double expected = 1.178;
    double v = 3.0;
    auto func = cosine_power_density(v);
    CHECK(round(func(0)*1000)/1000 == expected);

}
TEST_CASE("[Cosine-power density: mu=0, lambda=1, v=4]", "[statistics]") {
    double expected = 1.333;
    double v = 4.0;
    auto func = cosine_power_density(v);
    CHECK(round(func(0)*1000)/1000 == expected);

}
TEST_CASE("[Cosine-power density: mu=0, lambda=1, v=5]", "[statistics]") {
    double expected = 1.473;
    double v = 5.0;
    auto func = cosine_power_density(v);
    CHECK(round(func(0)*1000)/1000 == expected);

}
TEST_CASE("[Cosine-power density: mu=0, lambda=1, v=6]", "[statistics]") {
    double expected = 1.600;
    double v = 6.0;
    auto func = cosine_power_density(v);
    CHECK(round(func(0)*1000)/1000 == expected);

}

//--------------------------------------------------------------------------------
//Variance

TEST_CASE("[Cosine-power variance: mu=0, lambda=1, v=0]", "[statistics]") {
    double expected = 0.333;
    double v = 0.0;
    auto variance = cosine_power_mathematical_variance(v, 0.0, 1.0);
    CHECK(round(variance*1000)/1000 == expected);

}

TEST_CASE("[Cosine-power variance: mu=0, lambda=1, v=1]", "[statistics]") {
    double expected = 0.189;
    double v = 1.0;
    auto variance = cosine_power_mathematical_variance(v, 0.0, 1.0);
    CHECK(round(variance*1000)/1000 == expected);

}
TEST_CASE("[Cosine-power variance: mu=0, lambda=1, v=2]", "[statistics]") {
    double expected = 0.131;
    double v = 2.0;
    auto variance = cosine_power_mathematical_variance(v, 0.0, 1.0);
    CHECK(round(variance*1000)/1000 == expected);

}

TEST_CASE("[Cosine-power variance: mu=0, lambda=1, v=3]", "[statistics]") {
    double expected = 0.099;
    double v = 3.0;
    auto variance = cosine_power_mathematical_variance(v, 0.0, 1.0);
    CHECK(round(variance*1000)/1000 == expected);

}

TEST_CASE("[Cosine-power variance: mu=0, lambda=1, v=4]", "[statistics]") {
    double expected = 0.080;
    double v = 4.0;
    auto variance = cosine_power_mathematical_variance(v, 0.0, 1.0);
    CHECK(round(variance*1000)/1000 == expected);

}

TEST_CASE("[Cosine-power variance: mu=0, lambda=1, v=5]", "[statistics]") {
    double expected = 0.067;
    double v = 5.0;
    auto variance = cosine_power_mathematical_variance(v, 0.0, 1.0);
    CHECK(round(variance*1000)/1000 == expected);

}

TEST_CASE("[Cosine-power variance: mu=0, lambda=1, v=6]", "[statistics]") {
    double expected = 0.058;
    double v = 6.0;
    auto variance = cosine_power_mathematical_variance(v, 0.0, 1.0);
    CHECK(round(variance*1000)/1000 == expected);

}

//---------------------------------------------------------------------------------------------
//Excess

TEST_CASE("[Cosine-power excess: mu=0, lambda=1, v=0]", "[statistics]") {
    double expected = -1.200;
    double v = 0.0;
    auto variance = cosine_power_mathematical_excess(v, 0.0, 1.0);
    CHECK(round(variance*1000)/1000 == expected);

}

TEST_CASE("[Cosine-power excess: mu=0, lambda=1, v=1]", "[statistics]") {
    double expected = -0.806;
    double v = 1.0;
    auto variance = cosine_power_mathematical_excess(v, 0.0, 1.0);
    CHECK(round(variance*1000)/1000 == expected);

}


TEST_CASE("[Cosine-power excess: mu=0, lambda=1, v=2]", "[statistics]") {
    double expected = -0.594;
    double v = 2.0;
    auto variance = cosine_power_mathematical_excess(v, 0.0, 1.0);
    CHECK(round(variance*1000)/1000 == expected);

}


TEST_CASE("[Cosine-power excess: mu=0, lambda=1, v=3]", "[statistics]") {
    double expected = -0.466;
    double v = 3.0;
    auto variance = cosine_power_mathematical_excess(v, 0.0, 1.0);
    CHECK(round(variance*1000)/1000 == expected);

}


TEST_CASE("[Cosine-power excess: mu=0, lambda=1, v=4]", "[statistics]") {
    double expected = -0.381;
    double v = 4.0;
    auto variance = cosine_power_mathematical_excess(v, 0.0, 1.0);
    CHECK(round(variance*1000)/1000 == expected);

}


TEST_CASE("[Cosine-power excess: mu=0, lambda=1, v=5]", "[statistics]") {
    double expected = -0.322;
    double v = 5.0;
    auto variance = cosine_power_mathematical_excess(v, 0.0, 1.0);
    CHECK(round(variance*1000)/1000 == expected);

}


TEST_CASE("[Cosine-power excess: mu=0, lambda=1, v=6]", "[statistics]") {
    double expected = -0.278;
    double v = 6.0;
    auto variance = cosine_power_mathematical_excess(v, 0.0, 1.0);
    CHECK(round(variance*1000)/1000 == expected);

}


