#include "../statistics.h"
#include "../catch.h"

using namespace std;

//-------------------------------------------------------------------------
//Density

TEST_CASE("[Shift-scale density: mu=1, lambda=2, v=0]", "[statistics]"){
    double expected = 0.25;
    auto cosine_func = cosine_power_density(0);
    auto func = shift_scaled_density(cosine_func, 1, 2);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Shift-scale density: mu=1, lambda=2, v=1]", "[statistics]"){
    double expected = 0.278;
    auto cosine_func = cosine_power_density(1);
    auto func = shift_scaled_density(cosine_func, 1, 2);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Shift-scale density: mu=5, lambda=2, v=1]", "[statistics]"){
    double expected = 0.0;
    auto cosine_func = cosine_power_density(1);
    auto func = shift_scaled_density(cosine_func, 5, 2);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Shift-scale density: mu=1, lambda=2, v=2]", "[statistics]"){
    double expected = 0.25;
    auto cosine_func = cosine_power_density(2);
    auto func = shift_scaled_density(cosine_func, 1, 2);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Shift-scale density: mu=1, lambda=2, v=3]", "[statistics]"){
    double expected = 0.208;
    auto cosine_func = cosine_power_density(3);
    auto func = shift_scaled_density(cosine_func, 1, 2);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Shift-scale density: mu=1, lambda=2, v=4]", "[statistics]"){
    double expected = 0.167;
    auto cosine_func = cosine_power_density(4);
    auto func = shift_scaled_density(cosine_func, 1, 2);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Shift-scale density: mu=1, lambda=2, v=5]", "[statistics]"){
    double expected = 0.130;
    auto cosine_func = cosine_power_density(5);
    auto func = shift_scaled_density(cosine_func, 1, 2);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Shift-scale density: mu=1, lambda=2, v=6]", "[statistics]"){
    double expected = 0.1;
    auto cosine_func = cosine_power_density(6);
    auto func = shift_scaled_density(cosine_func, 1, 2);
    CHECK(round(func(0)*1000)/1000 == expected);
}

//-----------------------------------------------------------------------------------------------
//Variance

TEST_CASE("[Shift scaled variance: mu=1, lambda=2, v=0]", "[statistics]") {
    double expected = 1.333;
    double v = 0.0;
    auto variance = cosine_power_mathematical_variance(v, 1.0, 2.0);
    CHECK(round(variance*1000)/1000 == expected);

}
TEST_CASE("[Shift scaled variance: mu=1, lambda=2, v=1]", "[statistics]") {
    double expected = 0.758;
    double v = 1.0;
    auto variance = cosine_power_mathematical_variance(1, 1.0, 2.0);
    CHECK(round(variance*1000)/1000 == expected);

}
TEST_CASE("[Shift scaled variance: mu=1, lambda=2, v=2]", "[statistics]") {
    double expected = 0.523;
    double v = 2.0;
    auto variance = cosine_power_mathematical_variance(v, 1.0, 2.0);
    CHECK(round(variance*1000)/1000 == expected);

}
TEST_CASE("[Shift scaled variance: mu=1, lambda=2, v=3]", "[statistics]") {
    double expected = 0.397;
    double v = 3.0;
    auto variance = cosine_power_mathematical_variance(v, 1.0, 2.0);
    CHECK(round(variance*1000)/1000 == expected);

}
TEST_CASE("[Shift scaled variance: mu=1, lambda=2, v=4]", "[statistics]") {
    double expected = 0.320;
    double v = 4.0;
    auto variance = cosine_power_mathematical_variance(v, 1.0, 2.0);
    CHECK(round(variance*1000)/1000 == expected);

}
TEST_CASE("[Shift scaled variance: mu=1, lambda=2, v=5]", "[statistics]") {
    double expected = 0.268;
    double v = 5.0;
    auto variance = cosine_power_mathematical_variance(v, 1.0, 2.0);
    CHECK(round(variance*1000)/1000 == expected);

}
TEST_CASE("[Shift scaled variance: mu=1, lambda=2, v=6]", "[statistics]") {
    double expected = 0.230;
    double v = 6.0;
    auto variance = cosine_power_mathematical_variance(v, 1.0, 2.0);
    CHECK(round(variance*1000)/1000 == expected);

}

//---------------------------------------------------------------------------------------------
//Excess

TEST_CASE("[Shift scaled excess: mu=1, lambda=2, v=0]", "[statistics]") {
    double expected = -1.200;
    double v = 0.0;
    auto variance = cosine_power_mathematical_excess(v, 1.0, 2.0);
    CHECK(round(variance*1000)/1000 == expected);

}

TEST_CASE("[Shift scaled excess: mu=1, lambda=2, v=1]", "[statistics]") {
    double expected = -0.806;
    double v = 1.0;
    auto variance = cosine_power_mathematical_excess(v, 1.0, 2.0);
    CHECK(round(variance*1000)/1000 == expected);

}


TEST_CASE("[Shift scaled excess: mu=1, lambda=2, v=2]", "[statistics]") {
    double expected = -0.594;
    double v = 2.0;
    auto variance = cosine_power_mathematical_excess(v, 1.0, 2.0);
    CHECK(round(variance*1000)/1000 == expected);

}


TEST_CASE("[Shift scaled excess: mu=1, lambda=2, v=3]", "[statistics]") {
    double expected = -0.466;
    double v = 3.0;
    auto variance = cosine_power_mathematical_excess(v, 1.0, 2.0);
    CHECK(round(variance*1000)/1000 == expected);

}


TEST_CASE("[Shift scaled excess: mu=1, lambda=2, v=4]", "[statistics]") {
    double expected = -0.381;
    double v = 4.0;
    auto variance = cosine_power_mathematical_excess(v, 1.0, 2.0);
    CHECK(round(variance*1000)/1000 == expected);

}


TEST_CASE("[Shift scaled excess: mu=1, lambda=2, v=5]", "[statistics]") {
    double expected = -0.322;
    double v = 5.0;
    auto variance = cosine_power_mathematical_excess(v, 1.0, 2.0);
    CHECK(round(variance*1000)/1000 == expected);

}


TEST_CASE("[Shift scaled excess: mu=1, lambda=2, v=6]", "[statistics]") {
    double expected = -0.278;
    double v = 6.0;
    auto variance = cosine_power_mathematical_excess(v, 1.0, 2.0);
    CHECK(round(variance*1000)/1000 == expected);

}