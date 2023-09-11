#include "../statistics.h"
#include "../catch.h"

using namespace std;

//-------------------------------------------------------------------------
//Density

TEST_CASE("[Shift-scale density: mu=3, lambda=2, v=0]", "[statistics]"){
    double expected = 0.25;
    auto cosine_func = cosine_power_density(0);
    auto func = shift_scaled_density(cosine_func, 3, 2);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Shift-scale density: mu=3, lambda=2, v=1]", "[statistics]"){
    double expected = -0.278;
    auto cosine_func = cosine_power_density(1);
    auto func = shift_scaled_density(cosine_func, 3, 2);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Shift-scale density: mu=3, lambda=2, v=2]", "[statistics]"){
    double expected = 0.25;
    auto cosine_func = cosine_power_density(2);
    auto func = shift_scaled_density(cosine_func, 3, 2);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Shift-scale density: mu=3, lambda=2, v=3]", "[statistics]"){
    double expected = -0.208;
    auto cosine_func = cosine_power_density(3);
    auto func = shift_scaled_density(cosine_func, 3, 2);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Shift-scale density: mu=3, lambda=2, v=4]", "[statistics]"){
    double expected = 0.167;
    auto cosine_func = cosine_power_density(4);
    auto func = shift_scaled_density(cosine_func, 3, 2);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Shift-scale density: mu=3, lambda=2, v=5]", "[statistics]"){
    double expected = -0.130;
    auto cosine_func = cosine_power_density(5);
    auto func = shift_scaled_density(cosine_func, 3, 2);
    CHECK(round(func(0)*1000)/1000 == expected);
}

TEST_CASE("[Shift-scale density: mu=3, lambda=2, v=6]", "[statistics]"){
    double expected = 0.1;
    auto cosine_func = cosine_power_density(6);
    auto func = shift_scaled_density(cosine_func, 3, 2);
    CHECK(round(func(0)*1000)/1000 == expected);
}