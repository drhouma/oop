#pragma once

#include<cmath>
#include<stdexcept>

// Returns value of zeta function in point x (x > 1)
double zeta(double x);

// Return value of polygamma function of m-s order in point x (x = 1/2, 1, 3/2, 2, ...) if m = 0 and for positive x if m > 0
double polygamma(double x, int m);