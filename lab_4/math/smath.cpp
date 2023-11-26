#include "smath.h"

double EPS = 0.000000000001;

double zeta(double x) {
	if (x <= 1) {
		throw std::runtime_error("Parameter x must be higher then 1");
	}
	double output = 0, part = 1;
	for (int k = 1; abs(part) > EPS; ++k) {
		part = 1.0 / pow(k, x);
		output += part;
	}
	return output;
}

double polygamma(double x, int m) {
	if (m < 0) {
		throw std::invalid_argument("Parameter m must be integer >= 0");
	}
	if (x <= 0) {
		throw std::invalid_argument("Parameter x must be > 0");
	}
	double output = 0;
	if (m == 0 || abs(round(x*2)-x*2) < EPS) {
		double temp = x;
		for (; abs(temp - 1) > EPS && abs(temp - 0.5) > EPS && temp > 0; --temp) continue;
		if (temp < 0) {
			throw std::invalid_argument("Parameter x must be number in format k/2 where k is positive integer");
		}
		
		for (; abs(x - 1) > EPS && abs(x - 0.5) > EPS; --x) {
			output += pow(-1, m) * tgamma(m + 1) / double(pow(x - 1, m + 1));
		}
		if (abs(x - 1) <= EPS) {
			if (m > 0) {
				return output + pow(-1, m + 1) * tgamma(m + 1) * zeta(m + 1);
			}
			else if (m == 0) {
				return output - 0.57721566490153286060;
			}
		}
		if (abs(x - 0.5) <= EPS) {
			if (m > 0) {
				return output + pow(-1, m + 1) * tgamma(m + 1) * (pow(2, m + 1) - 1) * zeta(m + 1);
			}
			else if (m == 0) {
				return output - 0.57721566490153286060 - 2 * log(2.0);
			}
		}
	}
	else {
		double part = 1;
		for (int k = 0; abs(part) > EPS; ++k) {
			part = pow(-1, m + 1) * tgamma(m + 1) / pow(x + k, m + 1);
			output += part;
		}
		return output;
	}
}