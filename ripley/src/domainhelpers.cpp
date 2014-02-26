#include <ripley/domainhelpers.h>
#include <cmath>

void factorise(std::vector<int>& factors, int product) {
    int current = product;
    for (int p = 2; p <= sqrt(product); p++) {
        while (current % p == 0) {
            current /= p;
            factors.push_back(p);
        }
    }
    if (current != 1) {
	factors.push_back(current);
    }
}

