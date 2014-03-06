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

void doublyLink(std::vector<ripley::IndexVector>& va,
        std::vector<ripley::IndexVector>& vb, int a, int b) {
    va[a].push_back(b);
    vb[b].push_back(a);
}
