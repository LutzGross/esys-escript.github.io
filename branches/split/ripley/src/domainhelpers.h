#ifndef _DOMAINHELPERS_H_
#define _DOMAINHELPERS_H_

#include <vector>
#include <ripley/Ripley.h>

void factorise(std::vector<int>& factors, int product);

void doublyLink(std::vector<ripley::IndexVector>& va,
        std::vector<ripley::IndexVector>& vb, int a, int b);
#endif
