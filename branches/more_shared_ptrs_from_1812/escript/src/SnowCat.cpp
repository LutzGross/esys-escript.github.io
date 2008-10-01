
#include <iostream>

#include "SnowCat.h"

using namespace std;
using namespace escript;
using namespace boost;

namespace escript
{



SnowCat::SnowCat()
{
	v=0;
	cout << this << "Default SnowCat constructor called\n";
}

SnowCat::SnowCat(int value)
{
	v=value;
	cout << this << "Default SnowCat constructor called\n";
}


SnowCat::~SnowCat()
{
	cout << this << "SnowCat destructor called\n";
}

SnowCat::SnowCat(const SnowCat& other)
{
	v=other.v;
	cout << this << "Copy SnowCat constructor called\n";
}

SnowCat& SnowCat::operator=(const SnowCat& other)
{
	v=other.v;
	cout << this << "SnowCat::= called\n";
	return *this;
}

int SnowCat::getValue() const
{
	return v;
}

void SnowCat::setValue(int value)
{
	v=value;

}

SnowCat::shc 
SnowCat::makeShared() const
{
	return SnowCat::shc(new SnowCat());
}

SnowCat::shc
SnowCat::returnSelf(SnowCat::shc other) 
{
	return shared_from_this();
}

SnowCat::c2shc
SnowCat::returnSelfConst(SnowCat::shc other) const
{
	return shared_from_this();
}


SnowCat::shc
SnowCat::add(SnowCat::shc other, SnowCat::shc other2) const
{
	return shc(new SnowCat(other->v+other2->v));
}


SnowCat::cshc SnowCat::getConst1(SnowCat::shc other) const
{
// 	int h=other->v;
// 	(static_cast<SnowCat::cshc>(other))->setValue(-42);
// 	setValue(h);
	return other;		// or do we need to use their casting ops for that??

}

SnowCat::c2shc SnowCat::getConst2(SnowCat::shc other) const
{
// 	int h=other->v;
// 	(static_cast<SnowCat::c2shc>(other))->setValue(-42);
// 	setValue(h);
	return other;

} 






}
