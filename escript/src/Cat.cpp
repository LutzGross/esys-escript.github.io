
#include <iostream>

#include "Cat.h"

using namespace std;
using namespace escript;
using namespace boost;

namespace escript
{



Cat::Cat()
{
	v=0;
	cout << this << "Default Cat constructor called\n";
}

Cat::Cat(int value)
{
	v=value;
	cout << this << "Default Cat constructor called\n";
}


Cat::~Cat()
{
	cout << this << "Cat destructor called\n";
}

Cat::Cat(const Cat& other)
{
	v=other.v;
	cout << this << "Copy Cat constructor called\n";
}

Cat& Cat::operator=(const Cat& other)
{
	v=other.v;
	cout << this << "Cat::= called\n";
	return *this;
}

int Cat::getValue() const
{
	return v;
}

void Cat::setValue(int value)
{
	v=value;

}

Cat::shc 
Cat::makeShared() const
{
	return Cat::shc(new Cat());
}

Cat::shc
Cat::returnSelf(Cat::shc other) const
{
	return other;
}

Cat::shc
Cat::add(Cat::shc other, Cat::shc other2) const
{
	return shc(new Cat(other->v+other2->v));
}


Cat::cshc Cat::getConst1(Cat::shc other) const
{
// 	int h=other->v;
// 	(static_cast<Cat::cshc>(other))->setValue(-42);
// 	setValue(h);
	return other;		// or do we need to use their casting ops for that??

}

Cat::c2shc Cat::getConst2(Cat::shc other) const
{
// 	int h=other->v;
// 	(static_cast<Cat::c2shc>(other))->setValue(-42);
// 	setValue(h);
	return other;

} 






}
