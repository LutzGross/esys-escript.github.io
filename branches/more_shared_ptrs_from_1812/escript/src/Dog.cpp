
#include <iostream>

#include "Dog.h"

using namespace std;
using namespace escript;
using namespace boost;

namespace escript
{



Dog::Dog()
{
	v=0;
	refcount=0;
	cout << this << "Default Dog constructor called\n";
}

Dog::Dog(int value)
{
	refcount=0;
	v=value;
	cout << this << "Default Dog constructor called\n";
}


Dog::~Dog()
{
	cout << this << "Dog destructor called\n";
}

Dog::Dog(const Dog& other)
{
	v=other.v;
	refcount=0;
	cout << this << "Copy Dog constructor called\n";
}

Dog& Dog::operator=(const Dog& other)
{
	v=other.v;
	cout << this << "Dog::= called\n";
	return *this;
}

int Dog::getValue() const
{
	return v;
}

void Dog::setValue(int value)
{
	v=value;

}

Dog::shc 
Dog::makeShared() const
{
	return Dog::shc(new Dog());
}

Dog::shc
Dog::returnSelf(Dog::shc other) const
{
	return other;
}

Dog::shc
Dog::add(Dog::shc other, Dog::shc other2) const
{
	return shc(new Dog(other->v+other2->v));
}


Dog::cshc Dog::getConst1(Dog::shc other) const
{
	// some other testing
	c2shc p(other);


	return other;		// or do we need to use their casting ops for that??

}

Dog::c2shc Dog::getConst2(Dog::shc other) const
{
// 	int h=other->v;
// 	(static_cast<Dog::c2shc>(other))->setValue(-42);
// 	setValue(h);
	return other;

} 

void Dog::addRef() const{
	refcount++;
cout << "Refcount now=" << refcount << endl;
}

int Dog::decRef() const
{
	return --refcount;
}






}

namespace boost
{

// void intrusive_ptr_add_ref(escript::Dog* d)
// {
// cout << d << "Increasing ref on dog\n";
// 	d->addRef();
// }
// 
// void intrusive_ptr_release(escript::Dog* d)
// {
// int z=d->decRef();
// cout << d << "Decreasing ref -> " << z << endl;
// 	if (z==0)
// 	{
// cout << d << "About to delete\n";
// 		delete d;
// 	}
// }


void intrusive_ptr_add_ref(const escript::Dog* d)
{
cout << d << "Increasing ref on dog\n";
	d->addRef();
}

void intrusive_ptr_release(const escript::Dog* d)
{
int z=d->decRef();
cout << d << "Decreasing ref -> " << z << endl;
	if (z==0)
	{
cout << d << "About to delete\n";
		delete d;
	}
}


}
