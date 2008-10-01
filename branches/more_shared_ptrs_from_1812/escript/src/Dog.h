

#ifndef DOG_H
#define DOG_H

#include <boost/smart_ptr.hpp>


namespace escript
{


class Dog
{
public:
	typedef boost::intrusive_ptr<Dog> shc;
	typedef const boost::intrusive_ptr<Dog> cshc;
	typedef boost::intrusive_ptr<const Dog> c2shc;
	Dog();
	Dog(int value);
	~Dog();
	Dog(const Dog& other);
	Dog& operator=(const Dog& other);
	int getValue() const;
	void setValue(int value);
	shc makeShared() const;
	shc returnSelf(shc other) const;
	shc add(shc other, shc other2) const;
	cshc getConst1(shc other) const;

	// This does not work because python does not have a converter for the const form
	c2shc getConst2(shc other) const; 

	void addRef() const;
	int decRef() const;
private:
	int v;
	mutable int refcount;



};



}



namespace boost{

// void intrusive_ptr_add_ref(escript::Dog* d);
// void intrusive_ptr_release(escript::Dog* d);

void intrusive_ptr_add_ref(const escript::Dog* d);
void intrusive_ptr_release(const escript::Dog* d);


}



#endif

