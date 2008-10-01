
#include <boost/smart_ptr.hpp>

namespace escript
{


class Cat
{
public:
	typedef boost::shared_ptr<Cat> shc;
	typedef const boost::shared_ptr<Cat> cshc;
	typedef boost::shared_ptr<const Cat> c2shc;
	Cat();
	Cat(int value);
	~Cat();
	Cat(const Cat& other);
	Cat& operator=(const Cat& other);
	int getValue() const;
	void setValue(int value);
	shc makeShared() const;
	shc returnSelf(shc other) const;
	shc add(shc other, shc other2) const;
	cshc getConst1(shc other) const;
	c2shc getConst2(shc other) const; 
private:
	int v;



};





}
