
#include <boost/smart_ptr.hpp>
#include <boost/enable_shared_from_this.hpp>

namespace escript
{


class SnowCat : public boost::enable_shared_from_this<SnowCat>
{
public:
	typedef boost::shared_ptr<SnowCat> shc;
	typedef const boost::shared_ptr<SnowCat> cshc;
	typedef boost::shared_ptr<const SnowCat> c2shc;
	SnowCat();
	SnowCat(int value);
	~SnowCat();
	SnowCat(const SnowCat& other);
	SnowCat& operator=(const SnowCat& other);
	int getValue() const;
	void setValue(int value);
	shc makeShared() const;
	shc returnSelf(shc other);
	c2shc returnSelfConst(shc other) const;
	shc add(shc other, shc other2) const;
	cshc getConst1(shc other) const;
	c2shc getConst2(shc other) const; 
private:
	int v;



};





}
