#include "esys/SerialAccessor.h"

using namespace std;
using namespace boost::python;

SerialAccessor::SerialAccessor(const std::vector<int>& pointDataShape, 
			       int blockSize, int numBlocks):
  m_shape(pointDataShape),
  m_blockSize(blockSize),
  m_numBlocks(numBlocks)
{
}

SerialAccessor::~SerialAccessor()
{
}

void SerialAccessor::copy(const boost::python::numeric::array& pointData,ValueContainerType& values) const
{
  vector<int> shape=extract<vector<int> > (pointData.getshape());
  if (m_shape!=shape) {
    AccessorException e("Couldn't copy. Input data shape doesn't match accessor shape.");
    throw e;
  }
}

void SerialAccessor::evaluate(ValueContainerType& values, void (*operation)(double)) const
{
}
