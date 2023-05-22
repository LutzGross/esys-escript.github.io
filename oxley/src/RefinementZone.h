
#include <vector>

#include <oxley/RefinementType.h>

class RefinementZone
{
public:

private:

	/**
       \brief
       A queue of refinements 
    */
	std::vector<RefinementAlgorithm> queue;


	/**
       \brief
       Add to queue
    */
	void AddToQueue(RefinementAlgorithm R);

// protected:
};