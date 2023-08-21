
#ifndef _OXLEY_TICTOCCLOCK_
#define _OXLEY_TICTOCCLOCK_

#include <time.h>
#include <string>

class TicTocClock
{
	time_t starttime;

public:
	TicTocClock();
	~TicTocClock();

	void toc(std::string message);
	void setTime(time_t newtime);
	time_t getTime();
};

#endif //_OXLEY_TICTOCCLOCK_