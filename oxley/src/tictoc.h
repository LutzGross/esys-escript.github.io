
#ifndef _OXLEY_TICTOCCLOCK_
#define _OXLEY_TICTOCCLOCK_

#include <time.h>
#include <string>
#include <vector>

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

class incrementClock
{
	int numberOfIncrements;
	int currentIncrement;
	time_t starttime, currenttime, lasttime;
	std::vector<time_t> timer;

public:
	incrementClock(int n);
	~incrementClock();

	void toc(int n);
	void report();

};

#endif //_OXLEY_TICTOCCLOCK_