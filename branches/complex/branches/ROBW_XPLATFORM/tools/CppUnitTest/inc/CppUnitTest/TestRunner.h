#ifndef CPPUNIT_TEXTRUNNER_H
#define CPPUNIT_TEXTRUNNER_H

#include <iostream>
#include <vector>

#include "CppUnitTest/Test.h"

#include "CppUnitTest/CppUnitTestNamespace.h"
BEGIN_NAMESPACE_CPPUNITTEST

typedef std::pair<std::string, Test *>           mapping;
typedef std::vector<std::pair<std::string, Test *> >   mappings;

class TestRunner
{
protected:
    bool                                m_wait;
    std::vector<std::pair<std::string,Test *> >        m_mappings;

public:
	            TestRunner    () : m_wait (false) {}
                ~TestRunner   ();

    void        run           (int ac, char **av);
    void        addTest       (std::string name, Test *test)
    { m_mappings.push_back (mapping (name, test)); }

protected:
    void        run (Test *test);
    void        printBanner ();

};

END_NAMESPACE_CPPUNITTEST

#endif




