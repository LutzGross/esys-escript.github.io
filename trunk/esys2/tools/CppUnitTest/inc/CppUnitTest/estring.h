
#ifndef CPPUNIT_ESTRING_H
#define CPPUNIT_ESTRING_H

#include <cstdio>

#include "CppUnitTest/CppUnitTestNamespace.h"
BEGIN_NAMESPACE_CPPUNITTEST

  // Create a string from a const char pointer
  inline std::string estring (const char *cstring)
    { return std::string (cstring); }

  // Create a string from a string (for uniformities' sake)
  inline std::string estring (std::string& expandedString)
    { return expandedString; }

  // Create a string from an int
  inline std::string estring (int number)
    { char buffer [50]; sprintf (buffer, "%d", number); return std::string (buffer); }

  // Create a string from a long
  inline std::string estring (long number)
    { char buffer [50]; sprintf (buffer, "%ld", number); return std::string (buffer); }

  // Create a string from a double
  inline std::string estring (double number)
    { char buffer [50]; sprintf (buffer, "%f", number); return std::string (buffer); }

END_NAMESPACE_CPPUNITTEST

#endif




