
#ifndef CPPUNITTEST_NAMESPACE_H
#define CPPUNITTEST_NAMESPACE_H
#
//
// Assume all compilers support namespaces. Easy to redefine the defines
// if needed to support compilers that don't support namespaces.
//
#define COMPILER_SUPPORTS_NAMESPACES

#if defined COMPILER_SUPPORTS_NAMESPACES
#  define BEGIN_NAMESPACE_CPPUNITTEST namespace CppUnitTest {
#  define END_NAMESPACE_CPPUNITTEST }
#  define USING_NAMESPACE_CPPUNITTEST using namespace CppUnitTest ;

#else
//
// Empty defines for pre-namespace compilers.
//
#  define BEGIN_NAMESPACE_CPPUNITTEST
#  define END_NAMESPACE_CPPUNITTEST
#  define USING_NAMESPACE_CPPUNITTEST

# endif

#endif

