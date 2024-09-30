/*
// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER
*/

/**
 * This file shows examples of indentatation of the "thyra" emacs style that
 *  is ment to follow the style guidelines in "Code Complete", 2nd edition.
 */


//
// Prototypes
//


//
// Don't indent for namespaces enclosures
//

namespace NamespaceA {


// Indent arguments on continuation lines one offset from from the beginning
// of the function type+name line.  This style does deviate from the
// recommended style in "Code Complete" in that the paren lines up with the
// argument list, not the type+name line.
void func1( int a, int b, int c,
  int d, int e, int f,
  int g, int h, int i
  );


// Same as above, except the the argument list starts on the line below the
// opening paren.  This case can be handled differently in emacs.
void func2(
  int a, int b, int c,
  int d, int e, int f,
  int g, int h, int i
  );


} // namespace NamespaceA


//
// Defintions:
//


// The following function definitions shows a few things:
// 1) The defintions are indented two spaces from other entities
// 2) The function begin '{' and end '}' are both indented from
//    the rest of the code by one space.  This sets off the
//    boundaries for the function.

void NamespaceA::func1( int a, int b, int c,
  int d, int e, int f,
  int g, int h, int i
  )
{
  
  // Indent continuation lines on variable declarations one offset.
  double aa, bb, cc,
    dd;
  
  {
    std::vector<double> va(a);

    // Use "pure block emulation" for one-line control statements

    for ( int i = 0; i < a; ++i ) {
      if ( i*a < b ) {
        va[i] = 2.0;
      }
      else if ( i*b < c ) {
        va[i] = 2.5;
      }
      else {
        va[i] = 3.0;
      }
    }

    // Uses "unindented begin-end pairs" (not recommended, but see below).

    for ( int i = 0; i < a; ++i )
    {
      if ( i*a < b )
      {
        va[i] = 2.0;
      }
      else if ( i*b < c )
      {
        va[i] = 2.5;
      }
      else
      {
        va[i] = 3.0;
      }
    }

    // Above, not that (x)emacs shows the match for the opening '{' plus the
    // line above it when the '{' is not in the screen!
    
    // Indent case labels within switch statements

    switch(d) {
      case 0:
        aa = 4.0;
        break;
      case 1:
        aa = 5.0;
        break;
      case 2:
        aa = 6.0;
        break;
      default:
        TEUCHOS_TEST_FOR_EXCEPT("Should never get here!");
    }

    // For control statements that extend over one line, use "unindented
    // begin-end pairs".  This breaks with the advise of "Code Complete", 2nd
    // edition but this is more standard within the C++ community than the
    // recommended indented "begin/end pairs".  To be consistent with the
    // initial 'if' statement block within the same if/else if/else structure,
    // I put the '{' on the next line from the 'else' statement.

    if(
      a < b
      && c > d
      && f < g
      )
    {
      bb = 8.0;
    }
    else if( h < i ) {
      bb = 9.0;
    }
    else
    {
      cc = 10.0;
    }
    
  }
  
}


// Indented two spaces from above end of function '}'.
void NamespaceA::func2(
  int a, int b, int c,
  int d, int e, int f,
  int g, int h, int i
  )
{

  // The function arguments on continuation lines in a function call are
  // indented one offset instead of aligning them with the opening '('.
  func1( a, b, c, d, e,
    f, g, h, i );

  // Same as above, except that the arguments start on the next line from the
  // opening '('.  This is handled differently in emacs.  Also, note that the
  // closing ')' is aligned with the arguments and not the 'func2(' beginning
  // 
  func2(
    a, b, c, d, e,
    f, g, h, i
    );

}

