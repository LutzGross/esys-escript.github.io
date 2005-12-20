filter# $Id:$

#
#      COPYRIGHT ACcESS 2004 -  All Rights Reserved
#
#   This software is the property of ACcESS.  No part of this code
#   may be copied in any form or by any means without the expressed written
#   consent of ACcESS.  Copying, use or modification of this software
#   by any unauthorised person is illegal unless that
#   person has a software license agreement with ACcESS.
#

"""
A simple framework to run benchmarks under OPENMP and to summarize the results in tables for instance in HTML

@var __author__: name of author
@var __licence__: licence agreement
var __url__: url entry point on documentation
@var __version__: version
@var __date__: date of the version
"""

__author__="Lutz Gross, l.gross@uq.edu.au"
__licence__="contact: esys@access.uq.edu.au"
__url__="http://www.iservo.edu.au/esys/escript"
__version__="$Revision:$"
__date__="$Date:$"

import os,socket,time,sys

class BenchmarkSuite(object):
   """
   framework to run a bunch of L{Benchmark}s with the object to create a table of statistics.
   @var MAX_LEVEL: maximum number of level in headers for output
   """
   MAX_LEVEL=5
   def __init__(self,name=None):
       """
       sets up a suite of benchmarks
 
       @param name: name of the benchmark suite. If no name is given the class name is used.
       @type name: C{str}
       """
       super(BenchmarkSuite,self).__init__()
       self.__benchmarks=[]
       self.__scale=1
       if name==None:
           self.__name=self.__class__.__name__
       else:
          self.__name=name
          
   def __str__(self):
       """
       returns the name of the benchmark suite
       
       @return:  name
       @rtype: C{str}
       """
       return self.__name
   def addBenchmark(self,benchmark):
       """
       adds a L{Benchmark} to the suite

       @param benchmark: adds a new L{Benchmark} to the suite
       @type benchmark: L{Benchmark}
       """
       self.__benchmarks.append(benchmark)        
   def __len__(self):
       """
       returns the number of benchmarks in the suite

       @return:  number of benchmarks 
       @rtype: C{int}       
       """
       return len(self.__benchmarks)
   def __getitem__(self,i):
       """
       returns the i-th benchmark in the suite through self[i]
       
       @param i: index of the requested benchmark
       @type i: C{int}
       @return:  i-th benchmark 
       @rtype: L{Benchmark}       

       """
       return self.__benchmarks[i]
   def run(self,scale=1):
       """
       runs all benchmarks

       @param scale: defines the number of (OpenMP) threads to be used. If scale is a scalar all benchmarks 
                     are run with scale number of threads. If scale is a C{list}, the p-th problem in each of the benchmarks
                     in the suite is run with scale[p] threads. If scale[p]<1 teh p-th problem is omitted.
       @type scale: C{int} or C{list} of C{int}s. 
       """
       self.__scale=scale       
       for i in range(len(self)): self[i].run(scale=scale)
   def getHTML(self,filter,level=1):
       """
       returns the results of the last benchmark run in HTML format.

       @param filter: filter to be applied to the results
       @type filter: L{BenchmarkFilter}
       @param level: level used in header <H?> tags 
       @type level: C{int}
       @return: HTML document
       @rtype: C{str}
       """
       out=""
       if level==1: out+="<HTML><HEAD><TITLE>Benchmark: %s</TITLE></HEAD><BODY>\n"%str(self)
       out+="<H%s>%s</H%s>\n"%(level,str(self),level)
       if level==1: 
           m=""
           if isinstance(self.__scale,int):
              if self.__scale>1:
                  m=" (%s threads)"%self.__scale
           out+="<p>platform: %s%s</p>\n"%(socket.gethostname(),m)
       for i in range(len(self)):
           out+="<p>\n"
           out+=self[i].getHTML(filter=filter,level=min(level+1,self.MAX_LEVEL))
           out+="<p>\n"
       if level==1:
           out+="<hr><p align=\"center\">by %s at %s</p>\n"%(os.getlogin(),time.strftime('%X %x %Z'))
           out+="</BODY></HTML>\n"
       return out


class Benchmark(object):
   """
   runs a bunch of similar L{BenchmarkProblem}s with a bunch of L{Options}
   """
   def __init__(self,name=None,description=None):
       """
       sets up a benchmark 
 
       @param name: name of the benchmark. If no name is given the class name is used.
       @type name: C{str}
       @param description: description of the benchmark. 
       @type description: C{str} or C{None}      
       """
       super(Benchmark,self).__init__()
       self.__options=[]
       self.__problems=[]
       self.__results=[]
       self.__scale=1
       if name==None:
           self.__name=self.__class__.__name__
       else:
          self.__name=name
       self.__description=description
       
   def __str__(self):
       """
       returns the name of the benchmark suite
       
       @return:  name
       @rtype: C{str}
       """
       return self.__name
          
   def addProblem(self,problem):
       """
       adds a problem to the benchmark

       @param problem: adds a new problem to the bechmark
       @type problem: L{BenchmarkProblem}
       """
       self.__problems.append(problem)

   def addOptions(self,Options):
       """
       adds a options to the benchmark

       @param options: adds a new option to the bechmark
       @type problem: L{Options}
       """
       self.__options.append(Options)

   def run(self,scale=1):
       """
       runs all problems with all options. 


       @param scale: defines the number of (OpenMP) threads to be used. If scale is a scalar all benchmarks 
                     are run with scale number of threads. If scale is a C{list}, the p-th problem in each of the benchmarks
                     in the suite is run with scale[p] threads. If scale[p]<1 teh p-th problem is omitted.
       @type scale: C{int} or C{list} of C{int}s. 
       """
       if isinstance(scale,list):
           c_max=min(len(scale),len(self.__problems))
       else:
           c_max=len(self.__problems)
       self.__filter=filter
       self.__scale=scale
       self.__results=[]
       for c in range(c_max):
          r=self.__problems[c]
          if isinstance(scale,list):
             s=scale[c]
          else:
             s=scale
          row=[]
          if s>0:
              for p in self.__options:
                  os.environ['OMP_NUM_THREADS']=str(s)
                  row.append(r.run(p))
          self.__results.append(row)
   def getHTML(self,filter,level=1):
       """
       returns the results of the last benchmark run in HTML format.

       @param filter: filter to be applied to the results
       @type filter: L{BenchmarkFilter}
       @param level: level used in header <H?> tags 
       @type level: C{int}
       @return: HTML document
       @rtype: C{str}
       """
       out=""
       if level==1: out+="<HTML><HEAD><TITLE>Benchmark: %s</TITLE></HEAD><BODY>\n"%str(self)
       out+="<H%s>%s</H%s>\n"%(level,str(self),level)
       if level==1: 
         m=""
         if isinstance(self.__scale,int):
            if self.__scale>1:
                m=" (%s threads)"%self.__scale
         out+="<p>platform: %s%s</p>\n"%(socket.gethostname(),m)
       if self.__description: out+="<p>%s</p>\n"%str(self.__description)   
       if len(self.__problems)>0:
          out+="<TABLE ALIGN=\"center\" BORDER=3 CELLPADDING=5 CELLSPACING=1>\n"
          h1_seg=""
          rn=filter.getResultNames()
          if len(rn)==0:
             h1_seg+="<TD></TD>"
          else:
             for n in rn: h1_seg+="<TD ALIGN=\"center\">%s</TD>"%n        
          h0="<TR><TH ALIGN=\"center\" ROWSPAN=2>Case</TH>"
          h1="<TR>"
          if isinstance(self.__scale,list): h0+="<TH ALIGN=\"center\" ROWSPAN=2>Threads</TH>" 
          for o in self.__options:
                 if len(rn)==0:
                     h0+="<TH ALIGN=\"center\">%s</TH>"%str(o)
                 elif len(rn)==1:
                     h0+="<TH ALIGN=\"center\">%s</TH>"%str(o)
                     empty_h1=False
                 else:
                     h0+="<TH ALIGN=\"center\" COLSPAN=%s>%s</TH>"%(len(rn),str(o))
                 h1+=h1_seg
          out+=h0+"</TR>\n"+h1+"</TR>\n"
          c=0
          for r in range(len(self.__results)):
             out+="<TR><TH ALIGN=\"right\">%s</TH>"%str(self.__problems[r])
             if isinstance(self.__scale,list): out+="<TD ALIGN=\"right\">%s</TD>"%self.__scale[c] 
             for col in self.__results[r]:
                   for e in filter(col): out+="<TD ALIGN=\"right\">%s</TD>"%e
             out+="</TR>\n"
             c+=1
          out+="</TABLE>"
       if level==1:
          out+="<hr><p align=\"center\">by %s at %s</p>\n"%(os.getlogin(),time.strftime('%X %x %Z'))
          out+="</BODY></HTML>\n"
       return out 
       
class BenchmarkProblem(object):
   """
   something that can be run and returns a list of characteristics such as timing, Mflops, error, etc.
   """
   def __init__(self,name=None):
       """
       sets up a benchmark problem
 
       @param name: name of the problem. If no name is given the class name is used.
       @type name: C{str}
       """
       super(BenchmarkProblem,self).__init__()
       if name==None:
           self.__name=self.__class__.__name__
       else:
          self.__name=name

       
   def __str__(self):
       """
       returns the name of the benchmark suite
       
       @return:  name
       @rtype: C{str}
       """
       return self.__name

   def run(self,options=None):
       """
       runs the problem and returns a list of run characteristics


       @param options: the options that are used for the run. Note that the number of OpenMP threads is controlled 
                       by the L{Benchmark} the problem is run in.
       @type options: L{Options}
       @return: run characteristics
       @rtype: any type that can be read by the L{BenchmarkFilter} applied to it.
       @remark: this function has to overwritten by a particular problem
       """
       raise NotImplementedError
       return []
    
class BenchmarkFilter(object):
   """
   object to filter the characteristcs returned by Bechmark runs.
   
   """
   def __init__(self):
       """
       sets up a filter
       """
       pass
 

   def getResultNames(self):
       """
       return the names of the results produced when run() is called.
       
       @return: names the list of the names to be used when the results of the run() call are printed
       @rtype: C{list} of C{str}
       @remark: this function has to overwritten by a particular problem
       """
       raise NotImplementedError
       return []

   def __call__(self,result):
       """
       filters out values results returned as characteristcs of a problem run
       
       @param result: values to be filtered
       @type result: any type that is produced by the L{BenchmarkProblem} it is applied to
       @return: a list of strings selected from result
       @rtype: C{list} of C{str}
       @remark: this function has to overwritten by a particular problem
       """
       raise NotImplementedError
       return []


class Options(object):
    """
    defines a set of options to be used to run a L{BenchmarkProblem} 
    """
    def __init__(self,name=None):
       """
       sets up the options
 
       @param name: name of the option. If no name is given the class name is used.
       @type name: C{str}
       """
       super(Options,self).__init__()
       if name==None:
          self.__name=self.__class__.__name__
       else:
          self.__name=name
    def __str__(self):
       """
       returns the name of the benchmark suite
       
       @return:  name
       @rtype: C{str}
       """
       return self.__name
    
if __name__=="__main__":

    class OptionsTest1(Options):
        pass
    class OptionsTest2(Options):
        pass

    class BenchmarkProblemTest1(BenchmarkProblem):
       def __init__(self):
           super(BenchmarkProblemTest1,self).__init__(name="TEST1")
       def run(self,options=None):
           import time
           return time.time(),"A"

    class BenchmarkProblemTest2(BenchmarkProblem):
       def __init__(self):
           super(BenchmarkProblemTest2,self).__init__(name="TEST2")
       def run(self,options=None):
           import time
           return -time.time(),"B"

    class SimpleFilter(BenchmarkFilter):
       def getResultNames(self):
            return ["r0","r1"]  
       def __call__(self,result):
            return [str(result[0]),str(result[1])]

    bm=Benchmark("Example")
    bm.addProblem(BenchmarkProblemTest1())
    bm.addProblem(BenchmarkProblemTest2())
    bm.addOptions(OptionsTest1())
    bm.addOptions(OptionsTest2())

    bms=BenchmarkSuite("A Test")
    bms.addBenchmark(bm)

    bms.run()
    print bms.getHTML(filter=SimpleFilter())
    
    bms.run(scale=4)
    print bms.getHTML(filter=SimpleFilter())

    bms.run(scale=[1,2])
    print bms.getHTML(filter=SimpleFilter())
