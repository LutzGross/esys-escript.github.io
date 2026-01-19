
##############################################################################
#
# Copyright (c) 2003-2020 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
# Development from 2019 by School of Earth and Environmental Sciences
#
##############################################################################

"""
A simple framework to run benchmarks under OpenMP and to summarize the results
in tables for instance in HTML.

:var __author__: name of author
:var __license__: licence agreement
:var __copyright__: copyrights
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""


__copyright__="""Copyright (c) 2003-2020 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"
__author__="Lutz Gross, l.gross@uq.edu.au"

import os, socket, sys, time, traceback
from . import escriptcpp as esc

class BenchmarkSuite(object):
   """
   Framework to run a bunch of `Benchmark` s using the object and creating a
   table of statistics.

   :cvar MAX_LEVEL: maximum number of level in headers for output
   """
   MAX_LEVEL=5
   def __init__(self,name=None):
       """
       Sets up a suite of benchmarks.

       :param name: name of the benchmark suite. If no name is given the class
                    name is used.
       :type name: ``str``
       """
       super(BenchmarkSuite,self).__init__()
       self.__benchmarks=[]
       self.__scale=1
       if name is None:
           self.__name=self.__class__.__name__
       else:
          self.__name=name

   def __str__(self):
       """
       Returns the name of the benchmark suite.

       :return: the name
       :rtype: ``str``
       """
       return self.__name

   def addBenchmark(self,benchmark):
       """
       Adds a new `Benchmark` to the suite.

       :param benchmark: the benchmark to add
       :type benchmark: `Benchmark`
       """
       self.__benchmarks.append(benchmark)

   def __len__(self):
       """
       Returns the number of benchmarks in the suite.

       :return: number of benchmarks
       :rtype: ``int``
       """
       return len(self.__benchmarks)

   def __getitem__(self,i):
       """
       Returns the i-th benchmark in the suite through self[i].

       :param i: index of the requested benchmark
       :type i: ``int``
       :return: i-th benchmark
       :rtype: `Benchmark`

       """
       return self.__benchmarks[i]

   def run(self,scale=1):
       """
       Runs all benchmarks.

       :param scale: defines the number of (OpenMP) threads to be used. If
                     ``scale`` is a scalar all benchmarks are run with ``scale``
                     number of threads. If ``scale`` is a ``list``, the p-th
                     problem in each of the benchmarks in the suite is run with
                     ``scale[p]`` threads. If ``scale[p]`` <1 the p-th problem is
                     omitted.
       :type scale: ``int`` or ``list`` of ``int``
       """
       self.__scale=scale
       for i in range(len(self)): self[i].run(scale=scale)

   def getHTML(self,filter,level=1):
       """
       Returns the results of the last benchmark run in HTML format.

       :param filter: filter to be applied to the results
       :type filter: `BenchmarkFilter`
       :param level: level used in header <H?> tags
       :type level: ``int``
       :return: HTML document
       :rtype: ``str``
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
           try:
               name=os.getlogin()
               out+="<hr><p align=\"center\">by %s at %s</p>\n"%(name,time.strftime('%X %x %Z'))
           except OSError:
               out+="<hr><p align=\"center\">%s</p>\n"%(time.strftime('%X %x %Z'))

           out+="</BODY></HTML>\n"
       return out


class Benchmark(object):
   """
   Runs a bunch of similar `BenchmarkProblem` s with a bunch of `Options`.
   """
   def __init__(self,name=None,description=None):
       """
       Sets up a benchmark.

       :param name: name of the benchmark. If no name is given the class name
                    is used.
       :type name: ``str``
       :param description: description of the benchmark
       :type description: ``str`` or ``None``
       """
       super(Benchmark,self).__init__()
       self.__options=[]
       self.__problems=[]
       self.__results=[]
       self.__scale=1
       if name is None:
           self.__name=self.__class__.__name__
       else:
          self.__name=name
       self.__description=description

   def __str__(self):
       """
       Returns the name of the benchmark suite.

       :return: the name
       :rtype: ``str``
       """
       return self.__name

   def addProblem(self,problem):
       """
       Adds a problem to the benchmark.

       :param problem: the problem to be added
       :type problem: `BenchmarkProblem`
       """
       self.__problems.append(problem)

   def addOptions(self,options):
       """
       Adds options to the benchmark.

       :param options: the options to be added to the benchmark. If
                       options is None the options are left unchanged.
       :type options: `Options`
       """
       if options!=None: self.__options.append(options)

   def run(self,scale=1):
       """
       Runs all problems with all options.

       :param scale: defines the number of (OpenMP) threads to be used. If
                     ``scale`` is a scalar all benchmarks are run with ``scale``
                     number of threads. If ``scale`` is a ``list`` , the p-th
                     problem in each of the benchmarks in the suite is run with
                     ``scale[p]`` threads. If ``scale[p]`` <1 the p-th problem is
                     omitted.
       :type scale: ``int`` or ``list`` of ``int`` s
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
              t0=time.time()
              print(("%s with %s threads started."%(r.__class__,s)))
              for p in self.__options:
                  esc.setNumberOfThreads(s)
                  try:
                     row.append(r.run(p))
                  except:
                     traceback.print_exc(file=sys.stdout)
                     row.append(None)
              t0=time.time()-t0
              print(("%s with %s threads finished (walltime=%s sec)."%(r.__class__,s,t0)))
          self.__results.append(row)

   def getHTML(self,filter,level=1):
       """
       Returns the results of the last benchmark run in HTML format.

       :param filter: filter to be applied to the results
       :type filter: `BenchmarkFilter`
       :param level: level used in header <H?> tags
       :type level: ``int``
       :return: HTML document
       :rtype: ``str``
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
                     colspan=1
                 elif len(rn)==1:
                     h0+="<TH ALIGN=\"center\">%s</TH>"%str(o)
                     colspan=1
                     empty_h1=False
                 else:
                     colspan=len(rn)
                     h0+="<TH ALIGN=\"center\" COLSPAN=%s>%s</TH>"%(colspan,str(o))
                 h1+=h1_seg
          out+=h0+"</TR>\n"+h1+"</TR>\n"
          c=0
          for r in range(len(self.__results)):
             out+="<TR><TH ALIGN=\"right\">%s</TH>"%str(self.__problems[r])
             if isinstance(self.__scale,list):
                 out+="<TD ALIGN=\"right\">%s</TD>"%self.__scale[c]
             for col in self.__results[r]:
                   if col is None:
                      out+="<TD ALIGN=\"center\" COLSPAN=%s>failed.</TD>"%colspan
                   else:
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
   Represents a benchmark problem that can be run and which returns a list of
   characteristics such as timing, MFlops, error, etc.
   """
   def __init__(self,name=None):
       """
       Sets up a benchmark problem.

       :param name: name of the problem. If no name is given the class name
                    is used.
       :type name: ``str``
       """
       super(BenchmarkProblem,self).__init__()
       if name is None:
           self.__name=self.__class__.__name__
       else:
          self.__name=name

   def __str__(self):
       """
       Returns the name of the benchmark suite.

       :return: the name
       :rtype: ``str``
       """
       return self.__name

   def run(self,options=None):
       """
       Runs the problem and returns a list of run characteristics.

       :param options: the options that are used for the run. Note that the
                       number of OpenMP threads is controlled by the
                       `Benchmark` the problem is run in.
       :type options: `Options`
       :return: run characteristics
       :rtype: any type that can be read by the `BenchmarkFilter` applied
               to it
       :note: this function has to be overwritten by a particular problem
       """
       raise NotImplementedError
       return []

class BenchmarkFilter(object):
   """
   Object to filter the characteristics returned by Benchmark runs.

   """
   def __init__(self):
       """
       Sets up a filter.
       """
       pass

   def getResultNames(self):
       """
       Returns the names of the results produced when ``run()`` is called.

       :return: the list of the names to be used when the results of
                the ``run()`` call are printed
       :rtype: ``list`` of ``str``
       :note: this function has to overwritten by a particular problem
       """
       raise NotImplementedError
       return []

   def __call__(self,result):
       """
       Filters out results returned as characteristics of a problem run.

       :param result: values to be filtered
       :type result: any type that is produced by the `BenchmarkProblem`
                     it is applied to
       :return: a list of strings selected from result
       :rtype: ``list`` of ``str``
       :note: this function has to be overwritten by a particular problem
       """
       raise NotImplementedError
       return []


class Options(object):
    """
    Defines a set of options to be used to run a `BenchmarkProblem`.
    """
    def __init__(self,name=None):
       """
       Sets up the options.

       :param name: name of the option. If no name is given the class name
                    is used.
       :type name: ``str``
       """
       super(Options,self).__init__()
       if name is None:
          self.__name=self.__class__.__name__
       else:
          self.__name=name

    def __str__(self):
       """
       Returns the name of this options object.

       :return: the name
       :rtype: ``str``
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
    print((bms.getHTML(filter=SimpleFilter())))

    bms.run(scale=4)
    print((bms.getHTML(filter=SimpleFilter())))

    bms.run(scale=[1,2])
    print((bms.getHTML(filter=SimpleFilter())))

