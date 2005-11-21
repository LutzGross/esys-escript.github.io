esysroot = Dir('#.')
Export(["esysroot"])

libinstall = Dir('#lib')
Export(["libinstall"])

SConscript(['tools/CppUnitTest/SConstruct','tools/mmio/SConstruct','esysUtils/SConstruct','escript/SConstruct','paso/SConstruct','finley/SConstruct','bruce/SConstruct'], duplicate=0)
