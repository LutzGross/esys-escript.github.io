
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

import os
import sys
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.speckley import Rectangle, Brick, speckleycpp

class Test_Speckley_Assemblers(unittest.TestCase):
    TOLERANCE = 1e-10

    def test_Brick_XY_single(self):
        ranks = getMPISizeWorld()
        for expanded in (True, False):
            for order in range(2,11):
                dom = Brick(order, 3, 3*ranks, 3, l0=6, l1=6, l2=6, d1=ranks)
                Y = Data(3, Function(dom), expanded)
                X = Data(0, (3,), Function(dom), expanded)
                X[0] = dom.getX()[0]
                X[1] = dom.getX()[1]
                X[2] = dom.getX()[2]
                f = Data(0, Solution(dom), expanded)

                dom.addToRHS(f, [("Y",Y)],
                        dom.createAssembler("DefaultAssembler", []))
                dom.addToRHS(f, [("X", X)],
                        dom.createAssembler("DefaultAssembler", []))
                #nuke the boundary back to zero since it's not dealt with here
                for dim in range(3):
                    f *= whereNegative(dom.getX()[dim]-6)
                res = Lsup(f)
                self.assertLess(res, self.TOLERANCE,
                        ("assembly for {0}expanded order %d failed with %e >= %e"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))

    def test_Rectangle_XY_single(self):
        ranks = getMPISizeWorld()
        for expanded in (True, False):
            for order in range(2,11):
                dom = Rectangle(order, 3, 3*ranks, d1=ranks, l0=6, l1=6)
                Y = Data(2, Function(dom), expanded)
                X = Data(0, (2,), Function(dom), expanded)
                X[0] = dom.getX()[0]
                X[1] = dom.getX()[1]
                f = Data(0, Solution(dom), expanded)

                dom.addToRHS(f, [("Y",Y)],
                        dom.createAssembler("DefaultAssembler", []))
                dom.addToRHS(f, [("X", X)],
                        dom.createAssembler("DefaultAssembler", []))
                #nuke the boundary back to zero since it's not dealt with here
                for dim in range(2):
                    f *= whereNegative(dom.getX()[dim]-6)
                res = Lsup(f)
                self.assertLess(res, self.TOLERANCE,
                        ("assembly for {0}expanded order %d failed with %e >= %e"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))

    def test_Brick_XY_system(self):
        ranks = getMPISizeWorld()
        for expanded in (True, False):
            for order in range(2,11):
                dom = Brick(order, 3, 3*ranks, 3, l0=6, l1=6, l2=6, d1=ranks)
                Y = Data(1, (3,), Function(dom), expanded)
                X = Data(0, (3,3), Function(dom), expanded)
                X[0,0] = dom.getX()[0]
                X[1,1] = dom.getX()[1]
                X[2,2] = dom.getX()[2]

                f = Data(0, (3,), Solution(dom), expanded)

                dom.addToRHS(f, [("Y",Y)],
                        dom.createAssembler("DefaultAssembler", []))
                dom.addToRHS(f, [("X", X)],
                        dom.createAssembler("DefaultAssembler", []))
                #nuke the boundary back to zero since it's not dealt with here
                for dim in range(3):
                    f *= whereNegative(dom.getX()[dim]-6)
                res = Lsup(f)
                self.assertLess(res, self.TOLERANCE,
                        ("assembly for {0}expanded order %d failed with %e >= %e"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))

    def test_Rectangle_XY_system(self):
        ranks = getMPISizeWorld()
        for expanded in (True, False):
            for order in range(2,11):
                dom = Rectangle(order, 3, 3*ranks, l0=6, l1=6, d1=ranks)
                Y = Data(1, (2,), Function(dom), expanded)
                X = Data(0, (2,2), Function(dom), expanded)
                X[0,0] = dom.getX()[0]
                X[1,1] = dom.getX()[1]

                f = Data(0, (2,), Solution(dom), expanded)

                dom.addToRHS(f, [("Y",Y)],
                        dom.createAssembler("DefaultAssembler", []))
                dom.addToRHS(f, [("X", X)],
                        dom.createAssembler("DefaultAssembler", []))
                #nuke the boundary back to zero since it's not dealt with here
                f *= whereNegative(dom.getX()[0]-6)*whereNegative(dom.getX()[1]-6)
                res = Lsup(f)
                self.assertLess(res, self.TOLERANCE,
                        ("assembly for {0}expanded order %d failed with %e >= %e"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))

    def test_Brick_Du_Y_single(self):
        #solves for u in Du = Y, where D = 1, Y = 2
        ranks = getMPISizeWorld()
        for expanded in (True, False):
            for order in range(2,11):
                dom = Brick(order, 3, 3*ranks, 3, d1=ranks)
                D = Data(1, ContinuousFunction(dom), expanded)
                Y = Data(2, ContinuousFunction(dom), expanded)

                Y = interpolate(Y, Function(dom))
                D = interpolate(D, Function(dom))
                
                rhs = Data(0, Solution(dom), expanded)
                lhs = Data(0, Solution(dom), expanded)

                dom.addToRHS(rhs, [("Y",Y)],
                        dom.createAssembler("DefaultAssembler", []))
                dom.addToRHS(lhs, [("D",D)],
                        dom.createAssembler("DefaultAssembler", []))
                
                res = Lsup((rhs/lhs)-2)
                self.assertLess(res, self.TOLERANCE,
                        ("assembly for {0}expanded order %d failed with %e >= %e"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))
    
    def test_Rectangle_Du_Y_single(self):
        #solves for u in Du = Y, where D = 1, Y = 2
        ranks = getMPISizeWorld()
        for expanded in (True, False):
            for order in range(2,11):
                dom = Rectangle(order, 3, 3*ranks, d1=ranks)
                D = Data(1, ContinuousFunction(dom), expanded)
                Y = Data(2, ContinuousFunction(dom), expanded)

                Y = interpolate(Y, Function(dom))
                D = interpolate(D, Function(dom))
                
                rhs = Data(0, Solution(dom), expanded)
                lhs = Data(0, Solution(dom), expanded)

                dom.addToRHS(rhs, [("Y",Y)],
                        dom.createAssembler("DefaultAssembler", []))
                dom.addToRHS(lhs, [("D",D)],
                        dom.createAssembler("DefaultAssembler", []))
                
                res = Lsup((rhs/lhs)-2)
                self.assertLess(res, self.TOLERANCE,
                        ("assembly for {0}expanded order %d failed with %e >= %e"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))

    def test_Brick_Du_Y_system(self):
        #solves for u in Du = Y, where D = [1,2], Y = [2,4]
        ranks = getMPISizeWorld()
        for expanded in (True, False):
            for order in range(2,11):
                dom = Brick(order, 3, 3*ranks, 3, d1=ranks)
                D = Data(1, (2,), ContinuousFunction(dom), expanded)
                D[1] = 2
                Y = Data(2, (2,), ContinuousFunction(dom), expanded)
                Y[1] = 4

                Y = interpolate(Y, Function(dom))
                D = interpolate(D, Function(dom))
                
                rhs = Data(0, (2,), Solution(dom), expanded)
                lhs = Data(0, (2,), Solution(dom), expanded)

                dom.addToRHS(rhs, [("Y",Y)],
                        dom.createAssembler("DefaultAssembler", []))
                dom.addToRHS(lhs, [("D",D)],
                        dom.createAssembler("DefaultAssembler", []))
                
                res = Lsup((rhs/lhs)-2)
                self.assertLess(res, self.TOLERANCE,
                        ("assembly for {0}expanded order %d failed with %e >= %e"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))

    def test_Rectangle_Du_Y_system(self):
        #solves for u in Du = Y, where D = [1,2], Y = [2,4]
        ranks = getMPISizeWorld()
        for expanded in (True, False):
            for order in range(2,11):
                dom = Rectangle(order, 3, 3*ranks, d1=ranks)
                D = Data(1, (2,), ContinuousFunction(dom), expanded)
                D[1] = 2
                Y = Data(2, (2,), ContinuousFunction(dom), expanded)
                Y[1] = 4

                Y = interpolate(Y, Function(dom))
                D = interpolate(D, Function(dom))
                
                rhs = Data(0, (2,), Solution(dom), expanded)
                lhs = Data(0, (2,), Solution(dom), expanded)

                dom.addToRHS(rhs, [("Y",Y)],
                        dom.createAssembler("DefaultAssembler", []))
                dom.addToRHS(lhs, [("D",D)],
                        dom.createAssembler("DefaultAssembler", []))
                
                res = Lsup((rhs/lhs)-2)
                self.assertLess(res, self.TOLERANCE,
                        ("assembly for {0}expanded order %d failed with %e >= %e"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))

class Test_Speckley(unittest.TestCase):
    TOLERANCE = 5e-10
    def test_Rectangle_ReducedFunction(self):
        ranks = getMPISizeWorld()
        for order in range(2, 11):
            dom = Rectangle(order, 3, 3*ranks, l0=3, l1=3*ranks, d1=ranks)
            X = dom.getX()
            redData = interpolate(X, ReducedFunction(dom))
            data = [(interpolate(redData, ReducedFunction(dom)), "ReducedFunction"),
                    (interpolate(redData, Function(dom)), "Function"),
                    (interpolate(redData, ContinuousFunction(dom)), "ContinuousFunction")]
            for d, fs in data:
                self.assertLess(inf(d-[0.5]*2), self.TOLERANCE,
                        "reduced->%s failure with order %d: %e != 0"%(fs, order, inf(d-[0.5]*2)))
                self.assertLess(sup(d[0]+0.5) - 3, self.TOLERANCE,
                        "reduced->%s failure with order %d: %e != 3"%(fs, order, sup(d[0]+0.5)))
                self.assertLess(sup(d[1]+0.5) - 3*ranks, self.TOLERANCE,
                        "reduced->%s failure with order %d: %e >= %e"%(fs, order, sup(d[1]+0.5)-3*ranks, self.TOLERANCE))

    def test_Brick_ReducedFunction(self):
        ranks = getMPISizeWorld()
        for order in range(2, 11):
            dom = Brick(order, 3, 3*ranks, 3, l0=3, l1=3*ranks, l2=3, d1=ranks)
            X = dom.getX()
            redData = interpolate(X, ReducedFunction(dom))
            data = [(interpolate(redData, ReducedFunction(dom)), "ReducedFunction"),
                    (interpolate(redData, Function(dom)), "Function"),
                    (interpolate(redData, ContinuousFunction(dom)), "ContinuousFunction")]
            for d, fs in data:
                self.assertLess(inf(d-[0.5]*3), self.TOLERANCE,
                        "reduced->%s failure with order %d: %e != 0"%(fs, order, inf(d-[0.5]*3)))
                self.assertLess(sup(d[0]+0.5) - 3, self.TOLERANCE,
                        "reduced->%s failure with order %d: %e != 3"%(fs, order, sup(d[0]+0.5)))
                self.assertLess(sup(d[1]+0.5) - 3*ranks, self.TOLERANCE,
                        "reduced->%s failure with order %d: %e >= %e"%(fs, order, sup(d[1]+0.5)-3*ranks, self.TOLERANCE))
                self.assertLess(sup(d[2]+0.5) - 3, self.TOLERANCE,
                        "reduced->%s failure with order %d: %e != 3"%(fs, order, sup(d[2]+0.5)))
                
    def test_Rectangle_Function_gradient_AntiSimetric(self): #expanded and non-expanded
        ranks = getMPISizeWorld()
        for expanded in [True, False]:
            for order in range(2,11):
                dom = Rectangle(order, 3, 3*ranks, d1=ranks)
                x=ContinuousFunction(dom).getX()
                m=-2*1j+4
                b=3+8*1j
                f=Vector(0., Solution(dom))
                sample_data=[[0.000000000+0.000000000j, -4.000000000+2.000000000j],[4.000000000-2.000000000j, 0.000000000+0.000000000j]]
                sample=Data(sample_data, Function(dom), True)
                sample[0,0]=0.000000000+0.000000000j
                sample[0,1]=-4.000000000+2.000000000j
                sample[1,0]=4.000000000-2.000000000j
                sample[1,1]=0.000000000+0.000000000j
                for i in range(dom.getDim()):
                    for j in range(dom.getDim()):
                        f[i]+=m*(i-j)*x[j]+b-i*j
                g=grad(f)
                self.assertLess(Lsup(g-sample), 1e-9,"single component failure, order %d%s, %e >= 1e-10"%(order,(" expanded" if expanded else ""), Lsup(g-sample)))
                
    def test_Rectangle_Function_gradient_Simetric(self): #expanded and non-expanded
        ranks = getMPISizeWorld()
        for expanded in [True, False]:
            for order in range(2,11):
                dom = Rectangle(order, 3, 3*ranks, d1=ranks)
                x=ContinuousFunction(dom).getX()
                m=-2*1j+4
                b=3+8*1j
                f=Vector(0., Solution(dom))
                sample_data=[[0.000000000+0.000000000j, 4.000000000-2.000000000j],[4.000000000-2.000000000j, 8.000000000-4.000000000j]]
                sample=Data(sample_data, Function(dom), True)
                sample[0,0]=0.000000000+0.000000000j
                sample[0,1]=4.000000000-2.000000000j
                sample[1,0]=4.000000000-2.000000000j
                sample[1,1]=8.000000000-4.000000000j
                for i in range(dom.getDim()):
                    for j in range(dom.getDim()):
                        f[i]+=m*(i+j)*x[j]+b-i*j
                g=grad(f)
                self.assertLess(Lsup(g-sample), 1e-9,"single component failure, order %d%s, %e >= 1e-10"%(order,(" expanded" if expanded else ""), Lsup(g-sample)))         

    def test_Brick_Function_gradient_AntiSimetric(self): #expanded and non-expanded
        ranks = getMPISizeWorld()
        for expanded in [True, False]:
            for order in range(2,11):
                dom = Brick(order, 3, 3*ranks, 3, d1=ranks)
                x=ContinuousFunction(dom).getX()
                m=-2*1j+4
                b=3+8*1j
                f=Vector(0., Solution(dom))
                sample_data=[[0.000000000+0.000000000j, -4.000000000+2.000000000j, -8.000000000+4.000000000j],[4.000000000-2.000000000j, 0.000000000+0.000000000j, -4.000000000+2.000000000j],[8.000000000-4.000000000j, 4.000000000-2.000000000j, 0.000000000+0.000000000j]]
                sample=Data(sample_data, Function(dom), True)
                sample[0,0]=0.000000000+0.000000000j
                sample[0,1]=-4.000000000+2.000000000j
                sample[0,2]=-8.000000000+4.000000000j
                sample[1,0]=4.000000000-2.000000000j
                sample[1,1]=0.000000000+0.000000000j
                sample[1,2]=-4.000000000+2.000000000j
                sample[2,0]=8.000000000-4.000000000j
                sample[2,1]=4.000000000-2.000000000j
                sample[2,2]=0.000000000+0.000000000j
                for i in range(dom.getDim()):
                    for j in range(dom.getDim()):
                        f[i]+=m*(i-j)*x[j]+b-i*j
                g=grad(f)        
                self.assertLess(Lsup(g-sample), 1e-9,"single component failure, order %d%s, %e >= 1e-10"%(order,(" expanded" if expanded else ""), Lsup(g-sample)))  
                
    def test_Brick_Function_gradient_Simetric(self): #expanded and non-expanded
        ranks = getMPISizeWorld()
        for expanded in [True, False]:
            for order in range(2,11):
                dom = Brick(order, 3, 3*ranks, 3, d1=ranks)
                x=ContinuousFunction(dom).getX()
                m=-2*1j+4
                b=3+8*1j
                f=Vector(0., Solution(dom))
                sample_data=[[0.000000000+0.000000000j, 4.000000000-2.000000000j, 8.000000000-4.000000000j],[4.000000000-2.000000000j, 8.000000000-4.000000000j, 12.000000000-6.000000000j],[8.000000000-4.000000000j, 12.000000000-6.000000000j, 16.000000000-8.000000000j]]
                sample=Data(sample_data, Function(dom), True)
                sample[0,0]=0.000000000+0.000000000j
                sample[0,1]=4.000000000-2.000000000j
                sample[0,2]=8.000000000-4.000000000j
                sample[1,0]=4.000000000-2.000000000j
                sample[1,1]=8.000000000-4.000000000j
                sample[1,2]=12.000000000-6.000000000j
                sample[2,0]=8.000000000-4.000000000j
                sample[2,1]=12.000000000-6.000000000j
                sample[2,2]=16.000000000-8.000000000j
                for i in range(dom.getDim()):
                    for j in range(dom.getDim()):
                        f[i]+=m*(i+j)*x[j]+b-i*j
                g=grad(f)        
                self.assertLess(Lsup(g-sample), 1e-9,"single component failure, order %d%s, %e >= 1e-10"%(order,(" expanded" if expanded else ""), Lsup(g-sample)))                 
           

    def test_Rectangle_Function_gradient(self): #expanded and non-expanded
        ranks = getMPISizeWorld()
        for expanded in [True, False]:
            for order in range(2,11):
                dom = Rectangle(order, 3, 3*ranks, d1=ranks)
                x = Data(5, Function(dom), True)
                self.assertLess(Lsup(grad(x)), 1e-10,
                        "single component failure, order %d%s, %e >= 1e-10"%(order,
                        (" expanded" if expanded else ""), Lsup(grad(x))))
                for data in [[5,1], [-5,-1], [5,1,1e-5]]:
                    x = Data(data, Function(dom), True)
                    g = grad(x)
                    for n,d in enumerate(data):
                        self.assertLess(Lsup(g[n]), 1e-10,
                                "%d-component failure, order %d %sexpanded, %e >= 1e-10"%(len(data),
                                order, ("" if expanded else "un-"), Lsup(g[n])))

    def test_Brick_Function_gradient(self):
        ranks = getMPISizeWorld()
        for expanded in [True, False]:
            for order in range(2,11):
                dom = Brick(order, 3, 3*ranks, 3, d1=ranks)
                x = Data(5, Function(dom), True)
                self.assertLess(Lsup(grad(x)), 1e-10,
                        "single component failure, order %d%s, %e >= 1e-10"%(order,
                        (" expanded" if expanded else ""), Lsup(grad(x))))
                for data in [[5,1], [-5,-1], [5,1,1e-5]]:
                    x = Data(data, Function(dom), True)
                    g = grad(x)
                    for n,d in enumerate(data):
                        self.assertLess(Lsup(g[n]), 1e-10,
                                "%d-component failure, order %d %sexpanded, %e >= 1e-10"%(len(data),
                                order, ("" if expanded else "un-"), Lsup(g[n])))


    def test_Rectangle_ContinuousFunction_gradient(self):
        ranks = getMPISizeWorld()
        for order in range(2,11):
            dom = Rectangle(order, 3, 3*ranks, d1=ranks, l0=100, l1=100)
            X = dom.getX()
            u = X[0] + X[1] + 1
            v = Lsup(grad(u) - 1)
            self.assertLess(v, 1e-10, "order %d, %e >= 1e-10, %s"%(order, v, str(grad(u)-1)))
            for power in range(1, order+1):
                for power2 in range(1, order+1):
                    a = X[0]**power * X[1]**power2
                    da = grad(a)
                    first = Lsup(da[0] - power*X[0]**(power-1) * X[1]**power2) \
                            /Lsup(power*X[0]**(power-1) * X[1]**power2)
                    second = Lsup(da[1] - power2*X[1]**(power2-1) * X[0]**power) \
                            /Lsup(power2*X[1]**(power2-1) * X[0]**power)
                    self.assertLess(first, 1e-10,
                            "order %d and degree %d,%d, %e >= 1e-9"%(order,
                            power, power2, first))
                    self.assertLess(second, 1e-10,
                            "order %d and degree %d,%d, %e >= 1e-9"%(order,
                            power, power2, second))

    def test_Brick_ContinuousFunction_gradient(self):
        ranks = getMPISizeWorld()
        for order in range(2,11):
            dom = Brick(order, 3, 3*ranks, 3, l0=100, l1=100, l2=100, d1=ranks)
            X = dom.getX()
            u = X[0] + X[1] + X[2] + 1
            v = Lsup(grad(u) - 1)
            self.assertLess(v, 1e-10, "order %d, %e >= 1e-10, %s"%(order, v,
                    str(grad(u)-1)))
            for power1 in range(1, order+1, order//2):
                for power2 in range(1, order+1, order//2):
                    for power3 in range(1, order+1, order//2):
                        a = X[0]**power1 * X[1]**power2 * X[2]**power3
                        da = grad(a)

                        temp = power1*X[0]**(power1-1) * X[1]**power2*X[2]**power3
                        first = Lsup(da[0] - temp) / Lsup(temp)
                        temp = power2*X[1]**(power2-1) * X[0]**power1*X[2]**power3
                        second = Lsup(da[1] - temp) / Lsup(temp)
                        temp = power3*X[2]**(power3-1) * X[0]**power1*X[1]**power2
                        third = Lsup(da[2] - temp) / Lsup(temp)
                                
                        self.assertLess(first, 1e-10,
                            "order %d and degree %d,%d,%d, %e >= 1e-9"%(order,
                            power1, power2, power3, first))
                        self.assertLess(second, 1e-10,
                            "order %d and degree %d,%d,%d, %e >= 1e-9"%(order,
                            power1, power2, power3, second))
                        self.assertLess(third, 1e-10,
                            "order %d and degree %d,%d,%d, %e >= 1e-9"%(order,
                            power1, power2, power3, third))

    def test_Rectangle_interpolation_continuous_noncontinuous_and_back(self):
        ranks = getMPISizeWorld()
        for order in range(2,11):
            dom = Rectangle(order, 3, 3*ranks, d1=ranks, l0=6, l1=6)
            original = Data(5, Function(dom), True)
            cont = interpolate(original, ContinuousFunction(dom))
            func = interpolate(cont, Function(dom))
            self.assertEqual(Lsup(original-func), 0,
                    "interpolation of constant, order %d: original and final not equal, %e != 0"%(order, Lsup(original-func)))
            x = dom.getX()
            original = x[0] + x[1] + 1
            cont = interpolate(original, ContinuousFunction(dom))
            func = interpolate(cont, Function(dom))
            self.assertEqual(Lsup(original-func), 0,
                    "interpolation of expanded, order %d: original and final not equal, %e != 0"%(order, Lsup(original-func)))
            original = whereZero(x[0]-2) + whereZero(x[1]-2)
            cont = interpolate(original, ContinuousFunction(dom))
            func = interpolate(cont, Function(dom))
            self.assertEqual(Lsup(original-func), 0,
                    "interpolation of point, order %d: original and final not equal, %e != 0"%(order, Lsup(original-func)))

    def xtest_Brick_interpolation_continuous_noncontinuous_and_back(self):
        ranks = getMPISizeWorld()
        for order in range(2,11):
            dom = Brick(order, 3, 3*ranks, 3, l0=6, l1=ranks, l2=6, d1=ranks)
            original = Data(5, Function(dom), True)
            cont = interpolate(original, ContinuousFunction(dom))
            func = interpolate(cont, Function(dom))
            self.assertEqual(Lsup(original-func), 0,
                    "interpolation of constant, order %d: original and final not equal, %e != 0"%(order, Lsup(original-func)))
            x = dom.getX()
            original = x[0] + x[1] + x[2] + 1
            cont = interpolate(original, ContinuousFunction(dom))
            func = interpolate(cont, Function(dom))
            self.assertEqual(Lsup(original-func), 0,
                    "interpolation of expanded, order %d: original and final not equal, %e != 0"%(order, Lsup(original-func)))
            original = whereZero(x[0]-2) + whereZero(x[1]-2) + whereZero(x[2] - 2)
            cont = interpolate(original, ContinuousFunction(dom))
            func = interpolate(cont, Function(dom))
            self.assertEqual(Lsup(original-func), 0,
                    "interpolation of point, order %d: original and final not equal, %e != 0"%(order, Lsup(original-func)))

    def xtest_Rectangle_integration(self):
        ranks = getMPISizeWorld()
        for order in range(2,11):
            size = 6
            dom = Rectangle(order, 3, 3*ranks, d1=ranks, l0=size, l1=size)
            X = dom.getX()
            for k in range(1, order * 2):
                for l in range(1, order * 2):
                    powered = X[0]**k * X[1]**l
                    integral = integrate(powered)
                    actual = size**(k+1) * size**(l+1) / ((k+1.)*(l+1.))
                    self.assertLess(abs(integral - actual)/actual, 1e-11,
                            "too much variance in integral result (order %d, degrees %d %d)"%(order, k, l))

    def xtest_Brick_integration(self):
        ranks = getMPISizeWorld()
        for order in range(2,11):
            size = 6
            dom = Brick(order, 3, 3*ranks, 3, l0=6, l1=6, l2=6, d1=ranks)
            X = dom.getX()
            for k in [1, order, order*2 - 1]:
                for l in [1, order, order*2 - 1]:
                    for m in [1, order, order*2 - 1]:
                        powered = X[0]**k * X[1]**l * X[2]**m
                        integral = integrate(powered)
                        actual = size**(k+1) * size**(l+1) * size**(m+1) \
                                /((k+1.)*(l+1.)*(m+1.))
                        res = abs(integral - actual)/actual
                        self.assertLess(res, 1e-11,
                                "too much variance in integral result (order %d, degrees %d %d, %e >= 1e-11)"%(order,
                                k, l, res))

    @unittest.skipIf(getMPISizeWorld() == 1,
        "only works with more than one rank")
    def xtest_Rectangle_MPI_construction_single_dimensional(self):
        ranks = getMPISizeWorld()
        for order in range(2,11):
            dom = Rectangle(order, 2, 2*ranks, l1=ranks, d1=ranks)
            self.assertEqual(Lsup(dom.getX()[1]+dom.getX()[0]), ranks+1,
                    "invalid getX() for y-splits order %d"%order)

            filt = whereZero(dom.getX()[1] - 1)
            for i in range(2,ranks):
                filt += whereZero(dom.getX()[1] - i)
            if getMPIRankWorld() % 2:
                filt *= -1
            d = Vector(0, Function(dom))
            d[0] = 1 * filt
            d[1] = 10 * filt
            X = interpolate(d, Function(dom))
            res = interpolate(X, ContinuousFunction(dom))
            val = Lsup(res)
            self.assertEqual(val, 0,
                    "summation stage failure for y-splits in order %d"%order)

            dom = Rectangle(2, 2*ranks, 2, l0=ranks, d0=ranks)
            self.assertEqual(Lsup(dom.getX()[1]+dom.getX()[0]), ranks+1,
                    "invalid getX() for x-splits order %d"%order)
            filt = whereZero(dom.getX()[0] - 1)
            for i in range(2,getMPISizeWorld()):
                filt += whereZero(dom.getX()[0] - i)
            if getMPIRankWorld() % 2:
                filt *= -1
            d = Vector(0, Function(dom))
            d[0] = 1 * filt
            d[1] = 10 * filt
            X = interpolate(d, Function(dom))
            res = interpolate(X, ContinuousFunction(dom))
            val = Lsup(res)
            self.assertEqual(val, 0,
                    "summation stage failure for x-splits in order %d"%order)

            X = interpolate(d+2, Function(dom))
            res = interpolate(X, ContinuousFunction(dom))
            val = Lsup(res-2)
            self.assertEqual(val, 0,
                    "averaging stage failure for x-splits in order %d"%order)

    @unittest.skipIf(getMPISizeWorld() != 4, "requires 4 ranks")
    def xtest_Rectangle_MPI_construction_multi_dimensional(self):
        ranks = getMPISizeWorld()
        half = 2 #change if ranks != 4 (sqrt(ranks))
        for order in range(2, 11):
            dom = Rectangle(order, ranks, ranks, l0=half, l1=half,
                            d0=half, d1=half)
            self.assertEqual(Lsup(dom.getX()[0]), half,
                    "invalid getX() for multidimensional splits in order %d"%order)
            self.assertEqual(Lsup(dom.getX()[1]), half,
                    "invalid getX() for multidimensional splits in order %d"%order)
            xfilt = whereZero(dom.getX()[0] - 1) + whereZero(dom.getX()[1] - 1)
            for i in range(2,half):
                xfilt += whereZero(dom.getX()[0] - i)
                xfilt += whereZero(dom.getX()[1] - i)
            xfilt = whereNonZero(xfilt)
            if getMPIRankWorld() in [1,2]: #change if ranks != 4
                xfilt *= -1
            X = interpolate(xfilt, Function(dom))
            res = interpolate(X, ContinuousFunction(dom))
            val = Lsup(res)
            self.assertEqual(val, 0,
                    "summation failure for mixed-splits in order %d"%order)
            X = interpolate(xfilt+2, Function(dom))
            res = interpolate(X, ContinuousFunction(dom))
            val = Lsup(res-2)
            self.assertEqual(val, 0,
                    "averaging failure for mixed-splits in order %d"%order)

    @unittest.skipIf(getMPISizeWorld() == 1,
        "only works with more than one rank")
    def xtest_Brick_MPI_construction(self):
        for order in range(2,11):
            dom = Brick(order, 2*getMPISizeWorld(), 2, 2, d0=getMPISizeWorld())
            self.assertEqual(Lsup(dom.getX()[0] + dom.getX()[1] + dom.getX()[2]),
                    3, "invalid Lsup(getX()) for x-split order %d"%order)
            X = interpolate(dom.getX()[0], Function(dom)) \
                    + interpolate(dom.getX()[1], Function(dom))
            val = Lsup(interpolate(X, ContinuousFunction(dom)) \
                    - (dom.getX()[0] + dom.getX()[1]))
            self.assertLess(val, 1e-10,
                    "interpolation failure for x-split order %d"%order)

            dom = Brick(order, 2, 2*getMPISizeWorld(), 2, d1=getMPISizeWorld())
            self.assertEqual(Lsup(dom.getX()[0] + dom.getX()[1] + dom.getX()[2]),
                    3, "invalid Lsup(getX()) for y-split order %d"%order)
            X = interpolate(dom.getX()[0], Function(dom)) \
                    + interpolate(dom.getX()[1], Function(dom))
            val = Lsup(interpolate(X, ContinuousFunction(dom)) \
                    - (dom.getX()[0] + dom.getX()[1]))
            self.assertLess(val, 1e-10,
                    "interpolation failure for y-split order %d"%order)
            
            dom = Brick(order, 2, 2, 2*getMPISizeWorld(), d2=getMPISizeWorld())
            self.assertEqual(Lsup(dom.getX()[0] + dom.getX()[1] + dom.getX()[2]),
                    3, "invalid Lsup(getX()) for z-split order %d"%order)
            X = interpolate(dom.getX()[0], Function(dom)) \
                    + interpolate(dom.getX()[1], Function(dom))
            val = Lsup(interpolate(X, ContinuousFunction(dom)) \
                    - (dom.getX()[0] + dom.getX()[1]))
            self.assertLess(val, 1e-10,
                    "interpolation failure for z-split order %d"%order)
    
    @unittest.skipIf(getMPISizeWorld() == 1, "requires multiple MPI processes")
    def xtest_Brick_singledim_subdivision(self):
        ranks = getMPISizeWorld()
        for dim in range(0,3):
            label = ["x","y","z"][dim]
            size = [2,2,2]
            size[dim] *= ranks
            lengths = [1,1,1]
            lengths[dim] *= ranks
            splits = [1,1,1]
            splits[dim] *= ranks
            
            for order in range(2, 11):
                dom = Brick(order, size[0], size[1], size[2], 
                                   l0=lengths[0], l1=lengths[1], l2=lengths[2],
                                   d0=splits[0], d1=splits[1], d2=splits[2])
                self.assertEqual(Lsup(dom.getX()[1]+dom.getX()[0]+dom.getX()[2]),
                        ranks+2, "invalid getX() for %s-splits order %d"%(\
                        label, order))

                filt = whereZero(dom.getX()[dim] - 1)
                for i in range(2,ranks):
                    filt += whereZero(dom.getX()[0] - i)
                if getMPIRankWorld() % 2:
                    filt *= -1
                d = Vector(0, Function(dom))
                d[0] = 1 * filt
                d[1] = 10 * filt
                d[2] = 100 * filt
                X = interpolate(d, Function(dom))
                res = interpolate(X, ContinuousFunction(dom))
                val = Lsup(res)
                self.assertEqual(val, 0, 
                        "summation stage failure for %s-splits in order %d,"%(\
                        label, order))
                X = interpolate(d+2, Function(dom))
                res = interpolate(X, ContinuousFunction(dom))
                val = Lsup(res-2)
                self.assertEqual(val, 0, 
                        "averaging stage failure for %s-splits in order %d,"%(\
                        label, order))    

    @unittest.skipIf(getMPISizeWorld() != 4, "requires 4 ranks exactly")
    def xtest_Brick_multidim_subdivision(self):
        ranks = getMPISizeWorld()
        half = 2 #change if ranks != 4 (sqrt(ranks))
        for order in range(2, 11):
            for dom,dim1,dim2 in [
                (Brick(order, 2*half, 2*half, 2, l0=half, l1=half,
                                d0=half, d1=half), 0,1),
                (Brick(order, 2*half, 2, 2*half, l0=half, l2=half,
                                d0=half, d2=half), 0,2),
                (Brick(order, 2, 2*half, 2*half, l1=half, l2=half,
                                d1=half, d2=half), 1,2)]:
                self.assertEqual(ranks + 1,
                        Lsup(dom.getX()[0] + dom.getX()[1] + dom.getX()[2]),
                        "invalid getX() for multidimensional split " + \
                        "(dims %d,%d) in order %d"%(dim1,dim2,order))
                xfilt = whereZero(dom.getX()[dim1] - 1) \
                        + whereZero(dom.getX()[dim2] - 1)
                for i in range(2,half):
                    xfilt += whereZero(dom.getX()[dim1] - i)
                    xfilt += whereZero(dom.getX()[dim2] - i)
                xfilt = whereNonZero(xfilt)
                if getMPIRankWorld() in [1,2]: #change if ranks != 4
                    xfilt *= -1
                X = interpolate(xfilt, Function(dom))
                res = interpolate(X, ContinuousFunction(dom))
                val = Lsup(res)
                self.assertEqual(val, 0, 
                    "summation failure for mixed-splits " \
                    + "(dims %d,%d) in order %d"%(dim1,dim2,order))
                X = interpolate(xfilt+2, Function(dom))
                res = interpolate(X, ContinuousFunction(dom))
                val = Lsup(res-2)
                self.assertEqual(val, 0, 
                    "averaging failure for mixed-splits "\
                    + "(dims %d,%d) in order %d"%(dim1,dim2,order))

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

