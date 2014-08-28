
##############################################################################
#
# Copyright (c) 2003-2014 by University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################
from __future__ import print_function
from __future__ import division


__copyright__="""Copyright (c) 2003-2014 by University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import os
import sys
import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.speckley import Rectangle, Brick, speckleycpp

@unittest.skipIf(getMPISizeWorld() != 1, "Speckley doesn't support MPI currently")
class Test_Speckley_Assemblers(unittest.TestCase):
    TOLERANCE = 1e-10

    def test_Brick_XY_single(self):
        for expanded in (True, False):
            for order in range(2,11):
                dom = Brick(order, 3, 3, 3, l0=6, l1=6, l2=6)
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
                        ("assembly for {0}expanded order %d failed with %g >= %g"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))

    def test_Rectangle_XY_single(self):
        for expanded in (True, False):
            for order in range(2,11):
                dom = Rectangle(order, 3, 3, l0=6, l1=6)
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
                        ("assembly for {0}expanded order %d failed with %g >= %g"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))

    def test_Brick_XY_system(self):
        for expanded in (True, False):
            for order in range(2,11):
                dom = Brick(order, 3, 3, 3, l0=6, l1=6, l2=6)
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
                        ("assembly for {0}expanded order %d failed with %g >= %g"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))

    def test_Rectangle_XY_system(self):
        for expanded in (True, False):
            for order in range(2,11):
                dom = Rectangle(order, 3, 3, l0=6, l1=6)
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
                        ("assembly for {0}expanded order %d failed with %g >= %g"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))

    def test_Brick_Du_Y_single(self):
        #solves for u in Du = Y, where D = 1, Y = 2
        for expanded in (True, False):
            for order in range(2,11):
                dom = Brick(order, 3, 3, 3)
                D = Data(1, ContinuousFunction(dom), expanded)
                Y = Data(2, ContinuousFunction(dom), expanded)

                Y = interpolate(Y, Function(dom))
                D = interpolate(D, Function(dom))
                
                rhs = Data(0, Solution(dom), expanded)
                lhs = Data(0, Solution(dom), expanded)

                dom.addToRHS(rhs, [("Y",Y)], dom.createAssembler("DefaultAssembler", []))
                dom.addToRHS(lhs, [("D",D)], dom.createAssembler("DefaultAssembler", []))
                
                res = Lsup((rhs/lhs)-2)
                self.assertLess(res, self.TOLERANCE,
                        ("assembly for {0}expanded order %d failed with %g >= %g"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))
    
    def test_Rectangle_Du_Y_single(self):
        #solves for u in Du = Y, where D = 1, Y = 2
        for expanded in (True, False):
            for order in range(2,11):
                dom = Rectangle(order, 3, 3)
                D = Data(1, ContinuousFunction(dom), expanded)
                Y = Data(2, ContinuousFunction(dom), expanded)

                Y = interpolate(Y, Function(dom))
                D = interpolate(D, Function(dom))
                
                rhs = Data(0, Solution(dom), expanded)
                lhs = Data(0, Solution(dom), expanded)

                dom.addToRHS(rhs, [("Y",Y)], dom.createAssembler("DefaultAssembler", []))
                dom.addToRHS(lhs, [("D",D)], dom.createAssembler("DefaultAssembler", []))
                
                res = Lsup((rhs/lhs)-2)
                self.assertLess(res, self.TOLERANCE,
                        ("assembly for {0}expanded order %d failed with %g >= %g"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))

    def test_Brick_Du_Y_system(self):
        #solves for u in Du = Y, where D = [1,2], Y = [2,4]
        for expanded in (True, False):
            for order in range(2,11):
                dom = Brick(order, 3, 3, 3)
                D = Data(1, (2,), ContinuousFunction(dom), expanded)
                D[1] = 2
                Y = Data(2, (2,), ContinuousFunction(dom), expanded)
                Y[1] = 4

                Y = interpolate(Y, Function(dom))
                D = interpolate(D, Function(dom))
                
                rhs = Data(0, (2,), Solution(dom), expanded)
                lhs = Data(0, (2,), Solution(dom), expanded)

                dom.addToRHS(rhs, [("Y",Y)], dom.createAssembler("DefaultAssembler", []))
                dom.addToRHS(lhs, [("D",D)], dom.createAssembler("DefaultAssembler", []))
                
                res = Lsup((rhs/lhs)-2)
                self.assertLess(res, self.TOLERANCE,
                        ("assembly for {0}expanded order %d failed with %g >= %g"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))

    def test_Rectangle_Du_Y_system(self):
        #solves for u in Du = Y, where D = [1,2], Y = [2,4]
        for expanded in (True, False):
            for order in range(2,11):
                dom = Rectangle(order, 3, 3)
                D = Data(1, (2,), ContinuousFunction(dom), expanded)
                D[1] = 2
                Y = Data(2, (2,), ContinuousFunction(dom), expanded)
                Y[1] = 4

                Y = interpolate(Y, Function(dom))
                D = interpolate(D, Function(dom))
                
                rhs = Data(0, (2,), Solution(dom), expanded)
                lhs = Data(0, (2,), Solution(dom), expanded)

                dom.addToRHS(rhs, [("Y",Y)], dom.createAssembler("DefaultAssembler", []))
                dom.addToRHS(lhs, [("D",D)], dom.createAssembler("DefaultAssembler", []))
                
                res = Lsup((rhs/lhs)-2)
                self.assertLess(res, self.TOLERANCE,
                        ("assembly for {0}expanded order %d failed with %g >= %g"%(order,
                        res, self.TOLERANCE)).format("" if expanded else "un-"))

@unittest.skipIf(getMPISizeWorld() != 1, "Speckley doesn't support MPI currently")
class Test_Speckley(unittest.TestCase):
    def test_Rectangle_Function_gradient(self): #expanded and non-expanded
        for expanded in [True, False]:
            for order in range(2,11):
                dom = Rectangle(order, 3, 3)
                x = Data(5, Function(dom), True)
                self.assertLess(Lsup(grad(x)), 1e-10,
                        "single component failure, order %d%s, %g >= 1e-10"%(order,
                        (" expanded" if expanded else ""), Lsup(grad(x))))
                for data in [[5,1], [-5,-1], [5,1,1e-5]]:
                    x = Data(data, Function(dom), True)
                    g = grad(x)
                    for n,d in enumerate(data):
                        self.assertLess(Lsup(g[n]), 1e-10,
                                "%d-component failure, order %d %sexpanded, %g >= 1e-10"%(len(data),
                                order, ("" if expanded else "un-"), Lsup(g[n])))

    def test_Brick_Function_gradient(self):
        for expanded in [True, False]:
            for order in range(2,11):
                dom = Brick(order, 3, 3, 3)
                x = Data(5, Function(dom), True)
                self.assertLess(Lsup(grad(x)), 1e-10,
                        "single component failure, order %d%s, %g >= 1e-10"%(order,
                        (" expanded" if expanded else ""), Lsup(grad(x))))
                for data in [[5,1], [-5,-1], [5,1,1e-5]]:
                    x = Data(data, Function(dom), True)
                    g = grad(x)
                    for n,d in enumerate(data):
                        self.assertLess(Lsup(g[n]), 1e-10,
                                "%d-component failure, order %d %sexpanded, %g >= 1e-10"%(len(data),
                                order, ("" if expanded else "un-"), Lsup(g[n])))


    def test_Rectangle_ContinuousFunction_gradient(self):
        for order in range(2,11):
            dom = Rectangle(order, 3, 3, l0=100, l1=100)
            X = dom.getX()
            u = X[0] + X[1] + 1
            v = Lsup(grad(u) - 1)
            self.assertLess(v, 1e-10, "order %d, %g >= 1e-10, %s"%(order, v, str(grad(u)-1)))
            for power in range(1, order+1):
                for power2 in range(1, order+1):
                    a = X[0]**power * X[1]**power2
                    da = grad(a)
                    first = Lsup(da[0] - power*X[0]**(power-1) * X[1]**power2)/Lsup(power*X[0]**(power-1) * X[1]**power2)
                    second = Lsup(da[1] - power2*X[1]**(power2-1) * X[0]**power)/Lsup(power2*X[1]**(power2-1) * X[0]**power)
                    self.assertLess(first, 1e-10, "order %d and degree %d,%d, %g >= 1e-9"%(order, power, power2, first))
                    self.assertLess(second, 1e-10, "order %d and degree %d,%d, %g >= 1e-9"%(order, power, power2, second))

    def test_Brick_ContinuousFunction_gradient(self):
        for order in range(2,11):
            dom = Brick(order, 3, 3, 3, l0=100, l1=100, l2=100)
            X = dom.getX()
            u = X[0] + X[1] + X[2] + 1
            v = Lsup(grad(u) - 1)
            self.assertLess(v, 1e-10, "order %d, %g >= 1e-10, %s"%(order, v, str(grad(u)-1)))
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
                            "order %d and degree %d,%d,%d, %g >= 1e-9"%(order,
                            power1, power2, power3, first))
                        self.assertLess(second, 1e-10,
                            "order %d and degree %d,%d,%d, %g >= 1e-9"%(order,
                            power1, power2, power3, second))
                        self.assertLess(third, 1e-10,
                            "order %d and degree %d,%d,%d, %g >= 1e-9"%(order,
                            power1, power2, power3, third))

    def test_Rectangle_interpolation_continuous_noncontinuous_and_back(self):
        for order in range(2,11):
            dom = Rectangle(order, 3, 3, l0=6, l1=6)
            original = Data(5, Function(dom), True)
            cont = interpolate(original, ContinuousFunction(dom))
            func = interpolate(cont, Function(dom))
            self.assertEqual(Lsup(original-func), 0,
                    "interpolation of constant, order %d: original and final not equal, %g != 0"%(order, Lsup(original-func)))
            x = dom.getX()
            original = x[0] + x[1] + 1
            cont = interpolate(original, ContinuousFunction(dom))
            func = interpolate(cont, Function(dom))
            self.assertEqual(Lsup(original-func), 0,
                    "interpolation of expanded, order %d: original and final not equal, %g != 0"%(order, Lsup(original-func)))
            original = whereZero(x[0]-2) + whereZero(x[1]-2)
            cont = interpolate(original, ContinuousFunction(dom))
            func = interpolate(cont, Function(dom))
            self.assertEqual(Lsup(original-func), 0,
                    "interpolation of point, order %d: original and final not equal, %g != 0"%(order, Lsup(original-func)))

    def test_Brick_interpolation_continuous_noncontinuous_and_back(self):
        for order in range(2,11):
            dom = Brick(order, 3, 3, 3, l0=6, l1=6, l2=6)
            original = Data(5, Function(dom), True)
            cont = interpolate(original, ContinuousFunction(dom))
            func = interpolate(cont, Function(dom))
            self.assertEqual(Lsup(original-func), 0,
                    "interpolation of constant, order %d: original and final not equal, %g != 0"%(order, Lsup(original-func)))
            x = dom.getX()
            original = x[0] + x[1] + x[2] + 1
            cont = interpolate(original, ContinuousFunction(dom))
            func = interpolate(cont, Function(dom))
            self.assertEqual(Lsup(original-func), 0,
                    "interpolation of expanded, order %d: original and final not equal, %g != 0"%(order, Lsup(original-func)))
            original = whereZero(x[0]-2) + whereZero(x[1]-2) + whereZero(x[2] - 2)
            cont = interpolate(original, ContinuousFunction(dom))
            func = interpolate(cont, Function(dom))
            self.assertEqual(Lsup(original-func), 0,
                    "interpolation of point, order %d: original and final not equal, %g != 0"%(order, Lsup(original-func)))

    def test_Rectangle_integration(self):
        for order in range(2,11):
            size = 6
            dom = Rectangle(order, 3, 3, l0=size, l1=size)
            X = dom.getX()
            for k in range(1, order * 2):
                for l in range(1, order * 2):
                    powered = X[0]**k * X[1]**l
                    integral = integrate(powered)
                    actual = size**(k+1) * size**(l+1) / ((k+1.)*(l+1.))
                    self.assertLess(abs(integral - actual)/actual, 1e-11, "too much variance in integral result (order %d, degrees %d %d)"%(order, k, l))

    def test_Brick_integration(self):
        for order in range(2,11):
            size = 6
            dom = Brick(order, 3, 3, 3, l0=6, l1=6, l2=6)
            X = dom.getX()
            for k in [1, order, order*2 - 1]:
                for l in [1, order, order*2 - 1]:
                    for m in [1, order, order*2 - 1]:
                        powered = X[0]**k * X[1]**l * X[2]**m
                        integral = integrate(powered)
                        actual = size**(k+1) * size**(l+1) * size**(m+1)/ ((k+1.)*(l+1.)*(m+1.))
                        res = abs(integral - actual)/actual
                        self.assertLess(res, 1e-11, "too much variance in integral result (order %d, degrees %d %d, %g >= 1e-11)"%(order, k, l, res))

    def test_Rectangle_MPI_construction(self):
        for order in range(2,11):
            dom = Rectangle(order, 2, 2*getMPISizeWorld(), d1=getMPISizeWorld())
            self.assertEqual(Lsup(dom.getX()[0] + dom.getX()[1]), 2, "invalid Lsup(getX()) for y-split order %d"%order)
            X = interpolate(dom.getX()[0], Function(dom)) + interpolate(dom.getX()[1], Function(dom))
            val = Lsup(interpolate(X, ContinuousFunction(dom)) - (dom.getX()[0] + dom.getX()[1]))
            self.assertLess(val, 1e-10, "interpolation failure for y-split order %d"%order)

            dom = Rectangle(order, 2*getMPISizeWorld(), 2, d0=getMPISizeWorld())
            self.assertEqual(Lsup(dom.getX()[0] + dom.getX()[1]), 2, "invalid Lsup(getX()) for x-split order %d"%order)
            X = interpolate(dom.getX()[0], Function(dom)) + interpolate(dom.getX()[1], Function(dom))
            val = Lsup(interpolate(X, ContinuousFunction(dom)) - (dom.getX()[0] + dom.getX()[1]))
            self.assertLess(val, 1e-10, "interpolation failure for x-split order %d"%order)

        #TODO need some multi-dim division tests with associated checks

    def test_Brick_MPI_construction(self):
        for order in range(2,11):
            dom = Brick(order, 2*getMPISizeWorld(), 2, 2, d0=getMPISizeWorld())
            self.assertEqual(Lsup(dom.getX()[0] + dom.getX()[1] + dom.getX()[2]),
                    3, "invalid Lsup(getX()) for x-split order %d"%order)
            X = interpolate(dom.getX()[0], Function(dom)) + interpolate(dom.getX()[1], Function(dom))
            val = Lsup(interpolate(X, ContinuousFunction(dom)) - (dom.getX()[0] + dom.getX()[1]))
            self.assertLess(val, 1e-10, "interpolation failure for x-split order %d"%order)

            dom = Brick(order, 2, 2*getMPISizeWorld(), 2, d1=getMPISizeWorld())
            self.assertEqual(Lsup(dom.getX()[0] + dom.getX()[1] + dom.getX()[2]),
                    3, "invalid Lsup(getX()) for y-split order %d"%order)
            X = interpolate(dom.getX()[0], Function(dom)) + interpolate(dom.getX()[1], Function(dom))
            val = Lsup(interpolate(X, ContinuousFunction(dom)) - (dom.getX()[0] + dom.getX()[1]))
            self.assertLess(val, 1e-10, "interpolation failure for y-split order %d"%order)
            
            dom = Brick(order, 2, 2, 2*getMPISizeWorld(), d2=getMPISizeWorld())
            self.assertEqual(Lsup(dom.getX()[0] + dom.getX()[1] + dom.getX()[2]),
                    3, "invalid Lsup(getX()) for z-split order %d"%order)
            X = interpolate(dom.getX()[0], Function(dom)) + interpolate(dom.getX()[1], Function(dom))
            val = Lsup(interpolate(X, ContinuousFunction(dom)) - (dom.getX()[0] + dom.getX()[1]))
            self.assertLess(val, 1e-10, "interpolation failure for z-split order %d"%order)
        
        #TODO need some multi-dim division tests with associated checks

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)

