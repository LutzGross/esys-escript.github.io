
########################################################
#
# Copyright (c) 2003-2009 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2009 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

# ******NOTICE***************
# optimize.py module by Travis E. Oliphant
#
# You may copy and use this module as you see fit with no
# guarantee implied provided you keep this notice in all copies.
# *****END NOTICE************

# A collection of optimization algorithms.  Version 0.3.1

# Minimization routines
"""optimize.py

A collection of general-purpose optimization routines using numarrayeric

fmin        ---      Nelder-Mead Simplex algorithm (uses only function calls)
fminBFGS    ---      Quasi-Newton method (uses function and gradient)
fminNCG     ---      Line-search Newton Conjugate Gradient (uses function, gradient and hessian (if it's provided))

"""
import numarray
__version__="0.3.1"

def para(x):
    return (x[0]-2.)**2

def para_der(x):
    return 2*(x-2.)

def para_hess(x):
    return numarray.ones((1,1))*2.

def para_hess_p(x,p):
    return 2*p

def rosen(x):  # The Rosenbrock function
    return numarray.sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)

def rosen_der(x):
    xm = x[1:-1]
    xm_m1 = x[:-2]
    xm_p1 = x[2:]
    der = numarray.zeros(x.shape,x.typecode())
    der[1:-1] = 200*(xm-xm_m1**2) - 400*(xm_p1 - xm**2)*xm - 2*(1-xm)
    der[0] = -400*x[0]*(x[1]-x[0]**2) - 2*(1-x[0])
    der[-1] = 200*(x[-1]-x[-2]**2)
    return der

def rosen3_hess_p(x,p):
    assert(len(x)==3)
    assert(len(p)==3)
    hessp = numarray.zeros((3,),x.typecode())
    hessp[0] = (2 + 800*x[0]**2 - 400*(-x[0]**2 + x[1])) * p[0] \
               - 400*x[0]*p[1] \
               + 0
    hessp[1] = - 400*x[0]*p[0] \
               + (202 + 800*x[1]**2 - 400*(-x[1]**2 + x[2]))*p[1] \
               - 400*x[1] * p[2]
    hessp[2] = 0 \
               - 400*x[1] * p[1] \
               + 200 * p[2]
    
    return hessp

def rosen3_hess(x):
    assert(len(x)==3)
    hessp = numarray.zeros((3,3),x.typecode())
    hessp[0,:] = [2 + 800*x[0]**2 -400*(-x[0]**2 + x[1]), -400*x[0], 0]
    hessp[1,:] = [-400*x[0], 202+800*x[1]**2 -400*(-x[1]**2 + x[2]), -400*x[1]]
    hessp[2,:] = [0,-400*x[1], 200]
    return hessp
    
        
def fmin(func, x0, args=(), xtol=1e-4, ftol=1e-4, maxiter=None, maxfun=None, fulloutput=0, printmessg=1):
    """xopt,{fval,warnflag} = fmin(function, x0, args=(), xtol=1e-4, ftol=1e-4,
    maxiter=200*len(x0), maxfun=200*len(x0), fulloutput=0, printmessg=0)

    Uses a Nelder-Mead Simplex algorithm to find the minimum of function
    of one or more variables.
    """
    x0 = numarray.asarray(x0)
    assert (len(x0.shape)==1)
    N = len(x0)
    if maxiter is None:
        maxiter = N * 200
    if maxfun is None:
        maxfun = N * 200

    rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
    one2np1 = range(1,N+1)

    sim = numarray.zeros((N+1,N),x0.typecode())
    fsim = numarray.zeros((N+1,),'d')
    sim[0] = x0
    fsim[0] = apply(func,(x0,)+args)
    nonzdelt = 0.05
    zdelt = 0.00025
    for k in range(0,N):
        y = numarray.array(x0,copy=1)
        if y[k] != 0:
            y[k] = (1+nonzdelt)*y[k]
        else:
            y[k] = zdelt

        sim[k+1] = y
        f = apply(func,(y,)+args)
        fsim[k+1] = f

    ind = numarray.argsort(fsim)
    fsim = numarray.take(fsim,ind)     # sort so sim[0,:] has the lowest function value
    sim = numarray.take(sim,ind,0)
    
    iterations = 1
    funcalls = N+1
    
    while (funcalls < maxfun and iterations < maxiter):
        if (max(numarray.ravel(numarray.absolute(sim[1:]-sim[0]))) <= xtol \
            and max(numarray.absolute(fsim[0]-fsim[1:])) <= ftol):
            break

        xbar = numarray.add.reduce(sim[:-1],0) / N
        xr = (1+rho)*xbar - rho*sim[-1]
        fxr = apply(func,(xr,)+args)
        funcalls = funcalls + 1
        doshrink = 0

        if fxr < fsim[0]:
            xe = (1+rho*chi)*xbar - rho*chi*sim[-1]
            fxe = apply(func,(xe,)+args)
            funcalls = funcalls + 1

            if fxe < fxr:
                sim[-1] = xe
                fsim[-1] = fxe
            else:
                sim[-1] = xr
                fsim[-1] = fxr
        else: # fsim[0] <= fxr
            if fxr < fsim[-2]:
                sim[-1] = xr
                fsim[-1] = fxr
            else: # fxr >= fsim[-2]
                # Perform contraction
                if fxr < fsim[-1]:
                    xc = (1+psi*rho)*xbar - psi*rho*sim[-1]
                    fxc = apply(func,(xc,)+args)
                    funcalls = funcalls + 1

                    if fxc <= fxr:
                        sim[-1] = xc
                        fsim[-1] = fxc
                    else:
                        doshrink=1
                else:
                    # Perform an inside contraction
                    xcc = (1-psi)*xbar + psi*sim[-1]
                    fxcc = apply(func,(xcc,)+args)
                    funcalls = funcalls + 1

                    if fxcc < fsim[-1]:
                        sim[-1] = xcc
                        fsim[-1] = fxcc
                    else:
                        doshrink = 1

                if doshrink:
                    for j in one2np1:
                        sim[j] = sim[0] + sigma*(sim[j] - sim[0])
                        fsim[j] = apply(func,(sim[j],)+args)
                    funcalls = funcalls + N

        ind = numarray.argsort(fsim)
        sim = numarray.take(sim,ind,0)
        fsim = numarray.take(fsim,ind)
        iterations = iterations + 1

    x = sim[0]
    fval = min(fsim)
    warnflag = 0

    if funcalls >= maxfun:
        warnflag = 1
        if printmessg:
            print "Warning: Maximum number of function evaluations has been exceeded."
    elif iterations >= maxiter:
        warnflag = 2
        if printmessg:
            print "Warning: Maximum number of iterations has been exceeded"
    else:
        if printmessg:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % iterations
            print "         Function evaluations: %d" % funcalls

    if fulloutput:
        return x, fval, warnflag
    else:        
        return x


def zoom(a_lo, a_hi):
    pass

    

def line_search(f, fprime, xk, pk, gfk, args=(), c1=1e-4, c2=0.9, amax=50):
    """alpha, fc, gc = line_search(f, xk, pk, gfk, args=(), c1=1e-4, c2=0.9, amax=1)

    minimize the function f(xk+alpha pk) using the line search algorithm of
    Wright and Nocedal in 'numarrayerical Optimization', 1999, pg. 59-60
    """

    fc = 0
    gc = 0
    alpha0 = 1.0
    phi0  = apply(f,(xk,)+args)
    phi_a0 = apply(f,(xk+alpha0*pk,)+args)
    fc = fc + 2
    derphi0 = numarray.dot(gfk,pk)
    derphi_a0 = numarray.dot(apply(fprime,(xk+alpha0*pk,)+args),pk)
    gc = gc + 1

    # check to see if alpha0 = 1 satisfies Strong Wolfe conditions.
    if (phi_a0 <= phi0 + c1*alpha0*derphi0) \
       and (numarray.absolute(derphi_a0) <= c2*numarray.absolute(derphi0)):
        return alpha0, fc, gc

    alpha0 = 0
    alpha1 = 1
    phi_a1 = phi_a0
    phi_a0 = phi0

    i = 1
    while 1:
        if (phi_a1 > phi0 + c1*alpha1*derphi0) or \
           ((phi_a1 >= phi_a0) and (i > 1)):
            return zoom(alpha0, alpha1)

        derphi_a1 = numarray.dot(apply(fprime,(xk+alpha1*pk,)+args),pk)
        gc = gc + 1
        if (numarray.absolute(derphi_a1) <= -c2*derphi0):
            return alpha1

        if (derphi_a1 >= 0):
            return zoom(alpha1, alpha0)

        alpha2 = (amax-alpha1)*0.25 + alpha1
        i = i + 1
        alpha0 = alpha1
        alpha1 = alpha2
        phi_a0 = phi_a1
        phi_a1 = apply(f,(xk+alpha1*pk,)+args)

    

def line_search_BFGS(f, xk, pk, gfk, args=(), c1=1e-4, alpha0=1):
    """alpha, fc, gc = line_search(f, xk, pk, gfk, args=(), c1=1e-4, alpha0=1)

    minimize over alpha, the function f(xk+alpha pk) using the interpolation
    algorithm (Armiijo backtracking) as suggested by
    Wright and Nocedal in 'numarrayerical Optimization', 1999, pg. 56-57
    """

    fc = 0
    phi0 = apply(f,(xk,)+args)               # compute f(xk)
    phi_a0 = apply(f,(xk+alpha0*pk,)+args)     # compute f
    fc = fc + 2
    derphi0 = numarray.dot(gfk,pk)

    if (phi_a0 <= phi0 + c1*alpha0*derphi0):
        return alpha0, fc, 0

    # Otherwise compute the minimizer of a quadratic interpolant:

    alpha1 = -(derphi0) * alpha0**2 / 2.0 / (phi_a0 - phi0 - derphi0 * alpha0)
    phi_a1 = apply(f,(xk+alpha1*pk,)+args)
    fc = fc + 1

    if (phi_a1 <= phi0 + c1*alpha1*derphi0):
        return alpha1, fc, 0

    # Otherwise loop with cubic interpolation until we find an alpha which satifies
    #  the first Wolfe condition (since we are backtracking, we will assume that
    #  the value of alpha is not too small and satisfies the second condition.

    while 1:       # we are assuming pk is a descent direction
        factor = alpha0**2 * alpha1**2 * (alpha1-alpha0)
        a = alpha0**2 * (phi_a1 - phi0 - derphi0*alpha1) - \
            alpha1**2 * (phi_a0 - phi0 - derphi0*alpha0)
        a = a / factor
        b = -alpha0**3 * (phi_a1 - phi0 - derphi0*alpha1) + \
            alpha1**3 * (phi_a0 - phi0 - derphi0*alpha0)
        b = b / factor

        alpha2 = (-b + numarray.sqrt(numarray.absolute(b**2 - 3 * a * derphi0))) / (3.0*a)
        phi_a2 = apply(f,(xk+alpha2*pk,)+args)
        fc = fc + 1

        if (phi_a2 <= phi0 + c1*alpha2*derphi0):
            return alpha2, fc, 0

        if (alpha1 - alpha2) > alpha1 / 2.0 or (1 - alpha2/alpha1) < 0.96:
            alpha2 = alpha1 / 2.0

        alpha0 = alpha1
        alpha1 = alpha2
        phi_a0 = phi_a1
        phi_a1 = phi_a2

epsilon = 1e-8

def approx_fprime(xk,f,*args):
    f0 = apply(f,(xk,)+args)
    grad = numarray.zeros((len(xk),),'d')
    ei = numarray.zeros((len(xk),),'d')
    for k in range(len(xk)):
        ei[k] = 1.0
        grad[k] = (apply(f,(xk+epsilon*ei,)+args) - f0)/epsilon
        ei[k] = 0.0
    return grad

def approx_fhess_p(x0,p,fprime,*args):
    f2 = apply(fprime,(x0+epsilon*p,)+args)
    f1 = apply(fprime,(x0,)+args)
    return (f2 - f1)/epsilon


def fminBFGS(f, x0, fprime=None, args=(), avegtol=1e-5, maxiter=None, fulloutput=0, printmessg=1):
    """xopt = fminBFGS(f, x0, fprime=None, args=(), avegtol=1e-5, maxiter=None, fulloutput=0, printmessg=1)

    Optimize the function, f, whose gradient is given by fprime using the
    quasi-Newton method of Broyden, Fletcher, Goldfarb, and Shanno (BFGS)
    See Wright, and Nocedal 'numarrayerical Optimization', 1999, pg. 198.
    """

    app_fprime = 0
    if fprime is None:
        app_fprime = 1

    x0 = numarray.asarray(x0)
    if maxiter is None:
        maxiter = len(x0)*200
    func_calls = 0
    grad_calls = 0
    k = 0
    N = len(x0)
    gtol = N*avegtol
    I = numarray.identity(N)
    Hk = I

    if app_fprime:
        gfk = apply(approx_fprime,(x0,f)+args)
        func_calls = func_calls + len(x0) + 1
    else:
        gfk = apply(fprime,(x0,)+args)
        grad_calls = grad_calls + 1
    xk = x0
    sk = [2*gtol]
    while (numarray.add.reduce(numarray.absolute(gfk)) > gtol) and (k < maxiter):
        pk = -numarray.dot(Hk,gfk)
        alpha_k, fc, gc = line_search_BFGS(f,xk,pk,gfk,args)
        func_calls = func_calls + fc
        xkp1 = xk + alpha_k * pk
        sk = xkp1 - xk
        xk = xkp1
        if app_fprime:
            gfkp1 = apply(approx_fprime,(xkp1,f)+args)
            func_calls = func_calls + gc + len(x0) + 1
        else:
            gfkp1 = apply(fprime,(xkp1,)+args)
            grad_calls = grad_calls + gc + 1

        yk = gfkp1 - gfk
        k = k + 1

        rhok = 1 / numarray.dot(yk,sk)
        A1 = I - sk[:,numarray.NewAxis] * yk[numarray.NewAxis,:] * rhok
        A2 = I - yk[:,numarray.NewAxis] * sk[numarray.NewAxis,:] * rhok
        Hk = numarray.dot(A1,numarray.dot(Hk,A2)) + rhok * sk[:,numarray.NewAxis] * sk[numarray.NewAxis,:]
        gfk = gfkp1


    if printmessg or fulloutput:
        fval = apply(f,(xk,)+args)
    if k >= maxiter:
        warnflag = 1
        if printmessg:
            print "Warning: Maximum number of iterations has been exceeded"
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls
            print "         Gradient evaluations: %d" % grad_calls
    else:
        warnflag = 0
        if printmessg:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % func_calls
            print "         Gradient evaluations: %d" % grad_calls

    if fulloutput:
        return xk, fval, func_calls, grad_calls, warnflag
    else:        
        return xk


def fminNCG(f, x0, fprime, fhess_p=None, fhess=None, args=(), avextol=1e-5, maxiter=None, fulloutput=0, printmessg=1):
    """xopt = fminNCG(f, x0, fprime, fhess_p=None, fhess=None, args=(), avextol=1e-5, maxiter=None, fulloutput=0, printmessg=1)

    Optimize the function, f, whose gradient is given by fprime using the
    Newton-CG method.  fhess_p must compute the hessian times an arbitrary
    vector. If it is not given, finite-differences on fprime are used to
    compute it. See Wright, and Nocedal 'numarrayerical Optimization', 1999,
    pg. 140.
    """

    x0 = numarray.asarray(x0)
    fcalls = 0
    gcalls = 0
    hcalls = 0
    approx_hessp = 0
    if fhess_p is None and fhess is None:    # Define hessian product
        approx_hessp = 1
    
    xtol = len(x0)*avextol
    update = [2*xtol]
    xk = x0
    k = 0
    while (numarray.add.reduce(numarray.absolute(update)) > xtol) and (k < maxiter):
        # Compute a search direction pk by applying the CG method to
        #  del2 f(xk) p = - grad f(xk) starting from 0.
        b = -apply(fprime,(xk,)+args)
        gcalls = gcalls + 1
        maggrad = numarray.add.reduce(numarray.absolute(b))
        eta = min([0.5,numarray.sqrt(maggrad)])
        termcond = eta * maggrad
        xsupi = numarray.zeros((len(x0),))
        ri = -b
        psupi = -ri
        i = 0
        dri0 = numarray.dot(ri,ri)

        if fhess is not None:               # you want to compute hessian once.
            A = apply(fhess,(xk,)+args)
            hcalls = hcalls + 1

        while numarray.add.reduce(numarray.absolute(ri)) > termcond:
            if fhess is None:
                if approx_hessp:
                    Ap = apply(approx_fhess_p,(xk,psupi,fprime)+args)
                    gcalls = gcalls + 2
                else:
                    Ap = apply(fhess_p,(xk,psupi)+args)
                    hcalls = hcalls + 1
            else:
                # Ap = numarray.dot(A,psupi)
                Ap = numarray.matrixmultiply(A,psupi)
            # check curvature
            curv = numarray.dot(psupi,Ap)
            if (curv <= 0):
                if (i > 0):
                    break
                else:
                    xsupi = xsupi + dri0/curv * psupi
                    break
            alphai = dri0 / curv
            xsupi = xsupi + alphai * psupi
            ri = ri + alphai * Ap
            dri1 = numarray.dot(ri,ri)
            betai = dri1 / dri0
            psupi = -ri + betai * psupi
            i = i + 1
            dri0 = dri1          # update numarray.dot(ri,ri) for next time.
    
        pk = xsupi  # search direction is solution to system.
        gfk = -b    # gradient at xk
        alphak, fc, gc = line_search_BFGS(f,xk,pk,gfk,args)
        fcalls = fcalls + fc
        gcalls = gcalls + gc

        update = alphak * pk
        xk = xk + update
        k = k + 1

    if printmessg or fulloutput:
        fval = apply(f,(xk,)+args)
    if k >= maxiter:
        warnflag = 1
        if printmessg:
            print "Warning: Maximum number of iterations has been exceeded"
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % fcalls
            print "         Gradient evaluations: %d" % gcalls
            print "         Hessian evaluations: %d" % hcalls
    else:
        warnflag = 0
        if printmessg:
            print "Optimization terminated successfully."
            print "         Current function value: %f" % fval
            print "         Iterations: %d" % k
            print "         Function evaluations: %d" % fcalls
            print "         Gradient evaluations: %d" % gcalls
            print "         Hessian evaluations: %d" % hcalls
            
    if fulloutput:
        return xk, fval, fcalls, gcalls, hcalls, warnflag
    else:        
        return xk

    

if __name__ == "__main__":
    import string
    import time

    
    times = []
    algor = []
    x0 = [0.8,1.2,0.7]
    start = time.time()
    x = fmin(rosen,x0)
    print x
    times.append(time.time() - start)
    algor.append('Nelder-Mead Simplex\t')

    start = time.time()
    x = fminBFGS(rosen, x0, fprime=rosen_der, maxiter=80)
    print x
    times.append(time.time() - start)
    algor.append('BFGS Quasi-Newton\t')

    start = time.time()
    x = fminBFGS(rosen, x0, avegtol=1e-4, maxiter=100)
    print x
    times.append(time.time() - start)
    algor.append('BFGS without gradient\t')


    start = time.time()
    x = fminNCG(rosen, x0, rosen_der, fhess_p=rosen3_hess_p, maxiter=80)
    print x
    times.append(time.time() - start)
    algor.append('Newton-CG with hessian product')
    

    start = time.time()
    x = fminNCG(rosen, x0, rosen_der, fhess=rosen3_hess, maxiter=80)
    print x
    times.append(time.time() - start)
    algor.append('Newton-CG with full hessian')

    print "\nMinimizing the Rosenbrock function of order 3\n"
    print " Algorithm \t\t\t       Seconds"
    print "===========\t\t\t      ========="
    for k in range(len(algor)):
        print algor[k], "\t -- ", times[k]

    times = []
    algor=[]
    x0 = [1.,]
    start = time.time()
    x = fmin(para,x0)
    print x
    times.append(time.time() - start)
    algor.append('Nelder-Mead Simplex\t')

    start = time.time()
    x = fminBFGS(para, x0, fprime=para_der, maxiter=80)
    print x
    times.append(time.time() - start)
    algor.append('BFGS Quasi-Newton\t')

    start = time.time()
    x = fminBFGS(para, x0, avegtol=1e-4, maxiter=100)
    print x
    times.append(time.time() - start)
    algor.append('BFGS without gradient\t')


    start = time.time()
    x = fminNCG(para, x0, para_der, fhess_p=para_hess_p, maxiter=80)
    print x
    times.append(time.time() - start)
    algor.append('Newton-CG with hessian product')
    

    start = time.time()
    x = fminNCG(para, x0, para_der, fhess=para_hess, maxiter=80)
    print x
    times.append(time.time() - start)
    algor.append('Newton-CG with full hessian')

    print "\nMinimizing x^2\n"
    print " Algorithm \t\t\t       Seconds"
    print "===========\t\t\t      ========="
    for k in range(len(algor)):
        print algor[k], "\t -- ", times[k]
