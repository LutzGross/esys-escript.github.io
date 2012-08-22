
########################################################
#
# Copyright (c) 2003-2012 by University of Queensland
# Earth Systems Science Computational Center (ESSCC)
# http://www.uq.edu.au/esscc
#
# Primary Business: Queensland, Australia
# Licensed under the Open Software License version 3.0
# http://www.opensource.org/licenses/osl-3.0.php
#
########################################################

__copyright__="""Copyright (c) 2003-2012 by University of Queensland
Earth Systems Science Computational Center (ESSCC)
http://www.uq.edu.au/esscc
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Open Software License version 3.0
http://www.opensource.org/licenses/osl-3.0.php"""
__url__="https://launchpad.net/escript-finley"

import logging

from esys.escript import *
from esys.weipa import createDataset, saveSilo

from costfunctions import SimpleCostFunction
from forwardmodels import GravityModel
from mappings import *
from minimizers import *
from regularizations import Regularization

class InversionBase(object):
    """
    Base class for inversions.
    """
    def __init__(self):
        self.logger=logging.getLogger('inv.%s'%self.__class__.__name__)
        self._solver_opts = {}
        self._solver_tol = 1e-9
        self._solver_maxiter = 200
        self._output_dir = '.'
        # use identity mapping by default
        self.mapping=ScalingMapping(1)
        self.source=None
        self.solverclass=MinimizerLBFGS

    def setOutputDirectory(self, outdir):
        self._output_dir = outdir

    def setSolverMaxIterations(self, maxiter):
        self._solver_maxiter=maxiter

    def setSolverOptions(self, **opts):
        self._solver_opts.update(**opts)

    def setSolverTolerance(self, tol):
        self._solver_tol=tol

    def setSolverClass(self, solverclass):
        self.solverclass=solverclass

    def setDataSource(self, source):
        self.source=source

    def setMapping(self, mapping):
        self.mapping=mapping

    def run(self):
        raise NotImplementedError


class GravityInversion(InversionBase):
    """
    """
    def __init__(self):
        super(GravityInversion,self).__init__()
        self.setWeights()

    def setWeights(self, mu_reg=1., mu_model=1.):
        self._mu_reg=mu_reg
        self._mu_model=mu_model

    def solverCallback(self, k, x, fx, gfx):
        fn=os.path.join(self._output_dir, 'inv.%d'%k)
        ds=createDataset(rho=self.mapping.getValue(x))
        ds.setCycleAndTime(k,k)
        ds.saveSilo(fn)
        self.logger.debug("Jreg(m) = %e"%self.regularization.getValue(x))
        self.logger.debug("f(m) = %e"%fx)

    def run(self):
        if self.source is None:
            raise ValueError("No data source set!")

        self.logger.info('Retrieving domain...')
        domain=self.source.getDomain()
        DIM=domain.getDim()
        self.logger.info("Retrieving density mask...")
        rho_mask = self.source.getDensityMask()
        self.logger.info("Retrieving gravity and standard deviation data...")
        g, sigma=self.source.getGravityAndStdDev()
        chi=safeDiv(1., sigma*sigma)
        chi=interpolate(chi, Function(domain))
        g=interpolate(g, Function(domain))
        self.logger.debug("g = %s"%g)
        self.logger.debug("sigma = %s"%sigma)
        self.logger.debug("chi = %s"%chi)
        chi=chi*kronecker(DIM)[DIM-1]
        m_ref=self.mapping.getInverse(0.)
        self.regularization=Regularization(domain, m_ref=m_ref, w0=0, w=[1]*DIM, location_of_set_m=rho_mask)
        self.forwardmodel=GravityModel(domain, chi, g)
        self.f=SimpleCostFunction(self.regularization, self.mapping, self.forwardmodel)
        self.f.setWeights(mu_reg=self._mu_reg, mu_model=self._mu_model)
        solver=self.solverclass(self.f)
        solver.setTolerance(self._solver_tol)
        solver.setMaxIterations(self._solver_maxiter)
        solver.setOptions(**self._solver_opts)
        self.logger.info("Starting solver...")
        rho_init=Scalar(0, ContinuousFunction(domain))
        m_init=self.mapping.getInverse(rho_init)

        solver.setCallback(self.solverCallback)
        args={'rho_mask':rho_mask,'g':g[DIM-1],'chi':chi[DIM-1],'sigma':sigma}
        try:
            args['rho_ref']=self.source.getReferenceDensity()
        except:
            pass
        saveSilo(os.path.join(self._output_dir, 'ref'), **args)
        solver.run(m_init)
        m_star=solver.getResult()
        self.logger.info("m* = %s"%m_star)
        rho_star=self.mapping.getValue(m_star)
        self.logger.info("rho* = %s"%rho_star)
        solver.logSummary()
        return rho_star


if __name__=="__main__":
    from esys.escript import unitsSI as U
    from datasources import *

    p={
        'PADDING_L' : 5,
        'PADDING_H' : 0.2,
        'TOLERANCE' : 1e-9,
        'MAX_ITER'  : 200,
        'VERBOSITY' : 5,
        'OUTPUT_DIR': '.',
        'LOGFILE'   : '',
        'SOURCE'    : SyntheticDataSource,
        'SOLVER_OPTS': {\
            'initialHessian' : 100
        },
        'ARGS'      : {\
            'DIM'     : 2,
            'NE'      : 40,
            'l'       : 500*U.km,
            'h'       : 60*U.km,
            'features': [\
                SmoothAnomaly(lx=50*U.km, ly=20*U.km, lz=40*U.km, x=100*U.km, y=3*U.km, depth=25*U.km, rho_inner=200., rho_outer=1e-6),
                SmoothAnomaly(lx=50*U.km, ly=20*U.km, lz=40*U.km, x=400*U.km, y=1*U.km, depth=40*U.km, rho_inner=-200, rho_outer=1e-6)
            ]
        }
    }

    # 1..5 -> 50..10
    loglevel=60-10*max(1, min(5, p['VERBOSITY']))
    formatter=logging.Formatter('[%(name)s] \033[1;30m%(message)s\033[0m')
    logger=logging.getLogger('inv')
    logger.setLevel(loglevel)
    handler=logging.StreamHandler()
    handler.setFormatter(formatter)
    handler.setLevel(loglevel)
    logger.addHandler(handler)
    if len(p['LOGFILE'].strip())>0:
        handler=logging.FileHandler(os.path.join(p['OUTPUT_DIR'],p['LOGFILE']))
        formatter=logging.Formatter('%(asctime)s - [%(name)s] %(message)s')
        handler.setFormatter(formatter)
        handler.setLevel(loglevel)
        logger.addHandler(handler)

    if logger.isEnabledFor(logging.DEBUG):
        for k in sorted(p): logger.debug("%s = %s"%(k,p[k]))
    source=p['SOURCE'](**p['ARGS'])
    source.setPadding(p['PADDING_L'], p['PADDING_H'])
    inv=GravityInversion()
    inv.setDataSource(source)
    inv.setSolverTolerance(p['TOLERANCE'])
    inv.setSolverMaxIterations(p['MAX_ITER'])
    inv.setSolverOptions(**p['SOLVER_OPTS'])
    if p.has_key('MU'):
        mu=p['MU']
    else:
        logger.info('Generating domain...')
        x=source.getDomain().getX()
        l0=sup(x[0])-inf(x[0])
        l1=sup(x[1])-inf(x[1])
        l=max(l0,l1)
        G=6.6742e-11
        mu=0.5*(l**2*G)**2
        logger.debug("MU = %s"%mu)

    inv.setWeights(mu_reg=mu)
    #mapping=BoundedRangeMapping(-200, 200)
    #inv.setMapping(mapping)
    rho_new=inv.run()

