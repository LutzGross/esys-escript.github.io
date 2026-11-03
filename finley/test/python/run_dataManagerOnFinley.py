
##############################################################################
#
# Copyright (c) 2003-2026 by the esys.escript Group
# https://github.com/LutzGross/esys-escript.github.io
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# See CREDITS file for contributors and development history
#
##############################################################################


__copyright__="""Copyright (c) 2003-2026 by the esys.escript Group
https://github.com/LutzGross/esys-escript.github.io
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
Tests of the RESTART export of DataManager on finley (issue #137).
"""

import esys.escriptcore.utestselect as unittest
from esys.escriptcore.testing import *
from esys.escript import *
from esys.finley import Rectangle
import os
import pickle
import shutil

try:
     FINLEY_WORKDIR=os.environ['FINLEY_WORKDIR']
except KeyError:
     FINLEY_WORKDIR=os.path.join(os.getcwd(),"finley/test/")

class Test_DataManagerRestart(unittest.TestCase):

     def setUp(self):
        self.domain=Rectangle(n0=4, n1=4)
        self.workdir=os.path.join(FINLEY_WORKDIR, "datamanager_"+self._testMethodName)
        if getMPIRankWorld()==0 and os.path.isdir(self.workdir):
            shutil.rmtree(self.workdir)
        MPIBarrierWorld()

     def tearDown(self):
        del self.domain

     def test_unpicklable_value_is_named(self):
        dm=DataManager(formats=[DataManager.RESTART], work_dir=self.workdir,
                       do_restart=False)
        u=Scalar(1., Solution(self.domain))
        dm.addData(u=u, dt=0.1, mydomain=self.domain)
        with self.assertRaises(pickle.PicklingError) as cm:
            dm.export()
        msg=str(cm.exception)
        self.assertIn("'mydomain'", msg)
        self.assertIn("FinleyDomain", msg)
        self.assertNotIn("'dt'", msg)
        self.assertFalse(os.path.isdir(os.path.join(self.workdir, "restart_0000")),
                         "restart directory written despite the error")

     @unittest.skipIf(not loadIsConfigured(), "load not configured")
     def test_picklable_values_are_saved(self):
        dm=DataManager(formats=[DataManager.RESTART], work_dir=self.workdir,
                       do_restart=False)
        u=Scalar(1., Solution(self.domain))
        dm.addData(u=u, dt=0.1, label="run1")
        dm.export()
        restartdir=os.path.join(self.workdir, "restart_0000")
        stamp_file=os.path.join(restartdir, "stamp.%d"%getMPIRankWorld())
        with open(stamp_file, "rb") as f:
            stamp=pickle.load(f)
        self.assertEqual(stamp['dt'], 0.1)
        self.assertEqual(stamp['label'], "run1")
        # with more than one rank the file name carries the rank, u.h5.0000
        self.assertTrue(any(f.startswith("u.h5") for f in os.listdir(restartdir)))

     # work_dir is a subdirectory, never the current directory: restart
     # folders used to be looked up relative to the current directory
     @unittest.skipIf(not loadIsConfigured(), "load not configured")
     def test_restart_from_work_dir(self):
        dm=DataManager(formats=[DataManager.RESTART], work_dir=self.workdir,
                       do_restart=False)
        x=self.domain.getX()
        u=x[0]+2*x[1]
        dm.addData(u=u, dt=0.1)
        dm.export()
        dm2=DataManager(formats=[DataManager.RESTART], work_dir=self.workdir)
        self.assertTrue(dm2.hasData(), "restart files in work_dir not found")
        self.assertEqual(dm2.getValue("dt"), 0.1)
        u2=dm2.getValue("u")
        # the restart loads a new domain, so rebuild u on it to compare
        domain2=dm2.getDomain()
        self.assertEqual(u2.getDomain(), domain2)
        x2=domain2.getX()
        u_ref=interpolate(x2[0]+2*x2[1], u2.getFunctionSpace())
        self.assertTrue(Lsup(u2-u_ref) <= 1e-12*Lsup(u_ref), "restored data differ")
        self.assertEqual(dm2.getCycle(), 0)

     @unittest.skipIf(not loadIsConfigured(), "load not configured")
     def test_no_restart_removes_old_folders(self):
        dm=DataManager(formats=[DataManager.RESTART], work_dir=self.workdir,
                       do_restart=False)
        dm.addData(u=Scalar(1., Solution(self.domain)), dt=0.1)
        dm.export()
        restartdir=os.path.join(self.workdir, "restart_0000")
        self.assertTrue(os.path.isdir(restartdir))
        dm2=DataManager(formats=[DataManager.RESTART], work_dir=self.workdir,
                        do_restart=False)
        self.assertFalse(dm2.hasData())
        self.assertFalse(os.path.isdir(restartdir), "old restart folder not removed")

if __name__ == '__main__':
    run_tests(__name__, exit_on_failure=True)
