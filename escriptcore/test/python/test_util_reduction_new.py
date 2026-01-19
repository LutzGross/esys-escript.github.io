
##############################################################################
#
# Copyright (c) 2003-2018 by The University of Queensland
# http://www.uq.edu.au
#
# Primary Business: Queensland, Australia
# Licensed under the Apache License, version 2.0
# http://www.apache.org/licenses/LICENSE-2.0
#
# Development until 2012 by Earth Systems Science Computational Center (ESSCC)
# Development 2012-2013 by School of Earth Sciences
# Development from 2014 by Centre for Geoscience Computing (GeoComp)
#
##############################################################################


__copyright__="""Copyright (c) 2003-2018 by The University of Queensland
http://www.uq.edu.au
Primary Business: Queensland, Australia"""
__license__="""Licensed under the Apache License, version 2.0
http://www.apache.org/licenses/LICENSE-2.0"""
__url__="https://github.com/LutzGross/esys-escript.github.io"

"""
test for util operations for reduction operations with tagged data

:remark: use see `test_util`
:var __author__: name of author
:var __copyright__: copyrights
:var __license__: licence agreement
:var __url__: url entry point on documentation
:var __version__: version
:var __date__: date of the version
"""

__author__="Joel Fenwick, joelfenwick@uq.edu.au"

import esys.escriptcore.utestselect as unittest
import numpy
from esys.escript import *
from test_util_base import Test_util_base, Test_util_values

def zero_to_nan(obj):
    f=1./obj
    return f/f

def zero_to_inf(obj):
    return 1./obj

class Test_util_reduction_new(Test_util_base, Test_util_values):        
    def test_Lsup_new(self):
        supportcplx=True
        opstring="Lsup(a)"
        misccheck="isinstance(res,float)"
        oraclecheck="abs(ref).max()"
        opname="Lsup"
        update1="max(abs(r).max(),abs(r2).max())"
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)

    def test_sup_new(self):
        supportcplx=False
        opstring="sup(a)"
        misccheck="isinstance(res,float)"
        oraclecheck="ref.max()"
        opname="sup"
        update1="max(r.max(), r2.max())"
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)

    def test_inf_new(self):
        supportcplx=False
        opstring="inf(a)"
        misccheck="isinstance(res,float)"
        oraclecheck="ref.min()"
        opname="inf"
        update1="min(r.min(),r2.min())"
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1)
        
    @unittest.skipIf(not hasFeature('NAN_CHECK'), "test only fires if NAN_CHECK is enabled")
    def test_hasNaN(self):
        # Need to check for hasNaN as well
        supportcplx=True
        opstring="a.hasNaN()"
        misccheck=None
        oraclecheck="0 in ref"
        opname="hasNaN"
        update1="bool(numpy.isnan(r).max()) or bool(numpy.isnan(r2).max())"       # numpy.bool_ is not bool
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1, input_trans=zero_to_nan, no_scalars=True)   
        
    # It would be a bit tricky to reformulate this into the new form 
    # This will not test all possible type combinations 
    @unittest.skipIf(not hasFeature('NAN_CHECK'), "test only fires if NAN_CHECK is enabled")
    def test_NaNReduction_constData_rank4(self):
        oarg=Data(numpy.array([[[[0.50544713768476202, 0.96922321849050874, -0.81524480218696649, -0.36499730379849193], 
[-0.48131882706974372, 0.026812357207576465, 0.090903267401989618, -0.24742363369877829], [-0.51631372893805438, 
0.30410275437953183, -0.75149566289642533, -0.19930300338453599]], [[0.82034878499482788, -0.70904661587698792, 
-0.27637223434426073, -0.34818734117560401], [0.11686048779802416, -0.76746266142163178, -0.75578186306174833, 
0.14509316330390232], [0.1590050723141736, 0.69684384552537937, -0.58747105640080832, -0.28640840371441523]]], 
[[[0.14956532194045669, 0.081514192262221119, 0.32061383569406399, -0.2444346881437609], [0.79564139071785278, 
-0.5456680167461434, 0.24722978802719742, 0.28286130725068315], [0.10385207763921711, -0.064749181840278336, 
0.21325254547672734, -0.71875644540473838]], [[0.58552496009870802, 0.35472373485671338, -0.18411162994671826, 
0.71609038134967773], [-0.20966804574945064, -0.49286619989346314, 0.85116051808632553, -0.94417114370961075], 
[-0.40434528979823714, 0.62250343758157611, 0.64860074098639742, 0.0043146814280992096]]], [[[-0.14242849200713259, 
0.42551908502898095, 0.7691157770973962, -0.37595641162856674], [0.026655444032149589, -0.82186407521644167, 
0.40285091480648783, -0.53328831035315982], [-0.12887729257054481, 0.75610663428133451, 0.022049613835531723, 
0.59949338706293043]], [[-0.34506254315071772, 0.019719877473602043, 0.10216765908478709, 0.022681548062032153], 
[0.2228614880408597, 0.26944547311401901, -0.10122095357202965, -0.51019076850180589], [-0.081439546799124463, 
0.18829632566943544, 0.12366885442775377, 0]]]]),self.functionspace)
        arg=1/oarg       #will get us an inf
        arg=arg/arg     #will give a NaN in the last position, yes we could have just sqrt(arg) but I wanted last pos
        self.assertTrue(numpy.isnan(sup(arg)),"wrong result")
        self.assertTrue(numpy.isnan(inf(arg)),"wrong result")
        self.assertTrue(numpy.isnan(Lsup(arg)),"wrong result")
        arg=(1+0j)/oarg
        arg=arg/arg     #will give a NaN in the last position, yes we could have just sqrt(arg) but I wanted last pos
        self.assertRaises(RuntimeError,sup, arg)
        self.assertRaises(RuntimeError,inf, arg)
        self.assertTrue(numpy.isnan(Lsup(arg)),"wrong result")
        # Now testing tagged
        arg.tag()
        self.assertRaises(RuntimeError,sup, arg)
        self.assertRaises(RuntimeError,inf, arg)
        self.assertTrue(numpy.isnan(Lsup(arg)),"wrong result")        
                
    @unittest.skipIf(not hasFeature('NAN_CHECK'), "test only fires if NAN_CHECK is enabled")
    def test_NaNReduction_expandedData_rank4(self):
        oarg=Data(numpy.array([[[[0.50544713768476202, 0.96922321849050874, -0.81524480218696649, -0.36499730379849193], 
[-0.48131882706974372, 0.026812357207576465, 0.090903267401989618, -0.24742363369877829], [-0.51631372893805438, 
0.30410275437953183, -0.75149566289642533, -0.19930300338453599]], [[0.82034878499482788, -0.70904661587698792, 
-0.27637223434426073, -0.34818734117560401], [0.11686048779802416, -0.76746266142163178, -0.75578186306174833, 
0.14509316330390232], [0.1590050723141736, 0.69684384552537937, -0.58747105640080832, -0.28640840371441523]]], 
[[[0.14956532194045669, 0.081514192262221119, 0.32061383569406399, -0.2444346881437609], [0.79564139071785278, 
-0.5456680167461434, 0.24722978802719742, 0.28286130725068315], [0.10385207763921711, -0.064749181840278336, 
0.21325254547672734, -0.71875644540473838]], [[0.58552496009870802, 0.35472373485671338, -0.18411162994671826, 
0.71609038134967773], [-0.20966804574945064, -0.49286619989346314, 0.85116051808632553, -0.94417114370961075], 
[-0.40434528979823714, 0.62250343758157611, 0.64860074098639742, 0.0043146814280992096]]], [[[-0.14242849200713259, 
0.42551908502898095, 0.7691157770973962, -0.37595641162856674], [0.026655444032149589, -0.82186407521644167, 
0.40285091480648783, -0.53328831035315982], [-0.12887729257054481, 0.75610663428133451, 0.022049613835531723, 
0.59949338706293043]], [[-0.34506254315071772, 0.019719877473602043, 0.10216765908478709, 0.022681548062032153], 
[0.2228614880408597, 0.26944547311401901, -0.10122095357202965, -0.51019076850180589], [-0.081439546799124463, 
0.18829632566943544, 0.12366885442775377, 0]]]]),self.functionspace, True)
        arg=1/oarg       #will get us an inf
        arg=arg/arg     #will give a NaN in the last position, yes we could have just sqrt(arg) but I wanted last pos
        self.assertTrue(numpy.isnan(sup(arg)),"wrong result")
        self.assertTrue(numpy.isnan(inf(arg)),"wrong result")
        self.assertTrue(numpy.isnan(Lsup(arg)),"wrong result") 
        oarg.resolve()  # to prevent autolazy and complex interfering
        arg=(1+0j)/oarg
        arg.resolve() # to prevent autolazy and complex interfering
        arg=arg/arg     #will give a NaN in the last position, yes we could have just sqrt(arg) but I wanted last pos
        self.assertRaises(RuntimeError,sup, arg)
        self.assertRaises(RuntimeError,inf, arg)
        self.assertTrue(numpy.isnan(Lsup(arg)),"wrong result")        
        
    def test_hasInf(self):
        # Need to check for hasNaN as well
        supportcplx=True
        opstring="a.hasInf()"
        misccheck=None
        oraclecheck="0 in ref"
        opname="hasInf"
        update1="bool(numpy.isinf(r).max()) or bool(numpy.isinf(r2).max())"       # numpy.bool_ is not bool
        self.generate_operation_test_batch(supportcplx, opstring, misccheck, oraclecheck, opname, update1, input_trans=zero_to_inf, no_scalars=True)   
        
    # It would be a bit tricky to reformulate this into the new form 
    # This will not test all possible type combinations 
    def test_InfReduction_constData_rank4(self):
        oarg=Data(numpy.array([[[[0.50544713768476202, 0.96922321849050874, -0.81524480218696649, -0.36499730379849193], 
[-0.48131882706974372, 0.026812357207576465, 0.090903267401989618, -0.24742363369877829], [-0.51631372893805438, 
0.30410275437953183, -0.75149566289642533, -0.19930300338453599]], [[0.82034878499482788, -0.70904661587698792, 
-0.27637223434426073, -0.34818734117560401], [0.11686048779802416, -0.76746266142163178, -0.75578186306174833, 
0.14509316330390232], [0.1590050723141736, 0.69684384552537937, -0.58747105640080832, -0.28640840371441523]]], 
[[[0.14956532194045669, 0.081514192262221119, 0.32061383569406399, -0.2444346881437609], [0.79564139071785278, 
-0.5456680167461434, 0.24722978802719742, 0.28286130725068315], [0.10385207763921711, -0.064749181840278336, 
0.21325254547672734, -0.71875644540473838]], [[0.58552496009870802, 0.35472373485671338, -0.18411162994671826, 
0.71609038134967773], [-0.20966804574945064, -0.49286619989346314, 0.85116051808632553, -0.94417114370961075], 
[-0.40434528979823714, 0.62250343758157611, 0.64860074098639742, 0.0043146814280992096]]], [[[-0.14242849200713259, 
0.42551908502898095, 0.7691157770973962, -0.37595641162856674], [0.026655444032149589, -0.82186407521644167, 
0.40285091480648783, -0.53328831035315982], [-0.12887729257054481, 0.75610663428133451, 0.022049613835531723, 
0.59949338706293043]], [[-0.34506254315071772, 0.019719877473602043, 0.10216765908478709, 0.022681548062032153], 
[0.2228614880408597, 0.26944547311401901, -0.10122095357202965, -0.51019076850180589], [-0.081439546799124463, 
0.18829632566943544, 0.12366885442775377, 0]]]]),self.functionspace)
        arg=1/oarg       #will get us an inf
        self.assertTrue(numpy.isinf(sup(arg)),"wrong result")
        self.assertTrue(numpy.isinf(Lsup(arg)),"wrong result")
        arg=(1+1j)/oarg         # Why not just 1+0j? ... because that gives inf+nanJ ... Just don't
        self.assertRaises(RuntimeError,sup, arg)
        self.assertRaises(RuntimeError,inf, arg)
        self.assertTrue(numpy.isinf(Lsup(arg)),"wrong result")
        # Now testing tagged
        arg.tag()
        self.assertRaises(RuntimeError,sup, arg)
        self.assertRaises(RuntimeError,inf, arg)
        self.assertTrue(numpy.isinf(Lsup(arg)),"wrong result")        
                
    def test_InfReduction_expandedData_rank4(self):
        # Have not actually worked on this bit yet
        oarg=Data(numpy.array([[[[0.50544713768476202, 0.96922321849050874, -0.81524480218696649, -0.36499730379849193], 
[-0.48131882706974372, 0.026812357207576465, 0.090903267401989618, -0.24742363369877829], [-0.51631372893805438, 
0.30410275437953183, -0.75149566289642533, -0.19930300338453599]], [[0.82034878499482788, -0.70904661587698792, 
-0.27637223434426073, -0.34818734117560401], [0.11686048779802416, -0.76746266142163178, -0.75578186306174833, 
0.14509316330390232], [0.1590050723141736, 0.69684384552537937, -0.58747105640080832, -0.28640840371441523]]], 
[[[0.14956532194045669, 0.081514192262221119, 0.32061383569406399, -0.2444346881437609], [0.79564139071785278, 
-0.5456680167461434, 0.24722978802719742, 0.28286130725068315], [0.10385207763921711, -0.064749181840278336, 
0.21325254547672734, -0.71875644540473838]], [[0.58552496009870802, 0.35472373485671338, -0.18411162994671826, 
0.71609038134967773], [-0.20966804574945064, -0.49286619989346314, 0.85116051808632553, -0.94417114370961075], 
[-0.40434528979823714, 0.62250343758157611, 0.64860074098639742, 0.0043146814280992096]]], [[[-0.14242849200713259, 
0.42551908502898095, 0.7691157770973962, -0.37595641162856674], [0.026655444032149589, -0.82186407521644167, 
0.40285091480648783, -0.53328831035315982], [-0.12887729257054481, 0.75610663428133451, 0.022049613835531723, 
0.59949338706293043]], [[-0.34506254315071772, 0.019719877473602043, 0.10216765908478709, 0.022681548062032153], 
[0.2228614880408597, 0.26944547311401901, -0.10122095357202965, -0.51019076850180589], [-0.081439546799124463, 
0.18829632566943544, 0.12366885442775377, 0]]]]),self.functionspace, True)
        arg=1/oarg       #will get us an inf
        self.assertTrue(numpy.isinf(sup(arg)),"wrong result")
        self.assertTrue(numpy.isinf(Lsup(arg)),"wrong result") 
        oarg.resolve()  # to prevent autolazy and complex interfering
        arg=(1+1j)/oarg
        arg.resolve() # to prevent autolazy and complex interfering
        self.assertRaises(RuntimeError,sup, arg)
        self.assertRaises(RuntimeError,inf, arg)
        self.assertTrue(numpy.isinf(Lsup(arg)),"wrong result")        
    

        
