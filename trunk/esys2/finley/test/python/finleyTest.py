import sys
import unittest
import os

esys_root=os.getenv('ESYS_ROOT')
sys.path.append(esys_root+'/finley/lib')
sys.path.append(esys_root+'/escript/lib')
sys.path.append(esys_root+'/escript/py_src')

import escript
import finley

##      mesh=finley.ReadMesh('myMesh')
class MeshTestCase(unittest.TestCase):
  def testMeshIO(self):
      """Test reading and writing finley meshes."""
      print
      mesh=finley.Brick(1,1,1,1,1.,1.,1.,1,1,1,1,1)
      fileName='Junk.msh'
      mesh.write(fileName)
      mesh2=finley.ReadMesh(fileName)
  def testNonExistantMeshRead(self):
      """Test exception generated for attempt to read non existant file."""
      print
      try:
        mesh2=finley.ReadMesh("NonExistantFilename.msh")
        self.failIf(False,'Failed non existant mesh file exception test.')
      except Exception, e:
        print e
  def testSystemMatrix(self):
    """Test System Matrix."""
    print
    mesh=finley.Brick()
    systemMatrix=finley.MeshAdapter(mesh)

suite=unittest.TestSuite()
suite.addTest(unittest.makeSuite(MeshTestCase))
unittest.TextTestRunner(verbosity=2).run(suite)
