esys.downunder.apps Package
===========================

.. module:: esys.downunder.apps
   :synopsis: Our most general domain representation. Imports submodules into its namespace

Our most general domain representation. Imports submodules into its namespace

Classes
-------
* `ContinuousDomain`
* `DCResistivityModel`
* `DCResistivityModelNoPrimary`
* `Data`
* `DataManager`
* `Domain`
* `Evaluator`
* `FileWriter`
* `FunctionJob`
* `FunctionSpace`
* `GravityModel`
* `Job`
* `LinearPDE`
* `Locator`
* `MT2DTEModel`
* `MT2DTMModel`
* `MagneticModel2D`
* `MagneticModel3D`
* `NonlinearPDE`
* `Operator`
* `PMLCondition`
* `Reducer`
* `SolverBuddy`
* `SolverOptions`
* `SonicWaveInFrequencyDomain`
* `SplitWorld`
* `SubWorld`
* `Symbol`
* `TestDomain`
* `TransportProblem`

.. autoclass:: esys.downunder.apps.ContinuousDomain
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.DCResistivityModel
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.DCResistivityModelNoPrimary
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.Data
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.DataManager
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.Domain
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.Evaluator
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.FileWriter
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.FunctionJob
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.FunctionSpace
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.GravityModel
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.Job
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.LinearPDE
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.Locator
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.MT2DTEModel
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.MT2DTMModel
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.MagneticModel2D
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.MagneticModel3D
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.NonlinearPDE
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.Operator
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.PMLCondition
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.Reducer
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.SolverBuddy
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.SolverOptions
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.SonicWaveInFrequencyDomain
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.SplitWorld
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.SubWorld
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.Symbol
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.TestDomain
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.apps.TransportProblem
   :members:
   :undoc-members:

   .. automethod:: __init__


Functions
---------
.. autofunction:: Abs
.. autofunction:: C_GeneralTensorProduct
.. autofunction:: ComplexData
.. autofunction:: ComplexScalar
.. autofunction:: ComplexTensor
.. autofunction:: ComplexTensor3
.. autofunction:: ComplexTensor4
.. autofunction:: ComplexVector
.. autofunction:: ContinuousFunction
.. autofunction:: DiracDeltaFunctions
.. autofunction:: Function
.. autofunction:: FunctionOnBoundary
.. autofunction:: FunctionOnContactOne
.. autofunction:: FunctionOnContactZero
.. autofunction:: L2
.. autofunction:: LinearSinglePDE
.. autofunction:: Lsup
.. autofunction:: MPIBarrierWorld
.. autofunction:: NcFType
.. autofunction:: NumpyToData
.. autofunction:: RandomData
.. autofunction:: ReadMesh
.. autofunction:: ReducedContinuousFunction
.. autofunction:: ReducedFunction
.. autofunction:: ReducedFunctionOnBoundary
.. autofunction:: ReducedFunctionOnContactOne
.. autofunction:: ReducedFunctionOnContactZero
.. autofunction:: ReducedSolution
.. autofunction:: Scalar
.. autofunction:: Solution
.. autofunction:: Tensor
.. autofunction:: Tensor3
.. autofunction:: Tensor4
.. autofunction:: Vector
.. autofunction:: acos
.. autofunction:: acosh
.. autofunction:: antihermitian
.. autofunction:: antisymmetric
.. autofunction:: asin
.. autofunction:: asinh
.. autofunction:: atan
.. autofunction:: atan2
.. autofunction:: atanh
.. autofunction:: boundingBox
.. autofunction:: boundingBoxEdgeLengths
.. autofunction:: canInterpolate
.. autofunction:: clip
.. autofunction:: combineData
.. autofunction:: commonDim
.. autofunction:: commonShape
.. autofunction:: condEval
.. autofunction:: convertToNumpy
.. autofunction:: cos
.. autofunction:: cosh
.. autofunction:: delay
.. autofunction:: deviatoric
.. autofunction:: diameter
.. autofunction:: div
.. autofunction:: eigenvalues
.. autofunction:: eigenvalues_and_eigenvectors
.. autofunction:: erf
.. autofunction:: escript_generalTensorProduct
.. autofunction:: escript_generalTensorTransposedProduct
.. autofunction:: escript_generalTransposedTensorProduct
.. autofunction:: escript_inverse
.. autofunction:: exp
.. autofunction:: generalTensorProduct
.. autofunction:: generalTensorTransposedProduct
.. autofunction:: generalTransposedTensorProduct
.. autofunction:: getClosestValue
.. autofunction:: getEpsilon
.. autofunction:: getEscriptParamInt
.. autofunction:: getMPIRankWorld
.. autofunction:: getMPISizeWorld
.. autofunction:: getMPIWorldMax
.. autofunction:: getMPIWorldSum
.. autofunction:: getMachinePrecision
.. autofunction:: getMaxFloat
.. autofunction:: getNumberOfThreads
.. autofunction:: getNumpy
.. autofunction:: getRank
.. autofunction:: getShape
.. autofunction:: getTagNames
.. autofunction:: getTestDomainFunctionSpace
.. autofunction:: getTotalDifferential
.. autofunction:: getVersion
.. autofunction:: gmshGeo2Msh
.. autofunction:: grad
.. autofunction:: grad_n
.. autofunction:: hasFeature
.. autofunction:: hermitian
.. autofunction:: identity
.. autofunction:: identityTensor
.. autofunction:: identityTensor4
.. autofunction:: inf
.. autofunction:: inner
.. autofunction:: insertTagNames
.. autofunction:: insertTaggedValues
.. autofunction:: integrate
.. autofunction:: interpolate
.. autofunction:: interpolateTable
.. autofunction:: inverse
.. autofunction:: isSymbol
.. autofunction:: jump
.. autofunction:: kronecker
.. autofunction:: length
.. autofunction:: listEscriptParams
.. autofunction:: listFeatures
.. autofunction:: load
.. autofunction:: loadIsConfigured
.. autofunction:: log
.. autofunction:: log10
.. autofunction:: longestEdge
.. autofunction:: makeTagMap
.. autofunction:: matchShape
.. autofunction:: matchType
.. autofunction:: matrix_mult
.. autofunction:: matrix_transposed_mult
.. autofunction:: matrixmult
.. autofunction:: maximum
.. autofunction:: maxval
.. autofunction:: meanValue
.. autofunction:: minimum
.. autofunction:: minval
.. autofunction:: mkDir
.. autofunction:: mult
.. autofunction:: negative
.. autofunction:: nonsymmetric
.. autofunction:: normalize
.. autofunction:: outer
.. autofunction:: phase
.. autofunction:: pokeDim
.. autofunction:: polarToCart
.. autofunction:: positive
.. autofunction:: pprint
.. autofunction:: pretty_print
.. autofunction:: printParallelThreadCounts
.. autofunction:: releaseUnusedMemory
.. autofunction:: removeFsFromGrad
.. autofunction:: reorderComponents
.. autofunction:: resolve
.. autofunction:: resolveGroup
.. autofunction:: runMPIProgram
.. autofunction:: safeDiv
.. autofunction:: saveDataCSV
.. autofunction:: saveESD
.. autofunction:: saveSilo
.. autofunction:: saveVTK
.. autofunction:: setEscriptParamInt
.. autofunction:: setNumberOfThreads
.. autofunction:: showEscriptParams
.. autofunction:: sign
.. autofunction:: sin
.. autofunction:: sinh
.. autofunction:: sqrt
.. autofunction:: sup
.. autofunction:: swap_axes
.. autofunction:: symbols
.. autofunction:: symmetric
.. autofunction:: tan
.. autofunction:: tanh
.. autofunction:: tensor_mult
.. autofunction:: tensor_transposed_mult
.. autofunction:: tensormult
.. autofunction:: testForZero
.. autofunction:: trace
.. autofunction:: transpose
.. autofunction:: transposed_matrix_mult
.. autofunction:: transposed_tensor_mult
.. autofunction:: unitVector
.. autofunction:: vol
.. autofunction:: whereNegative
.. autofunction:: whereNonNegative
.. autofunction:: whereNonPositive
.. autofunction:: whereNonZero
.. autofunction:: wherePositive
.. autofunction:: whereZero
.. autofunction:: zeros

Others
------
* DBLE_MAX
* EPSILON
* HAVE_SYMBOLS

Packages
--------
.. toctree::


