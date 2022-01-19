esys.escript Package
====================

.. module:: esys.escript
   :synopsis: Classes and tools that form the basis of the escript system. Specific solvers and domains are found in their respective packages.

Classes and tools that form the basis of the escript system.
Specific solvers and domains are found in their respective packages.

Classes
-------
* `ContinuousDomain`
* `Data`
* `DataManager`
* `Domain`
* `Evaluator`
* `FileWriter`
* `FunctionJob`
* `FunctionSpace`
* `Internal_SplitWorld`
* `Job`
* `NonlinearPDE`
* `Operator`
* `Reducer`
* `SolverBuddy`
* `SolverOptions`
* `SplitWorld`
* `SubWorld`
* `Symbol`
* `TestDomain`
* `TransportProblem`

.. autoclass:: esys.escript.ContinuousDomain
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.Data
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.DataManager
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.Domain
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.Evaluator
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.FileWriter
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.FunctionJob
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.FunctionSpace
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.Internal_SplitWorld
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.Job
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.NonlinearPDE
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.Operator
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.Reducer
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.SolverBuddy
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.SolverOptions
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.SplitWorld
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.SubWorld
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.Symbol
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.TestDomain
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.TransportProblem
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
.. autofunction:: Lsup
.. autofunction:: MPIBarrierWorld
.. autofunction:: NcFType
.. autofunction:: NumpyToData
.. autofunction:: RandomData
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
.. autofunction:: internal_addJob
.. autofunction:: internal_addJobPerWorld
.. autofunction:: internal_addVariable
.. autofunction:: internal_buildDomains
.. autofunction:: internal_makeDataReducer
.. autofunction:: internal_makeLocalOnly
.. autofunction:: internal_makeScalarReducer
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


