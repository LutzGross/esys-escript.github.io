esys.escript.linearPDEs Package
===============================

.. module:: esys.escript.linearPDEs
   :synopsis: -

-

Classes
-------
* `ContinuousDomain`
* `Data`
* `Domain`
* `FileWriter`
* `FunctionSpace`
* `Helmholtz`
* `IllegalCoefficient`
* `IllegalCoefficientFunctionSpace`
* `IllegalCoefficientValue`
* `Internal_SplitWorld`
* `LameEquation`
* `LinearPDE`
* `LinearProblem`
* `Operator`
* `PDECoef`
* `Poisson`
* `Reducer`
* `SolverBuddy`
* `SolverOptions`
* `SubWorld`
* `TestDomain`
* `TransportPDE`
* `TransportProblem`
* `UndefinedPDEError`
* `WavePDE`

.. autoclass:: esys.escript.linearPDEs.ContinuousDomain
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.Data
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.Domain
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.FileWriter
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.FunctionSpace
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.Helmholtz
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.IllegalCoefficient
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.IllegalCoefficientFunctionSpace
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.IllegalCoefficientValue
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.Internal_SplitWorld
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.LameEquation
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.LinearPDE
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.LinearProblem
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.Operator
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.PDECoef
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.Poisson
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.Reducer
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.SolverBuddy
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.SolverOptions
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.SubWorld
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.TestDomain
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.TransportPDE
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.TransportProblem
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.UndefinedPDEError
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.escript.linearPDEs.WavePDE
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
.. autofunction:: LinearPDESystem
.. autofunction:: LinearSinglePDE
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
.. autofunction:: SingleTransportPDE
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
.. autofunction:: printParallelThreadCounts
.. autofunction:: releaseUnusedMemory
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

Packages
--------
.. toctree::


