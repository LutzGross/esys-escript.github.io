esys.downunder.forwardmodels.base Package
=========================================

.. module:: esys.downunder.forwardmodels.base
   :synopsis: Base classes for forward models

Base classes for forward models

Classes
-------
* `Data`
* `FileWriter`
* `ForwardModel`
* `ForwardModelWithPotential`

.. autoclass:: esys.downunder.forwardmodels.base.Data
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.forwardmodels.base.FileWriter
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.forwardmodels.base.ForwardModel
   :members:
   :undoc-members:

   .. automethod:: __init__

.. autoclass:: esys.downunder.forwardmodels.base.ForwardModelWithPotential
   :members:
   :undoc-members:

   .. automethod:: __init__


Functions
---------
.. autofunction:: Abs
.. autofunction:: C_GeneralTensorProduct
.. autofunction:: L2
.. autofunction:: LinearSinglePDE
.. autofunction:: Lsup
.. autofunction:: NumpyToData
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
.. autofunction:: getMPIRankWorld
.. autofunction:: getMPIWorldMax
.. autofunction:: getMaxFloat
.. autofunction:: getNumpy
.. autofunction:: getRank
.. autofunction:: getShape
.. autofunction:: getTagNames
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
.. autofunction:: jump
.. autofunction:: kronecker
.. autofunction:: length
.. autofunction:: listEscriptParams
.. autofunction:: log
.. autofunction:: log10
.. autofunction:: longestEdge
.. autofunction:: makeTagMap
.. autofunction:: makeTransformation
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
.. autofunction:: reorderComponents
.. autofunction:: resolve
.. autofunction:: safeDiv
.. autofunction:: saveDataCSV
.. autofunction:: saveESD
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


