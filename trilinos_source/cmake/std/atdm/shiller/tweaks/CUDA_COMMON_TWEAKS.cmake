# Disable test that runs over 30 min currently (#2446)
ATDM_SET_ENABLE(PanzerAdaptersSTK_MixedPoissonExample-ConvTest-Hex-Order-3_DISABLE ON)

# Disable test that consistantly times out at 10 minutes (#3579)
ATDM_SET_ENABLE(PanzerAdaptersSTK_CurlLaplacianExample-ConvTest-Quad-Order-4_DISABLE ON)
