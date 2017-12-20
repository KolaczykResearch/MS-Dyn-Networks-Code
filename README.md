
README

    Xinyu Kang
    xkang@bu.edu
    December 2017



This package contains MatLab(TM) scripts that implement the recursive partition based multi scale dynamic causal network model, introduced in  
   
    Dynamic Networks with Multi-scale Temporal Structure
    Xinyu Kang, Apratim Ganguly, Eric Kolaczyk

===== Implementation ==================================================================
The code is implemented using Matlab and depends on the econometrics toolbox and the cvx toolbox (http://cvxr.com/cvx/).

To implement the code, please run either RDPmain.m or RPmain.m to see a minimum working example with the recursive dyadic partitioning 
or the recursive partitioning respectively. This will call the other functions in this page.

===== Scripts ========================================================================+
- compute_lambda.m:	Compute the empirical penalty for the group lasso type of optimization.

- covar_design.m:	Reshape the multi-variate time series data into a form that is suitable for the group-lasso regression.

- data_gen.m:		Generate multi-variate time series data of length T using coeffcient matrix A

- GetChangepoints.m:	Extracts an optimal segmentation from probability arrays and mask calculated with TriangleMDN.m.

- g_lasso.m		Solve group lasso problem via ADMM minimize 1/2*|| Ax - b ||_2^2 + \lambda sum(norm(x_i))

- importfile.m		Import MEG data. (Data currently not available)

- MEGanalysis.m:		Brain image analysis (This replicates our applications in the paper. The data is confidential and is not available at this moment)

- MEG_breakpoints.m:    Detect change points of the MEG networks

- RDPmain.m:     	Simulation study using the recursive dyadic partitioning.

- RecoverCoef.m:	recover the piecewise stationary AR coefficients.

- RPmain.m:		Simulation study using recursive partitioning.

- plotbp.m:		Generate figure 3 in the multi-scale network paper.

- rdp_ts.m:		Fitting the model using recursive dyadic partitioning.

- shadedErrorBar.m	Makes a 2-d line plot with a shaded error bar. Copyright (c) 2016, Michele Giugliano.

- shadeplot.m		Generate figure 4 and figure 5 in the “Dynamic Networks with Multi-scale Temporal Structure” paper.

- TriangleMDN.m:	Performs a recursive, pyramidal series of local (in scale and position) penalized likelihood optimizations.
