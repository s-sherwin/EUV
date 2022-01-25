# EUV
Code for modeling, designing, and measuring EUV photomasks.

Apps \n\t Standalone MATLAB applications for reflectometry, scatterometry, and Zernike phase contrast imaging. 

Diffraction
	Reciprocal-space coordinate transforms; square-wave diffraction. 

Elemental nk
	Calculating refractive indices of materials for EUV wavelengths.

Film models
	Empirically determined film models for 3 absorbers: 131 (60nm TaN), 073 (2-layer aPSM), 074 (3-layer aPSM).
	
Fourier
	Fourier transforms.

Fresnel
	Fresnel reflection and transmission coefficients w/transfer matrix method. 

Imaging
	Partially coherent imaging.

Optimization
	Optimization algorithm: coordinate descent with golden-section search, exponentially-weighted moving average (EWMA); i.e. each coordinate 1D line search, update coordinate reducing step by factor alpha, equivalent to "forgetting factor" in EWMA.
	
PanoramicAPIFuncs
	Code to interface with Panoramic EM-Suite API, to run RCWA and FDTD simulations through EM-Suite. 

PhaseLift
	Solver for PhaseLift.

Reflectivityapp
	Outputs from App for 3 masks referenced in "Film models"

Scattering
	Fresnel transfer functions and near-field double-scattering approximation 
	
Utilities
	Misc functions