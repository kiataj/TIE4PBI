This MATLAB script estimates fringe visibility based on the Transport of Intensity Equation (TIE) in X-ray Propagation-Based Imaging (PBI). It accounts for the imaging system’s point spread function (PSF) by incorporating the full-width at half-maximum (FWHM) of both the source and detector PSFs.

By adjusting these parameters, the script enables pre-scan optimization of the imaging geometry, helping to determine the optimal source-to-object (SOD) and source-to-detector (SDD) distances before performing a PBI scan.

Below is an example of the code output, where fringe visibilities across a range of geometries are visualized using a contour plot.

If you used this code to optimize your PBI geometry please cite this [paper](https://ieeexplore.ieee.org/abstract/document/10542137).


![Fringe visibility](images/fringe_visibility.png)
