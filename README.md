<div align="left">
    <img width=100 src="https://c1.staticflickr.com/9/8543/8982514796_309c8553b0_b.jpg" alt="ascidian" title="ascidian"</img>
</div>
# ascidian: parallel 2D wind driven ocean circulation wave model

**ascidian** is intended to simulate regional-scale wind driven ocean circulation by simulating:
- 2D ocean currents, and
- ocean waves dynamic

The model will be merged within **badlands** to compute seabed transport over geological time-scale (several thousands of years).

In comparison to other approaches used for geological problem, the model solves the ocean dynamic in a more rigourous way. The method consists in defining several seasonal (fair-weather / storms) inputs and to compute for each of these scenarios the associated stresses applied to the seabed. The implications on the seabed evolution are computed on a yearly to decadal basis, assuming the changes on the seabed morphology are small and have minor effects on ocean/wave dynamic. This assumption is considered valid for the spatial scale that we intend to use.  

The model is built around 2 codes:
- a simplified version of the spectral wave model **SWAN** (**S**imulating **WA**ve **N**earshore)
- a updated version of the tide and wind driven circulation model (**TAWIC**) 

## SWAN

The model ([code](http://swanmodel.sourceforge.net/download/download.htm)) is a spectral wave model developed at the Delft University of Technology, The Netherlands.  SWAN models the energy contained in waves as they travel over the ocean surface towards the shore. 

## TAWIC

The model solves the nonlinear two-dimensional hydrodynamic equations by a finite differencing method on a constructed zigzag boundary. 

Reference:

- Caviglia F. J. & Dragani W. C., 1996. *An improved 2-d finite-difference circulation model for tide- and wind-induced flows*. Computers & Geosciences, 22(10), 1083-1096.  

## Input/output

**ascidian** uses a XmL input file and output are generate in hdf5. The results are directly visualised in Paraview or MayaVI. 

## Example


