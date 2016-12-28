Utilities for flow analysis using Proper Orthogonal Decomposition
=================================================================

Utilities designed for OpenFOAM, for saving field values at certain timesteps as plain text and interpolating fields at the pre-post-processing stages.

Requirements:
* [OpenFOAM](http://www.openfoam.org)
* Git (actually, this one is optional)

How to get:
```
git clone https://github.com/engenegr/foamPODAnalysis.git

```

then if you need setPODFields utility
```
cd setPODFields
wmake

```

If you don't have Git installed, you can use the ZIP option shown in the project page: https://github.com/engenegr/foamPODAnalysis

The reconstructFromComponents utility was developed using code (https://github.com/wyldckat/reconstruct-interpolate-fields) has been written by Bruno Santos (wyldckat@github working at [blueCAPE Lda](http://www.bluecape.com.pt)).