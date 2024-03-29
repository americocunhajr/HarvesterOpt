<img src="logo/HarvesterOpt.png" width="40%">

**HarvesterOpt: Harvester Optimizer** is a Matlab code to solve an optimization problem to maximize the energy recovered by a bistable piezo-magneto-elastic energy harvesting system. The implementation follows an educational style, to make its use very intuitive. 

This package includes the following files:
* CrossEntropyOpt_Example.m.  -- script with an one dimensional optimization with the CE method
* HarvesterOpt_DirectSearch.m -- script with a program to optimize the harveters via direct serach
* HarvesterOpt_CrossEntropy.m -- script with a program to optimize the harveters via CE method
* PiezoMagBeam_ObjFunc.m.     -- function that defines the objetive function
* PiezoMagBeam_ObjFuncNoise.m -- function that defines the objetive function with noise
* PiezoMagBeam_RHS.m          -- function that defines the harvester equations
* test01chaos.m               -- function to compute the 0-1 test for chaos
* trandn.m                    -- function to draw samples from truncated distribution
* plot_ce_animation.m         -- function to plot an animation with CE search for the optimum
* graph_contour_pnt.m         -- function to plot a contourmap and a given point
* graph_binarymap.m           -- function to plot a binary contourmap

## Software History

This code was developed to optimize the amount of energy recovered by a nonlinear bistable piezoelectric energy harvester. The simulation results associated with this optimization process are reported in the following paper:
- *A. Cunha Jr, Enhancing the performance of a bistable energy harvesting device via the cross-entropy method, Nonlinear Dynamics, 103:137-155, 2021 https://doi.org/10.1007/s11071-020-06109-0*

Preprint available at: https://hal.archives-ouvertes.fr/hal-01531845v3

## Videos
Animations illustrating the operation of the CE method can be seen in the following playlist:
- https://www.youtube.com/playlist?list=PLjNDdMKtfqVnU4QIubzpeE25szNEbEdBZ


## Author
- Americo Cunha Jr

## Citing HarvesterOpt
We ask the code users to cite the following manuscript in any publications reporting work done with **HarvesterOpt**:
- *A. Cunha Jr, Enhancing the performance of a bistable energy harvesting device via the cross-entropy method, Nonlinear Dynamics, 103:137-155, 2021 https://doi.org/10.1007/s11071-020-06109-0*

```
@article{CunhaJr2021p137,
   author  = {A. {Cunha~Jr}},
   title   = {Enhancing the performance of a bistable energy harvesting device via the cross-entropy method},
   journal = {Nonlinear Dynamics},
   year    = {2021},
   volume  = {103},
   pages   = {137-155},
   note    = {https://doi.org/10.1007/s11071-020-06109-0},
}
```

## License
**HarvesterOpt** is released under the MIT license. See the LICENSE file for details. All new contributions must be made under the MIT license.

<img src="logo/mit_license_red.png" width="10%"> 

## Institutional support

<img src="logo/logo_uerj_color.jpeg" width="10%">

## Funding

<img src="logo/faperj.jpg" width="20%"> &nbsp; &nbsp; <img src="logo/cnpq.png" width="20%"> &nbsp; &nbsp; <img src="logo/capes.png" width="10%">
