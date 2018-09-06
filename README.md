# GPisMap 

This repository contains source code and demo files for our paper **Online
Continuous Mapping using Gaussian Process Implicit Surfaces (GPIS)**, which is
currently submitted to [IEEE ICRA 2019](https://www.icra2019.org/).

The representation of the environment strongly affects how robots can move and
interact with it. The paper presents an online approach for continuous mapping
using Gaussian Process Implicit Surfaces (GPIS). Compared with grid-based
methods, GPIS better utilizes sparse measurements to represent the world
seamlessly. It provides direct access to the signed-distance function (SDF) and
its derivatives which are invaluable for other robotic tasks and incorporates
uncertainty in the sensor measurements. Our approach incrementally and
efficiently updates GPIS by employing a regressor on observations and a spatial
tree structure.
 
## License

Licensed under [GNU General Public License version 3](https://www.gnu.org/licenses/gpl-3.0.html).

## Requirements: Software

1. [Eigen3](http://eigen.tuxfamily.org/)

2. [MATLAB](https://www.mathworks.com/products/matlab.html)

## Compiling and Running

1. Clone this repository
```
git clone https://github.com/leebhoram/GPisMap.git
```

2. Cd to the mex directive in MATLAB
```
cd GPisMap/mex
```

3. Compile the mex functions by executing the make script.
    * Setup mex 
    ```
    mex -setup
    ```
    * Run the make scripts
    ```
    make_GPisMap
    make_GPisMap3
    ```

4. Run the demo scripts

    * For 2D 
    ```
    run('../matlab/demo_gpisMap.m')
    ```
    * For 3D 
    ```
    run('../matlab/demo_gpisMap3.m')
    ```

5. Trouble shooting
    * If mex complains about not finding eigen, configure the eigen path
        appropriately in both `make_GPisMap.m` and `make_GPisMap3.m`

## Video  
[![](http://img.youtube.com/vi/_EqeoLeHzXU/0.jpg)](http://www.youtube.com/watch?v=_EqeoLeHzXU "Online Continuous Mapping using GPIS")

## Contributors

The major contributors of this work include
[Bhoram Lee](https://github.com/leebhoram),
[Clark Zhang](https://github.com/chickensouple) and
[HUANG Zonghao](https://github.com/huangzonghao).

## Misc.

Code has been tested under:

- Ubuntu 16.04 with MATLAB R2018a on Intel Core i7-4900MQ @ 2.90GHz

<!-- ## Citation
   - 
   - If you find GPisMap useful in your research, please consider citing:
   - ```
   - @article{<++>,
   -     Author = {Bhoram Lee, Clark Zhang, Zonghao Huang, Daniel D. Lee},
   -     Title = {Oneline Continuous Mapping using Gaussian Process Implicit Surfaces},
   -     Journal = {<++>},
   -     Year = {<++>}
   - }
   - 
   - ```
   -->
