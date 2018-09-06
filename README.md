# GPisMap -- Gaussian Process Implicit Surface Map

The major contributors of this repository include [Bhoram Lee](https://github.com/leebhoram) and [HUANG Zonghao](https://github.com/huangzonghao).

## Introduction

<+Introduction to the paper and implementation+>
  
## License

Licensed under [GNU General Public License version 3](https://www.gnu.org/licenses/gpl-3.0.html).

## Main Results

<+Results and video link+>

## Requirements: Software

1. [Eigen](http://eigen.tuxfamily.org/)

2. [MATLAB](https://www.mathworks.com/products/matlab.html)

## Compiling and Running

1. Clone this repository
```
git clone https://github.com/leebhoram/GPisMap.git
```

2. Cd to the mex directive in MATLAB
```
cd mex
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

## Misc.

Code has been tested under:

- Ubuntu 16.04 with Intel Core i7-4900MQ @ 2.90GHz

## Note

This work has been submitted to [ICRA 2019](https://www.icra2019.org/) for
review.

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
