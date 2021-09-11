# LexLS
A fast C++ solver for lexicographic least-squares problems based on the *lexicographic QR* (l-QR) presented in the paper [Efficient resolution of potentially conflicting
linear constraints in robotics](https://hal.inria.fr/hal-01183003/document), by D. Dimitrov, A. Sherikov and P.-B. Wieber.

## Lexicographic Least-Squares
Let's consider <img src="https://render.githubusercontent.com/render/math?math=\color{gray}p"> sets of constraints <img src="https://render.githubusercontent.com/render/math?math=\color{gray} b_k \leq C_k x \leq u_k,\ k=1..p">, potentially conflicting. We consider that each of these sets have a different priority level, where `1` has the higher priority and `p` the lower. We are looking for a solution that first minimize the violation of <img src="https://render.githubusercontent.com/render/math?math=\color{gray} b_1 \leq C_1 x \leq u_1">, in the least-square sense (that is minimizing <img src="https://render.githubusercontent.com/render/math?math=\color{gray}\left\|v_1\right\|^2"> such that <img src="https://render.githubusercontent.com/render/math?math=\color{gray} b_1 \leq C_1 x - v_1\leq u_1"> is verified), then minimizing the violation of <img src="https://render.githubusercontent.com/render/math?math=\color{gray} b_2 \leq C_2 x \leq u_2"> without increasing <img src="https://render.githubusercontent.com/render/math?math=\color{gray}\left\|v_1\right\|^2"> and so on.

Formally this can be written as <br />
<img src="https://render.githubusercontent.com/render/math?math=\color{gray}\begin{align*}\mathrm{lex.\ min.}_{x,v}\ \ \left(\left\|v_1\right\|^2, \ldots, \left\|v_p\right\|^2\right)\\ \mathrm{s.t.} \ \ \ \ \  b_k \leq C_k x - v_k \leq u_k,\ \ k=1..p \end{align*}">

`LexLS` solves the above problem with a primal, active-set approach, where all priority levels are tackled together.

## Scope
`LexLS` was optimized for problems with few changes of active set, leveraging warm-start to limit the number of iterations between the resolution of closely related problems as can be found in e.g succesive inverse kinematics problems. The l-QR decomposition is highly optimized to run from scratch, with speed comprised between Eigen's LU and QR, depending on the number of priority levels. No update mechanism has been implemented yet, so that each iteration of the solver performs a full decomposition.

## Build and use with CMake

### Dependencies

- [Eigen](https://eigen.tuxfamily.org/)

For documentation generation:
- [Doxygen](https://www.doxygen.nl/index.html) (`sudo apt install doxygen` on Debian-based systems)
- pdflatex with extra packages and extra fonts (`sudo apt install texlive-base texlive-science texlive-fonts-extra` on Debian-based systems)

### Building

This package can be built with CMake:
```bash
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=RelWithDebInfo -GNinja
ninja && sudo ninja install
```

You can customize the build with the following options:
- `BUILD_TESTING` build unit tests (`ON` by default)
- `INSTALL_HTML_DOCUMENTATION` build and install the HTML documentation (`ON` by default, automatically skipped if doxygen is missing)
- `INSTALL_PDF_DOCUMENTATION` build and install the PDF documentation (`ON` by default, automatically skipped if pdflatex is missing, does not check for missing texlive packages)

### Using it in your project

Once the project is installed, you can use the following to use lexls in your project:

```cmake
find_package(lexls REQUIRED)
target_link_libraries(MyProject PUBLIC lexls::lexls)
```

## Credits
The solver was mainly implemented by Dimitar Dimitrov, with additionnal code from Alexander Sherikov and Nestor Bohorquez, with supervision from Pierre-Brice Wieber.
