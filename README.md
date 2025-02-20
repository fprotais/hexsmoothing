# hexsmoothing
Hexahedral and mixed elements smoothing.

# Prerequisite

If you want to run the parallel version of the code, you'll need Eigen (>3.1) and OpenMP (> 4.0).
```sh
brew install eigen
```
```sh
sudo apt install libeigen3-dev
```
```sh
sudo dnf install eigen3-devel
```
# building

You can build the library and its corresponding executable doing: 
```
git clone https://github.com/fprotais/hexsmoothing
cd hexsmoothing
mkdir build
cd build
cmake -DPARALLEL_MIXED_SMOOTHING=ON .. 
make -j
```
There might be some linking issues, redoing `cmake .. && make -j` usually solve those. The lib was made to be easy to link against. It is poorly documented though. Don't hesitate to raise an issue if you have any questions. 

# Mixed elements smoothing

You can run a mixed elements smoother/untangler (Tets, Pyramids, Wedges and Hexes) using the following executable, from the build directory: 
```
./parallelMixedSmoothing ../mixedElementsMesh.vtk result.vtk 
```
An optional parameter is the number of iteration of smoothing. Input and output files can be .vtk or .mesh, for me to be able to read mixed meshes. Boundary of the model will be fixed. For more precise control, I recommend looking into the C++ code and directly using the library.  

The code will first do an untangling loop, then try to optimize the quality further. I am still tweaking the parameters so they may not be perfect, particularly on large meshes. 
# Cite the repo
You can directly cite this repository: 
```
@misc{Hexsmoothing2023,
  author = {Francois, Protais},
  title = {Ph.D. smoothing codes},
  year = {2023},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/fprotais/hexsmoothing}},
}
```
The best explanation about how the code works is available in my Ph.D. thesis (_in French_):  
```
@phdthesis{protaisPHD2023,
  TITLE = {{Maillage {\`a} dominante Polycube}},
  AUTHOR = {Protais, Fran{\c c}ois},
  URL = {https://hal.science/tel-03775686},
  NUMBER = {2022LORR0144},
  SCHOOL = {{Universit{\'e} de Lorraine}},
  YEAR = {2022},
  MONTH = Oct,
  KEYWORDS = {Meshing ; hexahedra ; deformation ; Polycube ; Maillage ; hexa{\`e}dres ; d{\'e}formation ; Polycube},
  TYPE = {Theses},
  PDF = {https://hal.science/tel-03775686/file/finalVersion.pdf},
  HAL_ID = {tel-03775686},
  HAL_VERSION = {v1},
}
```
