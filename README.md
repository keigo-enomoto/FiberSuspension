# FiberSuspension
This is the code used in my [graduate thesis](http://rheology.jp/nagoya/employees/keigo-enomoto/).

Watch the movie on [YouTube](https://www.youtube.com/watch?v=f9K-1ByvgVI)

## What is needed

- gcc (if you want use OpenMP)
- clang or gcc
- [PETSc](https://www.mcs.anl.gov/petsc/index.html)


## Prepare PETSc

1. Download zip file form HP
2. move zip file your favorite place (Here, put zip file in `lib` directory )
   ```
    mkdir lib
    mv path/to/file/petsc-3.xx.x.tar.gz lib
    cd lib
    tar xzf petsc-3.xx.x.tar.gz
   ```
3. configure (change the options according to your environment. refer to [Instllation page](https://www.mcs.anl.gov/petsc/documentation/installation.html))

   **example**
   - don't use OpenMP
   ```
    export PETSC_DIR=/path/to/lib/petsc-3.xx.x

    ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran -- download-fblaslapack --with-mpi=0 PETSC_DIR=$PETSC_DIR
   ```
   - use OpenMP
    ```
    export PETSC_DIR=/path/to/lib/petsc-3.xx.x

    ./configure --with-cc=gcc --with-cxx=g++ --with-fc=gfortran -- download-fblaslapack --with-openmp --with-mpi=0 PETSC_DIR=$PETSC_DIR
    ```

    After configure, do `make all`, `make check` as the log said.
4. set Environmental variables `PETSC_DIR` and `PETSC_ARCH`
   set `PETSC_DIR` and `PETSC_ARCH` in `~/.bashrc` or `~/.zshrc`




## make initial coordinate

Use `src/mk_fiber` and `fp_fiber3/init.c`

1. make
    ```
    cd src/mk_fiber
    make 2d
    cd ../fp_fiber3
    make init
    cd ../..
    ```

2. set parameter
   ```
   cd output/mk_fiber
   ```
    and edit `mk_particle.sh`

    ```sh
    # change here ==========
    xx_arr=( 1 )
    yy_arr=( 60 ) # the number of particle in y axis
    ap_arr=( 10 ) # aspect ratio of fiber
    nf_arr=( 42 ) # number of fiber
    ss_arr=( 98 ) # seed for initial coordinate

    # ========================
    ```

3. ```
   ./mk_particle.sh
    ```

## calculate

1. make (don't use OpenMP)
   ```
    cd src/fp_fiber3
    make
    cd ..
   ```
   make (use OpenMP)
   ```
    cd src/fp_fiber3
    ./omp.sh
    cd ..
   ```
2. edit parameter
   ```
   cd output
   cp -r example test1
   cd test1
   ```
   and edit `input_cluster.txt`
   ```txt
    DIM 2           # dimension(Now, only 2 is available)
    Rho_s 1.0       # mass density of fiber
    Nu_f 10000      # viscosity of fluid
    Nu_s 10000      # viscosity of fiber
    strain_rate 10
    Gravity 0.0
    DT 0.00001
    FINISH_TIME 0.2
    OUTPUT_INTERVAL 20
    Init_place ./../mk_fiber/18f20-10_3600_60-60_2d.dat # the file for initial coordinate 
    output_dir .    # where to output
    ks 1.0E+09      # spring constant
    kb 5.0E+08
    VTU_flag 0     # 1 -> output vtu file, 0 -> don't output
   ```
3. ```
   ./execute.sh
   ```
