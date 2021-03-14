
# FP_FIBER3

The difference from fp_fiber2
- Use OpenMP directive
- The timing of updating wall particle's position (labnote vol.2 p193)

[object to calculate]
- poiseuille flow
- couette flow (both side of wall will move)
- with fiber

[feature]
- 2D and 3D
- cell list
- NonDimension
- viscous term --> explicit(mps : main.c), implicit(mps_i : main_i.c)
- pressure     --> implicit(with PETSc and single node)
    |            | poiseuille | couette |
    |:----------:|:----------:|:-------:|
    |  Gravity   |    > 0.0   |  = 0.0  |
    | wall_speed |    = 0.0   |  > 0.0  |

- fiber particle is treated as a fluid particle with bond constraint
- Acceleration[] of wall particle is \rho * Du/Dt --> calculate stress


[input]
1. input.txt
2. Initial arrangement
3. output directory

4. ----------------------------------------
```
DIM 2
Rho_s 1.0
Nu_f 1.0
Nu_s 1.0
strain_rate 0.01
Gravity 0.0
DT 0.1
FINISH_TIME 100.0
OUTPUT_INTERVAL 10
output_dir .
Init_place ./../../couette_1760.dat
ks 1.0
kd 1.0
VTU_flag 1 or 0
```

ParticleDistance = 1.0
Rho_f = 1.0 ;                                          # (s --> solid, f --> fluid)
t = t*/(ParticleDistance * ParticleDistance / Nu_f)
  = Nu_f * t* (t :dimensionless, t* :dimension)        # Nu -> Kinematic Viscosity
g = 1.0

2. -----------------------------------------------


3. -----------------------------------------------
If the path to output directory is written uncorrectoly in input file, it will cause error.

[output]
1. input_parameters.txt
2. output_*.dat : record vx and y
3. output_analysis.txt

[how to execute]

```command line
  path/to/mps input.txt
```
[change note]
