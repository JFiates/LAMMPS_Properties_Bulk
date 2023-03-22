# LAMMPS_Properties_Bulk
This repository contains codes written in the C programming language for calculating properties based on bulk simulation with output from LAMMPS. Additionally, a parallelized code for calculating viscosity is also available, which can run on both CPUs and GPUs.

To run the code using CUDA, it is necessary to compile it using nvcc, and the GPU toolkit should be installed beforehand.

To run the serial codes in C, use the following command:
gcc name_of_code.c -lm -O3
To execute:
echo input_file | ./a.out

To run the CPU/GPU code, use the following command:
bash cudacompile.sh input_file
