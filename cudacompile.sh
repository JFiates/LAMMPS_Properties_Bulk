#! /bin/bash
nvcc -O3 -prec-div false -prec-sqrt false -fmad false -Xcompiler -fopenmp -std=c++11 -gencode arch=compute_61,code=sm_61  -M -o "viscosityvf-cuda.d" "viscosityvf-cuda.cu"

nvcc -O3 -prec-div false -prec-sqrt false -fmad false -Xcompiler -fopenmp -std=c++11 --compile --relocatable-device-code=false -gencode arch=compute_61,code=compute_61 -gencode arch=compute_35,code=sm_35  -x cu -o  "viscosityvf-cuda.o" "viscosityvf-cuda.cu"

nvcc --cudart static --relocatable-device-code=false -gencode arch=compute_35,code=compute_35 -gencode arch=compute_61,code=sm_61 -link -o  "viscosity-cuda"  viscosityvf-cuda.o   -lgomp

