# ColorBlock: Parallel Image Style Transformation
This is a course project of CMU 15-618: Parallel Computer Architecture and Programming. 
It uses CUDA accelerated K-means and OpenMP accelerated BFS to transform the image into color-block style. 
Its image processing engine is written in C++ without using any image processing library.
Sequential codes and MPI accelerated codes are also provided in order to make validations and performance comparisons.

## Authors
Dongsheng Yang, Xuyang Fang

## Requirements
```sh
pillow
MPI
CUDA
```

## Running
First, move the image you want to process to ``image/upload.jpg``
```sh
make
sh run.sh
```
