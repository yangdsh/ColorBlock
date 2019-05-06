![alt text](image/result1.bmp)
# ColorBlock: Parallel Image Style Transformation
This is a course project of CMU 15-618: Parallel Computer Architecture and Programming. 

It uses CUDA accelerated K-means, OpenMP accelerated BFS, and CUDA accelerated edge smoothing to transform the image into color-block/glass-painting style. 
Its image processing engine is written in C++ without using any image processing library.
Sequential codes are also provided in order to make validations and performance comparisons.

## Authors
Dongsheng Yang, Xuyang Fang

## Requirements
```sh
pillow
MPI
CUDA
```

## Running
```sh
make
sh oilPainting.sh
```
You can change the input image and the parameters in the script.
