python script/jpg2bmp.py
./release/kmeans_cuda image/input.bmp image/kmeans.bmp 100
./release/colorBlock_omp image/kmeans.bmp image/result.bmp 15 3
python script/bmp2jpg.py
