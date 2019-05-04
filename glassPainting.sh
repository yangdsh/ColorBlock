wget "https://cdn3.tropicalsky.co.uk/images/800x600/lake-kawaguchiko-autumn.jpg" -O image/upload.jpg
python script/jpg2bmp.py
./release/glassPainting image/input.bmp image/result.bmp 3 3
python script/bmp2jpg.py
