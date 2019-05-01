from PIL import Image
	
im = Image.open("../image/result.bmp")
im.save("../image/download.jpg")
im.show()
