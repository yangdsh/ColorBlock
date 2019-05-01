#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>

using namespace std;
#define DIM 3

int ReadBMP(string strFile, int &size, double *& pixels, char *& head) {
	
	int width, height ;
	FILE *fin ;
	fin=fopen(strFile.c_str(),"rb");
	
	//check file pointer
	if(fin == NULL) {
		cout<<"file open error!"<<endl;
		return 0;
	}
	//check file type
	short bfType;
	fread(&bfType,1,sizeof(short),fin);
	if(0x4d42!=bfType) {
		cout<<"the file is not a bmp file!"<<endl;
		return 0;
	}
	//get the number of pixels 
	fseek(fin,18,SEEK_SET) ;
	fread(&width,1,sizeof(int),fin);
	fread(&height,1,sizeof(int),fin);
	size = width * height ;
	//check the color map
	fseek(fin,28,SEEK_SET) ;
	unsigned short colors ;
	fread(&colors,1,sizeof(unsigned short),fin);
	if (colors != 24 ) {
		cout << "The color map must be 24 bits" << endl ;
		return 0 ;
	}
	//get the file header
	fseek(fin,0,SEEK_SET);
	head = (char *)malloc(54 * sizeof(char));
	fread(head,54,sizeof(char),fin);
		
	//read the pixels
	fseek(fin,54,SEEK_SET);
	pixels = (double *)malloc(size * DIM * sizeof(double));
	for (int i = 0; i < size; i ++) {
		for (int j = 0; j < DIM; ++j) {
			unsigned char color;
			fread(&color, 1, sizeof(char), fin);
			pixels[i*DIM + j] = double(color);
		}
	}
	fclose(fin);
  return 1;
}

int compare_image(double *image1, double *image2, int size) {
  int cnt = 0;
  for (int i = 0 ; i < size; i ++)
    cnt += image1[i] == image2[i];
  return cnt;
}
	
int main(int argc, char ** argv) {
	double *image1, *image2;
	char * head;	
	string inputFile1 = "kmeans.bmp";
	string inputFile2 = "kmeans_p.bmp";
	int size = 0;
  
  cout << "-----------------compare image------------------" << endl;
	
	if (argc == 1) {
		cout << "input file1 name: ";
		char line[100];
		cin.getline(line, 100);
		if (strlen(line) > 0)
			inputFile1 = string(line);
		else
      cout << inputFile1 << endl;
		
		cout << "input file2 name: ";
		cin.getline(line, 100);
		if (strlen(line) > 0)
			inputFile2 = string(line);
		else
      cout << inputFile2 << endl;
	}
  ReadBMP(inputFile1, size, image1, head);
  ReadBMP(inputFile2, size, image2, head);
 
  double cnt = compare_image(image1, image2, size);
  free(image1);
  free(image2);
  free(head);
  cout << "similarity: " << cnt/size * 100 << " %" << endl;
  
	cout << "-----------------compare done------------------" << endl;
}
