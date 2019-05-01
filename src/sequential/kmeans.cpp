#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>
#include "../util/cycletimer.h"

using namespace std;

typedef unsigned int uint32_t;

#define INF (1e20)
#define THRESHOLD	(1e-2)
#define DIM 3

#define INIT1
#define DEBUG

inline double sqr(double d) {
	return d * d;
}

inline double abs(double d) {
	if (d > 0) return d; else return -d;
}

inline double CalcDistance(const double *p, const double *q) {
	uint32_t i;
	double r = 0;
	for (i = 0; i < DIM; ++i)
		r += sqr(p[i] - q[i]);
	r = sqrt(r);
	return r;
}

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
	
int WriteBMP(string strFile, int size, double *pixels, char *&head) {	
	FILE *fout ;
	fout=fopen(strFile.c_str(),"wb");
	if (fout==NULL) {
		cout<<"create the bmp file error!"<<endl;
		return 0;
	}
	fwrite(head, sizeof(char), 54, fout);

	for (int i = 0; i < size; i++) {
	    for (int j = 0; j < DIM ; j ++)	{
			unsigned char temp = (unsigned char) pixels[i*DIM+j];
			fwrite(&temp, sizeof(char), 1, fout);
		}
	}
	fclose(fout);
  return 1;
}

void KMeans(uint32_t n, uint32_t k, double * pixels) {
	uint32_t i; 
	uint32_t j;
	uint32_t c = 0;
	double *center;
	uint32_t *clst;
	uint32_t *clst_size;
	double min_cr = 50, max_cr = 200;
	double old_sum_dis, new_sum_dis;
	double min_dis, dis;

	int cstart=0;
	int cend=0;
	int assign_time=0;
	int update_time=0;

	center = (double *)malloc(k * DIM * sizeof(double));
  srand(0);
  #ifdef INIT0
	for (j = 0; j < DIM; ++j) {
		min_cr = INF; max_cr = -INF;
		for (i = 0; i < n; ++i) {
			dis = pixels[i * DIM + j];
			if (dis < min_cr) min_cr = dis;
			if (dis > max_cr) max_cr = dis;
		}
		for (i = 0; i < k; ++i) {
			center[i * DIM + j] = min_cr + (double)(rand() * (max_cr - min_cr)) / (double)RAND_MAX;
		}
	}
  #endif
  #ifdef INIT1
  for (i = 0; i < k; ++i)
    for (j = 0; j < DIM; ++j) {
			center[i * DIM + j] = min_cr + (double)(rand() * (max_cr - min_cr)) / (double)RAND_MAX;
		}
  #endif


	clst_size = (uint32_t *)malloc(k * sizeof(uint32_t));

	clst = (uint32_t *)malloc(n * sizeof(uint32_t));
  
  #ifdef INIT2
  c = 0;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < DIM; ++j)
			center[c * DIM + j] += pixels[i];
		clst_size[c] ++;
		c = (c+1) % k;
	}
	for (i = 0; i < k; ++i) {
	  for (j = 0; j < DIM; ++j)
			center[i * DIM + j] /= clst_size[i];
	}
  #endif
	
	old_sum_dis = INF;
  int iter = 0;
	while (1) {
		iter ++;
    new_sum_dis = 0;
    cstart = currentSeconds(); 
		for (i = 0; i < n; ++i) {
			min_dis = INF;
			for (j = 0; j < k; ++j) {
				dis = CalcDistance(&pixels[i * DIM], &center[j * DIM]);
				if (dis < min_dis) {
					min_dis = dis;
					c = j;
				}
			}
			clst[i] = c;
			new_sum_dis += min_dis;
		}
  cend = currentSeconds();
  assign_time+=cend-cstart;
    #ifdef DEBUG
    //cout << "k-means iteration: " << iter  << ", total distance: " << new_sum_dis << endl;
    #endif
		if (abs(old_sum_dis - new_sum_dis) < THRESHOLD || old_sum_dis - new_sum_dis < 0) break;

		for (i = 0; i < k; ++i) {
			for (j = 0; j < DIM; ++j)
				center[i * DIM + j] = 0;
			clst_size[i] = 0;
		}

	  cstart = currentSeconds();	
		for (i = 0; i < n; ++i) {
			c = clst[i];
			for (j = 0; j < DIM; ++j)
				center[c * DIM + j] += pixels[i * DIM + j];
			++ clst_size[c];
		}
		for (i = 0; i < k; ++i) {
			for (j = 0; j < DIM; ++j)
				center[i * DIM + j] /= clst_size[i];
		}
  cend = currentSeconds();
  update_time+=cend-cstart;
		old_sum_dis = new_sum_dis;
	}

	for (i = 0; i < n; ++i) {
		c = clst[i];
		for (j = 0; j < DIM; ++j)
			pixels[i * DIM + j] = center[c * DIM + j];
	}

	free(clst_size);
	free(center);
	free(clst);
	cout<<assign_time<<" "<<update_time<<endl;
}

int main(int argc, char ** argv) {
	double * pixels;
	char * head;	
	string inputFile = "input.bmp";
	string outputFile = "kmeans.bmp";
	int size = 0;
	int k = 20;	
  double cstart, cend;
  
  cout << "-----------------K-means part------------------" << endl;
	
	if (argc == 1) {
		cout << "input file name: ";
		char line[100];
		cin.getline(line, 100);
		if (strlen(line) > 0)
			inputFile = string(line);
		else
      cout << inputFile << endl;
		
		cout << "output file name: ";
		cin.getline(line, 100);
		if (strlen(line) > 0)
			outputFile = string(line);
		else
      cout << outputFile << endl;
		
		cout << "number of colors for k-means: ";
		cin.getline(line, 100);
		if (strlen(line) > 0)
			k = atoi(line);
		else
      cout << k << endl;
	}
  cstart = currentSeconds(); 
	ReadBMP(inputFile, size, pixels, head);
  cend = currentSeconds();
  cout << "read image time: " << cend-cstart << " s" <<endl;  
  cstart = currentSeconds(); 
	KMeans(size, k, pixels);
  cend = currentSeconds();
  cout << "K-means time: " << cend-cstart << " s" <<endl;  
  cstart = currentSeconds(); 
  WriteBMP(outputFile, size, pixels, head);
  cend = currentSeconds();
  cout << "write image time: " << cend-cstart << " s" <<endl;  
	free(pixels);
	free(head);
	cout << "-----------------K-means done------------------" << endl;
}
