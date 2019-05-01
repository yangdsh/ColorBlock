#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>
#include <math.h>
#include <algorithm>
using namespace std;

typedef unsigned int uint32_t;

#define INF (100000000.0)
#define THRESHOLD	(1e-7)
#define dim 3
#define PI  3.14159265f

int ReadBMP(string strFile, int &size,int & width, int & height , double *&pixels, char *&head) {
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
	head = (char *)malloc(54* sizeof(char));
	fread(head,54,sizeof(char),fin);
	fseek(fin,54,SEEK_SET);
	//read the pixels
	pixels = (double *)malloc(size * dim * sizeof(double));
	for (int i = 0; i < size; i ++) {
		for (int j = 0; j < dim; ++j) {
			unsigned char color;
			fread(&color, 1, sizeof(char), fin);
			pixels[i*dim + j] = double(color);
		}
	}
	fclose(fin);
}
	
int WriteBMP(string strFile, int size, double *pixels, char *head) {	
	FILE *fout ;

	if ((fout=fopen(strFile.c_str(),"wb"))==NULL) {
		cout<<"create the bmp file error!"<<endl;
		return 0;
	}
	fwrite(head,54,sizeof(char),fout);

	for (int i = 0; i < size; ++i) {
	    for (int j = 0 ; j < dim; j ++)	{
	    	char temp = (char) pixels[i*dim + j];
			fwrite(&temp, 1, sizeof(char), fout);
		}
	}
	fclose(fout);

	free(pixels);
}

void conv(int size, double* weight, int WIDTH,
int HEIGHT,
double* input_nopad,
double* output){
	double * input = (double *)calloc((WIDTH+size-1) * (HEIGHT+size-1), sizeof(double));
	int pad_size=int(size/2);
	for (int j=0; j<HEIGHT; j++) {
	 	for (int i=0; i<WIDTH; i++) 
	 		input[(j+pad_size)*(WIDTH+size-1)+i+pad_size]=input_nopad[j*WIDTH+i];
 	}
/*	for (int j=0; j<pad_size; j++) {
	 	for (int i=0; i<pad_size; i++) {
	 		input[j*(WIDTH+size-1)+i]=input_nopad[0];
	 		input[j*(WIDTH+size-1)+i]=input_nopad[0];
	 		input[j*(WIDTH+size-1)+i]=input_nopad[0];
	 		input[j*(WIDTH+size-1)+i]=input_nopad[0];
	 	}
 	}*/	

	for (int j=0; j<HEIGHT; j++) {
	 	for (int i=0; i<WIDTH; i++) {
			double tmp = 0;
	 		for (int jj=0; jj<size; jj++)
	 			for (int ii=0; ii<size; ii++)
	 				tmp += input[(j+jj)*(WIDTH+size-1) + (i+ii)] * weight[jj*size + ii];
	 	output[j*WIDTH + i] = int(tmp);
	 	}
	}
}

void gray(int size, double* input, double* output){
	for (int i = 0; i < size; ++i)
	{
		output[i]=int((input[3*i]+input[3*i+1]+input[3*i+2])/3);
	}
}

void colorful(int size, double* input, double* output){
	for (int i = 0; i < size; ++i)
	{
		output[3*i]=int(input[i]);
		output[3*i+1]=int(input[i]);	
		output[3*i+2]=int(input[i]);	
	}
}



void sobel_suppression(int width,int height,double* sobel_ximg,
	double* sobel_yimg,double* max_suppressed){
	int size=width*height;
	double max_sobel=0;
	double angle=0;
	double *sobel_img=(double *)malloc(size  * sizeof(double));	
	double *sobel_direction=(double *)malloc(size  * sizeof(double));
	for (int i = 0; i < size; ++i){
		sobel_img[i]=sqrt(sobel_ximg[i]*sobel_ximg[i]+sobel_yimg[i]*sobel_yimg[i]);
		max_sobel=sobel_img[i]>max_sobel ? sobel_img[i] : max_sobel;

		if ((sobel_ximg[i] != 0.0) || (sobel_yimg[i] != 0.0)) {
			angle = atan2(sobel_yimg[i], sobel_ximg[i]) * 180.0 / PI;
		} else {
			angle = 0.0;
		}
		if(angle<0)
			angle+=180;
		sobel_direction[i] =angle;
	}
/*	for (int i = 0; i < size; ++i){
		sobel_img[i]=int(255.0f * sobel_img[i] / max_sobel);
	}*/

	int maxmum=0;
	for (int i = 10; i < width-10; ++i)
	{
		for (int j = 10; j < height-10; ++j)
		{
			if(sobel_img[j*width+i]>maxmum){
				maxmum=sobel_img[j*width+i];
			}
		}
	}
	for (int i = 0; i < size; ++i)
	{
		sobel_img[i]=int(255.0f * sobel_img[i] / maxmum);
		if(sobel_img[i]>=255){
			sobel_img[i]=255;
		}
	}

	double p=0;
	double q=0;

	for (int i = 1; i < width-1; ++i){
		for (int j = 1; j < height-1; ++j){
			double angle=sobel_direction[j*width+i];
			if (angle < 22.5 || angle >= 157.5) {
				p=sobel_img[j*width+i-1];
				q=sobel_img[j*width+i+1];
			} else if (angle >= 22.5 && angle < 67.5) {
				p=sobel_img[(j+1)*width+i-1];
				q=sobel_img[(j-1)*width+i+1];
			} else if (angle >= 67.5 && angle < 112.5) {
				p=sobel_img[(j+1)*width+i];
				q=sobel_img[(j-1)*width+i];
			} else{
				p=sobel_img[(j+1)*width+i+1];
				q=sobel_img[(j-1)*width+i-1];
			}
			if(sobel_img[j*width+i]>=p && sobel_img[j*width+i]>=q )
				max_suppressed[j*width+i]=sobel_img[j*width+i];
		}
	}
}

void hysteresis(int width,int height,double* max_suppressed,double* hysteresis_output,int thre1,int thre2){
	int size=width*height;
	int sum=0;
	for(int i=0;i<size;i++){
		if(max_suppressed[i]>=thre2){
			hysteresis_output[i]=255;
			sum++;
		}else if (max_suppressed[i]<thre1){
			hysteresis_output[i]=0;			
		}else{
			hysteresis_output[i]=thre1;
		}
	}
	//printf("%d\n",sum );

	for (int i = 1; i < width-1; i++) {
		for (int j = 1; j < height - 1; j++) {
			int src_pos = i + (j * width);
			if (hysteresis_output[src_pos] == thre1) {
				if (hysteresis_output[src_pos - 1] == 255 || 
					hysteresis_output[src_pos + 1] == 255 ||
					hysteresis_output[src_pos - 1 - width] == 255 ||
					hysteresis_output[src_pos + 1 - width] == 255 || 
					hysteresis_output[src_pos + width] == 255 ||
					hysteresis_output[src_pos - width] == 255 ||
					hysteresis_output[src_pos + width - 1] == 255 ||
					hysteresis_output[src_pos + width + 1] == 255) {
						hysteresis_output[src_pos] = 255;
				} else {
					hysteresis_output[src_pos] = 0;
				}
			} 
		}
	}
}
void canny(double* pixels,double *canny_output, int size, int width,int height, int thre1,int thre2 ){
 	double Gaussian[]={2,4,5,4,2,4,9,12,9,4,5,12,15,12,5,4,9,12,9,4,2,4,5,4,2};
	for (int i=0;i<25;i++){
		Gaussian[i]=Gaussian[i]/159;
	}
	double sobel_x[]={-1,0,1,-2,0,2,-1,0,1};	
	double sobel_y[]={1,2,1,0,0,0,-1,-2,-1};	

	double *grayimg=(double *)malloc(size* sizeof(double));	
	gray(size,pixels,grayimg);

	double *gaussianimg=(double *)malloc(size* sizeof(double));	
	conv(5, Gaussian, width,height, grayimg, gaussianimg);
	double *sobel_ximg=(double *)malloc(size  * sizeof(double));
	double *sobel_yimg=(double *)malloc(size  * sizeof(double));
	conv(3, sobel_x, width,height, gaussianimg, sobel_ximg);
	conv(3, sobel_y, width,height, gaussianimg, sobel_yimg);
	double *max_suppressed=(double *) calloc(size,sizeof(double));	
	sobel_suppression(width,height,sobel_ximg,sobel_yimg, max_suppressed);

/*	double *hysteresis_output=(double *) malloc(size*sizeof(double));*/
	hysteresis(width,height,max_suppressed,canny_output,thre1,thre2);
}

void dilate(double* input, double* output, int kernel_size,int width,int height){
	double* weight=(double *)malloc(kernel_size *kernel_size * sizeof(double));	
	for(int i=0;i<kernel_size*kernel_size;i++)
		weight[i]=1;
	conv(kernel_size, weight, width,height,input,output);
	for(int i=0;i<width*height;i++){
		if(output[i]>0)
			output[i]=255;
	}
}

void blur(double* input, double* output, int kernel_size,int width,int height){
	double* weight=(double *)malloc(kernel_size *kernel_size * sizeof(double));	
	for(int i=0;i<kernel_size*kernel_size;i++)
		weight[i]=float(1)/float(kernel_size*kernel_size);
	conv(kernel_size, weight, width,height,input,output);
	for(int i=0;i<width*height;i++){
		if(output[i]>127){
			output[i]=255;
		}else{
			output[i]=0;			
		}
	}

	int sum=0;
	for (int i = 0; i < width*height; ++i)
	{
		if(output[i]==255)
			sum++;
	}
	//printf("%d\n", sum );
}

void edge2color(double* origin, double* input, double* output,int width,int height){
	int size=width*height;
	colorful(size,input, output);
	double* tmp=(double *)malloc(size*dim * sizeof(double));	
	for(int i=0;i<size*dim;i++){
		output[i]+=origin[i];
		if(output[i]>255)
			output[i]=255;
	}
	while(1){
		for (int i = 0; i < size*dim; ++i)
		{
			tmp[i]=output[i];
		}
		for(int i=0;i<width;i++){
			for (int j = 0; j < height; j++)
			{
				if(tmp[3*(j*width+i)]!=255 || tmp[3*(j*width+i)+1]!=255 || tmp[3*(j*width+i)+2]!=255)
					continue;
				int color1=0;
				int color2=0;
				int color3=0;
				int sum_color=255*3;				
				for(int p=max(0,i-1);p<min(width,i+2);p++){
					for(int q=max(0,j-1);q<min(height,j+2);q++){
						int tmp_sum=tmp[3*(q*width+p)]+tmp[3*(q*width+p)+1]+tmp[3*(q*width+p)+2];
						if (tmp_sum<sum_color){
							sum_color=tmp_sum;
							output[3*(j*width+i)]=tmp[3*(q*width+p)];
							output[3*(j*width+i)+1]=tmp[3*(q*width+p)+1];
							output[3*(j*width+i)+2]=tmp[3*(q*width+p)+2];
						}
					}
				}
			}
		}
		int flag=false;
		for(int i=0;i<size;i++){
			if(output[3*i]==255 && output[3*i+1]==255 && output[3*i+2]==255){
				flag=true;
				break;
			}
		}
		if(!flag)
			break;
	}
}


int main(int argc, char ** argv){
	double * pixels;
	char *head;
	string inputFile = "colorBlock.bmp";
	string outputFile = "finalResult.bmp";
	int size = 0;
	int width=0;
	int height=0;
  int dilate_size=5;
  int blur_size=5;
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
		
		cout << "dilate size: ";
		cin.getline(line, 100);
		if (strlen(line) > 0)
			dilate_size = atoi(line);
		else
      cout << dilate_size << endl;
    
    cout << "blur size: ";
		cin.getline(line, 100);
		if (strlen(line) > 0)
			blur_size = atoi(line);
		else
      cout << blur_size << endl;
	}
	ReadBMP(inputFile, size,width,height, pixels, head);
  clock_t cstart=clock();
	double *canny_output=(double *)malloc(dim*size * sizeof(double));		
	canny(pixels,canny_output,size,width,height,5,6);

	double *dilate_output=(double *)malloc(size * sizeof(double));	
	dilate(canny_output,dilate_output,dilate_size,width,height);

	double* final_output=(double *)malloc(dim*size * sizeof(double));		

	double *blur_output=(double *)malloc(size * sizeof(double));	
	blur(dilate_output,blur_output,blur_size,width,height);
	
	clock_t cstart_e2c=clock();
	double *output=(double *)malloc(dim*size * sizeof(double));	
	edge2color(pixels,blur_output,output,width,height);
	clock_t cend_e2c=clock();
	clock_t cend=clock();

/*	colorful(size,blur_output, final_output);*/
	WriteBMP(outputFile, size, output, head);

	cout<<"total smoothing time: "<<(cend-cstart)/1000000.0 << " s" <<endl;
	cout<<"edge2color time: "<<(cend_e2c-cstart_e2c)/1000000.0 << " s" <<endl;
	return 0;
}
