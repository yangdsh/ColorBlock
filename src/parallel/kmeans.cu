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
#define BLOCK_SIZE 512
#define c_num 50

__constant__ unsigned int d_size;
__constant__ unsigned int d_k;

inline double CalcDistance(const double *p, const double *q) {
	uint32_t i;
	double r = 0;
	for (i = 0; i < DIM; ++i)
		r += sqrt(p[i] - q[i]);
	r = sqrt(r);
	return r;
}

int ReadBMP(string strFile, int &size,int &width,int &height, double *& pixels, char *& head) {
	
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
	return 0;
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
	return 0;
}

__device__ double distance(double x1_x,double x1_y,double x1_z,double x2_x,double x2_y,double x2_z){
	return sqrt((x2_x-x1_x)*(x2_x-x1_x)+(x2_y-x1_y)*(x2_y-x1_y)+(x2_z-x1_z)*(x2_z-x1_z));
}


__global__ void Assign_center(double *data, unsigned int *assign, double *center,double *all_dis){
	//get idx for this datapoint
	int idx=blockIdx.x*blockDim.x+threadIdx.x;
	__shared__ double shared_center[c_num*3];
	if (threadIdx.x<d_k){
		shared_center[DIM*threadIdx.x]=center[DIM*threadIdx.x];
		shared_center[DIM*threadIdx.x+1]=center[DIM*threadIdx.x+1];
		shared_center[DIM*threadIdx.x+2]=center[DIM*threadIdx.x+2];
	}
	//if (threadIdx.x==0)
	//	all_dis[blockIdx.x]=0;
	if (idx<2)
		all_dis[idx]=0;
	__syncthreads();

	if (idx >= d_size) return;

	double min_distance = INFINITY;
	unsigned int mycenter = 0;
	dim3 mydata(data[DIM*idx],data[DIM*idx+1],data[DIM*idx+2]);

	for(int i = 0; i<d_k;i++)
	{
		double mydistance = distance(mydata.x,mydata.y,mydata.z,shared_center[i*DIM],shared_center[i*DIM+1],shared_center[i*DIM+2]);
		if(mydistance < min_distance){
			min_distance = mydistance;
			mycenter=i;
		}
	}
	assign[idx]=mycenter;
	//atomicAdd(&all_dis[int(blockIdx.x)],min_distance);
	atomicAdd(&all_dis[0],min_distance);
}

__global__   void Update_center(double *data, unsigned int  *assign,double *center_dist,double *center_num){
	int idx=blockIdx.x*blockDim.x+threadIdx.x;
	__shared__ double local_center_dist[c_num*3];
	__shared__ double local_center_num[c_num];	
	if (threadIdx.x<d_k){
		local_center_dist[DIM*threadIdx.x]=0;
		local_center_dist[DIM*threadIdx.x+1]=0;
		local_center_dist[DIM*threadIdx.x+2]=0;
		local_center_num[threadIdx.x]=0;
	}
	__syncthreads();

	if (idx >= d_size) return;

	int  mycenter=assign[idx];
	atomicAdd(&local_center_num[mycenter],1);
	atomicAdd(&local_center_dist[DIM*mycenter],data[DIM*idx]);
	atomicAdd(&local_center_dist[DIM*mycenter+1],data[DIM*idx+1]);
	atomicAdd(&local_center_dist[DIM*mycenter+2],data[DIM*idx+2]);	
	__syncthreads();

	if(threadIdx.x==0){
		for(int i=0;i<d_k;i++){
			atomicAdd(&center_num[i],local_center_num[i]);
			atomicAdd(&center_dist[DIM*i],local_center_dist[DIM*i]);
			atomicAdd(&center_dist[DIM*i+1],local_center_dist[DIM*i+1]);
			atomicAdd(&center_dist[DIM*i+2],local_center_dist[DIM*i+2]);				
		}	
	}

}

void KMeans(uint32_t size, uint32_t k, double * pixels) {
	uint32_t i; 
	uint32_t j;
	uint32_t c;
	double *center;
	uint32_t *clst;
	//uint32_t *clst_size;
	double min_cr = 50, max_cr = 200;
	double old_sum_dis, new_sum_dis;
	double assign_time=0;
	double update_time=0;
	cudaEvent_t start1;
	cudaEvent_t stop1;
	cudaEventCreate(&start1);
	cudaEventCreate(&stop1);
	float tmptime;

	center = (double *)malloc(k * DIM * sizeof(double));
	double *center_num = (double *)malloc(k * sizeof(double));
  srand(0);
  for (i = 0; i < k; ++i)
    for (j = 0; j < DIM; ++j) {
		  center[i * DIM + j] = min_cr + (double)(rand() * (max_cr - min_cr)) / (double)RAND_MAX;
	}

	//clst_size = (uint32_t *)malloc(k * sizeof(uint32_t));

	clst = (uint32_t *)malloc(size * sizeof(uint32_t));
	
	old_sum_dis = INF;
  	int iter = 0;

	int numBlocks = ceil(size/BLOCK_SIZE);
	double *d_pixels,*d_center,*d_center_num;
	unsigned int *d_clst;
	cudaMalloc(&d_pixels, sizeof(double)*size*DIM);	
	cudaMalloc(&d_clst, sizeof(unsigned int)*size);	
	cudaMalloc(&d_center, sizeof(double)*k*DIM);	
	cudaMalloc(&d_center_num, sizeof(double)*k);	

	double *test;
	double *local_test=(double *)malloc(sizeof(double)*size);
	cudaMalloc(&test, sizeof(double)*size);	

/*
	double *local_all_distance=(double *)malloc(sizeof(double)*numBlocks);	
	double *all_dis;
	cudaMalloc(&all_dis, sizeof(double)*numBlocks);	
*/
	double *local_all_distance=(double *)malloc(sizeof(double)*2);	
	double *all_dis;
	cudaMalloc(&all_dis, sizeof(double)*2);	

	cudaMemcpyToSymbol(d_size, &size, sizeof(unsigned int), 0, cudaMemcpyHostToDevice);
	cudaMemcpyToSymbol(d_k, &k, sizeof(unsigned int), 0, cudaMemcpyHostToDevice);

	cudaMemcpy(d_pixels,pixels, sizeof(double)*size*DIM,cudaMemcpyHostToDevice);
	cudaMemcpy(d_center,center, sizeof(double)*k*DIM,cudaMemcpyHostToDevice);		

	while (1) {
		iter ++;
    //#ifdef DEBUG
		//cout << "k-means iteration: " << iter  << ", total distance: " << old_sum_dis << endl;
    //#endif
    new_sum_dis = 0;

    cudaEventRecord(start1, NULL);
    Assign_center<<<numBlocks, BLOCK_SIZE>>>(d_pixels, d_clst, d_center,all_dis);
	cudaEventRecord(stop1, NULL);
	cudaEventSynchronize(stop1);
	cudaEventElapsedTime(&tmptime, start1, stop1);

    assign_time+=tmptime;

		cudaMemcpy(clst,d_clst, sizeof(unsigned int)*size,cudaMemcpyDeviceToHost);      		
			// very time-consuming
		//cudaMemcpy(local_all_distance,all_dis, sizeof(double)*numBlocks,cudaMemcpyDeviceToHost);
cudaMemcpy(local_all_distance,all_dis, sizeof(double)*2,cudaMemcpyDeviceToHost);

//		for(i=0;i< numBlocks;i++)
for(i=0;i< 1;i++)
			new_sum_dis+=local_all_distance[i];

		if (old_sum_dis - new_sum_dis < THRESHOLD) break;

		old_sum_dis = new_sum_dis;

		cudaMemset(d_center_num,0,k*sizeof(double));
		//cudaMemset(all_dis,0,numBlocks*sizeof(double));
cudaMemset(all_dis,0,2*sizeof(double));

		cudaMemset(d_center,0,DIM*k*sizeof(double));		

		cudaEventRecord(start1, NULL);
    		Update_center<<<numBlocks, BLOCK_SIZE>>>(d_pixels, d_clst, d_center,d_center_num);
		cudaEventRecord(stop1, NULL);
		cudaEventSynchronize(stop1);
		cudaEventElapsedTime(&tmptime, start1, stop1);
		update_time+=tmptime;

		cudaMemset(d_clst,0,size*sizeof(unsigned int));
		cudaMemcpy(center,d_center, sizeof(double)*k*DIM,cudaMemcpyDeviceToHost);  
		cudaMemcpy(center_num,d_center_num, sizeof(double)*k,cudaMemcpyDeviceToHost);  

		for (i = 0; i < k; ++i){
			center[i*DIM]/=center_num[i];
			center[i*DIM+1]/=center_num[i];
			center[i*DIM+2]/=center_num[i];
		}
		cudaMemcpy(d_center,center, sizeof(double)*k*DIM,cudaMemcpyHostToDevice);  		
	}

	for (i = 0; i < size; ++i) {
		c = clst[i];
		for (j = 0; j < DIM; ++j)
			pixels[i * DIM + j] = center[c * DIM + j];
	}
    cout<<"assign center time:"<< assign_time/1000 <<" s"<<endl;
    cout<<"update center time:"<< update_time/1000 <<" s"<<endl;
/*	free(clst_size);
	free(center);
	free(clst);*/
}

int main(int argc, char ** argv) {
	double * pixels;
	char * head;	
	string inputFile = "input.bmp";
	string outputFile = "kmeans_p.bmp";
	int size = 0;
	int width=0;
	int height=0;
	int k = 100;	
	
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
	if (argc>=2) {
		inputFile = argv[1];
	}
	if (argc>=3) {
		outputFile = argv[2];
  	}
	ReadBMP(inputFile, size, width, height, pixels, head);

    double cstart = currentSeconds();
	KMeans(size, c_num, pixels);
    double cend = currentSeconds();
    cout<<"total k-means time: "<< (cend-cstart) << " s" <<endl;   
	WriteBMP(outputFile, size, pixels, head);
	free(pixels);
	free(head);
	// cout << "K-means done." << endl;

}
