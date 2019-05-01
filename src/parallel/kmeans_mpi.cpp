#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include "../util/cycletimer.h"

//#define OPENMP1
//#define OPENMP2
//#define OPENMP3
//#define OPENMP4
//#define OPENMP5
#define INIT1
//#define DEBUG

using namespace std;

typedef unsigned int uint32_t;

#define INF (1e20)
#define THRESHOLD	(1e-2)
#define DIM 3
#define MASTER 0
#define MAX_K 256
#define MAX_SIZE 6000*4000

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

int ReadBMPHead(string strFile, int &size, char *& head) {
  
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

	fclose(fin);
  return 1;
}

int ReadBMP(string strFile, int &avgSize, int &mysize, double *& pixels) { 
  FILE *fin ;
	fin=fopen(strFile.c_str(),"rb");
	
	//check file pointer
	if(fin == NULL) {
		cout<<"file open error!"<<endl;
		return 0;
	}

  int rank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//read the pixels
	fseek(fin,54+rank*avgSize*DIM,SEEK_SET);
	for (int i = 0; i < mysize; i ++) {
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
  int p = 0;
  int rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
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

void KMeans(uint32_t mysize, uint32_t avgSize, uint32_t size, uint32_t k, double * pixels) {
	uint32_t i; 
	uint32_t j;
	uint32_t c;
	double *center;
	uint32_t *clst_id;
	uint32_t *clst_size;
  uint32_t *clst_id_buf = NULL;
  double *center_buf;
  uint32_t *clst_size_buf;
	double old_sum_dis, new_sum_dis, new_sum_dis_buf;
	double min_dis, dis;
	double min_cr = 50, max_cr = 200;
  int iter = 0;
  bool cont;
  int p = 0;
  int rank = 0;
  double ustart, astart, utime = 0, atime = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //malloc centers
  center = (double *)calloc(k * DIM, sizeof(double));
	clst_size = (uint32_t *)calloc(k, sizeof(uint32_t));
  clst_id = (uint32_t *)malloc(avgSize * sizeof(uint32_t));
 	
	center_buf = (double *)malloc(k * DIM * sizeof(double));
	clst_size_buf = (uint32_t *)calloc(k, sizeof(uint32_t));	
  #ifdef INIT1
  //randomly initialize centers
  if(rank==MASTER) {
    srand(0);
    for (i = 0; i < k; ++i)
      for (j = 0; j < DIM; ++j) {
			  center[i * DIM + j] = min_cr + (double)(rand() * (max_cr - min_cr)) / (double)RAND_MAX;
		  }
  }
  #else
  //round robin initialize centers
	c = 0;
	for (i = 0; i < mysize; ++i) {
		for (j = 0; j < DIM; ++j)
			center[c * DIM + j] += pixels[i];
		clst_size[c] ++;
		c = (c + 1) % k;
	}
  //gather center sum and calculate center location
  MPI_Reduce(center, center_buf, k * DIM, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
	MPI_Reduce(clst_size, clst_size_buf, k, MPI_INT32_T, MPI_SUM, MASTER, MPI_COMM_WORLD);
	if (rank == MASTER) { 
    #pragma omp parallel for private(j)
		for (i = 0; i < k; ++i) {
			for (j = 0; j < DIM; ++j) {
				center[i * DIM + j] = center_buf[i * DIM + j] / clst_size_buf[i];
      }
    }
	}  
  #endif

	old_sum_dis = INF;
	while (1) {
    astart = currentSeconds();
    //broadcast centers
		MPI_Bcast(center, k * DIM, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

    //find new center for each pixel
    new_sum_dis = 0;
    #ifdef OPENMP1
    #pragma omp parallel for private(j, dis, min_dis, c) reduction(+: new_sum_dis)
    #endif
		for (i = 0; i < mysize; ++i) {
			min_dis = INF;
      c = 0;
			for (j = 0; j < k; ++j) {
				dis = CalcDistance(&pixels[i * DIM], &center[j * DIM]);
				if (dis < min_dis) {
					min_dis = dis;
					c = j;
				}
			}
			clst_id[i] = c;
			new_sum_dis += min_dis;
		}
    atime += currentSeconds() - astart;

    //check if the calculation can terminate
    MPI_Reduce(&new_sum_dis, &new_sum_dis_buf, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
		if (rank == MASTER) {
			cont = (old_sum_dis - new_sum_dis_buf > THRESHOLD);
			old_sum_dis = new_sum_dis_buf;
      iter ++;
      #ifdef DEBUG
      cout << "k-means iteration: " << iter  << ", total distance: " << old_sum_dis << endl;
      #endif
		}
		MPI_Bcast(&cont, 1, MPI_C_BOOL, MASTER, MPI_COMM_WORLD);	
		if (!cont) break;
   
    ustart = currentSeconds();
    //clear centers (local copy)
    #ifdef OPENMP3
    #pragma omp parallel for private(j)
    #endif
		for (i = 0; i < k; ++i) {
			for (j = 0; j < DIM; ++j)
				center[i * DIM + j] = 0;
			clst_size[i] = 0;
		}

    //add up to centers (local copy)
    #ifdef OPENMP2
    #pragma omp parallel
    {
        double center_temp[MAX_K*DIM] = {0};
        int clst_size_temp[MAX_K] = {0};
        #pragma omp for private(j, c)
    		for (i = 0; i < mysize; ++i) {
    			c = clst_id[i];
    			for (j = 0; j < DIM; ++j)
    				center_temp[c * DIM + j] += pixels[i * DIM + j];
    			clst_size_temp[c] += 1;
    		}
        #pragma omp critical
        {
            for (i = 0; i < k*DIM; i++) {
                center[i] += center_temp[i];
            }
            for (i = 0; i < k; i++) {
                clst_size[i] += clst_size_temp[i];
            }
        }
    }
    #else
    for (i = 0; i < mysize; ++i) {
			c = clst_id[i];
			for (j = 0; j < DIM; ++j)
				center[c * DIM + j] += pixels[i * DIM + j];
			clst_size[c] += 1;
		}
    #endif
    //gather sums of clusters and calculate new cluster centers
    MPI_Reduce(center, center_buf, k * DIM, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
		MPI_Reduce(clst_size, clst_size_buf, k, MPI_INT32_T, MPI_SUM, MASTER, MPI_COMM_WORLD);
		if (rank == MASTER) { 
      #ifdef OPENMP4
      #pragma omp parallel for private(j)
      #endif
			for (i = 0; i < k; ++i) {
				for (j = 0; j < DIM; ++j) {
					center[i * DIM + j] = center_buf[i * DIM + j] / clst_size_buf[i];
        }
      }
		}
    utime += currentSeconds() - ustart;
	} //end k-means loop
 
  //gather the cluster id for each pixel
  if (rank != MASTER)
    MPI_Gather(clst_id, avgSize, MPI_UINT32_T, clst_id_buf, avgSize, MPI_UINT32_T, MASTER, MPI_COMM_WORLD);
  else {	
    clst_id_buf = (uint32_t *)malloc((avgSize * p) * sizeof(uint32_t));
	  MPI_Gather(clst_id, avgSize, MPI_UINT32_T, clst_id_buf, avgSize, MPI_UINT32_T, MASTER, MPI_COMM_WORLD);
  
    //update the color of each pixel
    #ifdef OPENMP5
    #pragma omp parallel for private(j, c)
    #endif
		for (i = 0; i < size; ++i) {
			c = clst_id_buf[i];
			for (j = 0; j < DIM; ++j)
				pixels[i * DIM + j] = center[c * DIM + j];
		}
		free(clst_id_buf);
    cout << "update time: " << utime << " s" <<endl;  
    cout << "assign time: " << atime << " s" <<endl;
	}
	free(center_buf);
	free(clst_size_buf);
	free(clst_size);
	free(center);
	free(clst_id);
  MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char ** argv) {
	double * pixels;
	char * head;	
	string inputFile = "input.bmp";
	string outputFile = "kmeans_p.bmp";
  double cstart = 0, cend;
  int size = 0;
  int avgSize = 0;
	int mysize = 0;
	int k = 100;	
  int p = 0;
  int rank = 0;
  
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (argc>=2) {
    inputFile = argv[1];
  }
  if (argc>=3) {
    outputFile = argv[2];
  }
  if (argc>=4) {
    k=atoi(argv[3]);
  }
	//read parameters
  if (rank == MASTER) {
    cout << "-----------------K-means part------------------" << endl;
    ReadBMPHead(inputFile, size, head);
  }
  if (k > MAX_K) {
    if (rank == MASTER)
      cout << "k exceeds the limitation: " << MAX_K << endl;
    MPI_Finalize();
    return 0;
  }
  if (size > MAX_SIZE) {
    if (rank == MASTER)
      cout << "image size exceeds the limitation: " << MAX_SIZE << endl;
    MPI_Finalize();
    return 0;
  }
  
  //broadcast and calculate sizes
  MPI_Bcast(&k, 1, MPI_INT32_T, MASTER, MPI_COMM_WORLD);
 	MPI_Bcast(&size, 1, MPI_INT32_T, MASTER, MPI_COMM_WORLD);
	avgSize = size / p;
	if (size % p != 0) ++avgSize;
	if (rank == p - 1)
		mysize = size - (p - 1) * avgSize;
  else
    mysize = avgSize;
   
  //allocate memory for pixels
  if (rank == MASTER) {
    pixels = (double *)malloc(size * DIM * sizeof(double));
  } else {
    pixels = (double *)malloc(avgSize * DIM * sizeof(double));
  }
  
  //read image
  if (rank == MASTER) {
    cstart = currentSeconds();
  }
  ReadBMP(inputFile, avgSize, mysize, pixels);
  if (rank == MASTER) {
    cend = currentSeconds();
    cout << "read image time: " << cend-cstart << " s" <<endl; 
    cstart = currentSeconds();
  }

  //perform k-means
	KMeans(mysize, avgSize, size, k, pixels);
  
  if (rank == MASTER) {
    cend = currentSeconds();
    cout << "K-means time: " << cend-cstart << " s" <<endl;  
    
    //write image
    cstart = currentSeconds(); 
    WriteBMP(outputFile, size, pixels, head);
    cend = currentSeconds();
    cout << "write image time: " << cend-cstart << " s" <<endl;  
   	free(head);
  }
	free(pixels);
  MPI_Finalize();
	if (rank == MASTER) {
    cout << "-----------------K-means done------------------" << endl;
  }
}
