#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <string.h>
#include <map>
#include <cmath>
#include <vector>
#include <queue>
#include <omp.h>
#include <mpi.h>
#include "../util/cycletimer.h"

//#define RAND_MODE
#define MASTER 0
#define NUM_THREADS 16
//#define ADVANCE

using namespace std;

struct KStruct {
	char a, b, c;
	KStruct(char a_, char b_, char c_): a(a_), b(b_), c(c_) {}
	friend bool operator < (const struct KStruct &k1, const struct KStruct &k2);
};
 
inline bool operator < (const struct KStruct &k1, const struct KStruct &k2) {
	return k1.a < k2.a || (k1.a == k2.a && k1.b < k2.b) || (k1.a == k2.a && k1.b == k2.b && k1.c < k2.c) ;
}

int cmp(const pair<KStruct, int>& x, const pair<KStruct, int>& y)    
{    
    return x.second > y.second;    
}

void get_pic(string strFile, int &width, int &height, int &n, int &mysize, int &avgSize, 
     char** &pixels, char* &head, bool* &visited, int *&id, int* &size, char** &color, bool mode = 0){
	FILE *fin ;
	fin=fopen(strFile.c_str(),"rb");
	
	if (fin == NULL){
    cout<<"file open error!"<<endl;
  }
  int rank = 0, p = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (mode == 0) {
		//check file type
		short bfType;
		fread(&bfType,1,sizeof(short),fin);
		if(0x4d42!=bfType)
		{
			cout<<"the file is not a bmp file!"<<endl;
		}
		//get the n of BMP 
		fseek(fin,18,SEEK_SET) ;
		fread(&width,1,sizeof(int),fin);
		fread(&height,1,sizeof(int),fin);
		n = width * height ;
		//check the number of colors
		fseek(fin,28,SEEK_SET) ;
		unsigned short colors ;
		fread(&colors,1,sizeof(unsigned short),fin);
		if( colors != 24 )
		{
			cout << "The length of colors must be 24 bits" << endl ;
		}
		//get the file header
		fseek(fin,0,SEEK_SET);
    head = (char *)malloc(54 * sizeof(char));
		fread(head,54,sizeof(char),fin);
		
    avgSize = height / p;
	  if (height % p != 0) ++avgSize;
    avgSize = avgSize * width;
	  if (rank == p-1)
 		  mysize = n - (p - 1) * avgSize;
    else
      mysize = avgSize;
    // cout << width << ' ' << height << ' ' << avgSize << endl;
    pixels = (char **)malloc(avgSize * sizeof(char *));
		visited = (bool *) malloc(avgSize * sizeof(bool));
    id = (int *) malloc(avgSize * sizeof(int) * NUM_THREADS);
    size = (int *) malloc(avgSize * sizeof(int) * NUM_THREADS);
    color = (char **)malloc(avgSize * sizeof(char *) * NUM_THREADS);  
	  for (int i = 0; i < avgSize; i ++)
      pixels[i] = (char *)malloc(3 * sizeof(char));
      
    //read the pixels
	  fseek(fin,54+rank*avgSize*3,SEEK_SET);
	  for (int i = 0; i < mysize; i ++) {
		  fread(pixels[i], 3, sizeof(char), fin) ;
	  }
	} else {
    pixels = (char **)malloc(avgSize * sizeof(char *));
    for (int i = 0; i < avgSize; i ++)
      pixels[i] = (char *)malloc(3 * sizeof(char));
    if (rank != p-1) {
		  fseek(fin,54+height/p/2*width*3+rank*avgSize*3,SEEK_SET);
	    for (int i = 0; i < mysize; i ++) {
		    fread(pixels[i], 3, sizeof(char), fin) ;
	    }
    }
	}
  fclose(fin);
}

void update_pixels(int n, char** pixels, char** color, int *id)
{	
	#pragma omp parallel for num_threads(NUM_THREADS)
  for (int i = 0 ; i < n ; ++i) {
	    for(int k = 0 ; k < 3; k ++ ) {
            pixels[i][k] = color[id[i]][k];
	    }
	}
}
	
void put_pic(string file_out, int n, char* &head, char* &pixels)
{	
	FILE *fpw ;
	if((fpw=fopen(file_out.c_str(),"wb"))==NULL) {
		cout<<"create the bmp file error!"<<endl;
	}
	fwrite(head,sizeof(char),54,fpw);

	for (int i = 0 ; i < n ; ++i) {
 	   for(int k = 0 ; k < 3; k ++ ) {
            fwrite(&pixels[i*3+k], sizeof(char), 1, fpw);
     }
	}
	fclose(fpw);
}

bool same(char * a, char * b) {
    if(a[0] == b[0] && a[1] == b[1] && a[2] == b[2])
        return true;
    else
        return false;
}

int dx[4] = {0, -1, 0, 1};
int dy[4] = {-1, 0, 1, 0};

void bfs(int threshold, int width, int height, int n, char** pixels, bool* visited, int *id, int* size, char** color) {
	/* 
  ** find the colorblock id each pixel belongs to,
  ** and get the color and the size of each colorblock
  */
  height = n/width;

#pragma omp parallel num_threads(NUM_THREADS)
{
  bool * visited = (bool *) calloc(n, sizeof(bool));
  int p = omp_get_thread_num();
  int cnt = p*n;
  #pragma omp for
  for (int i = 0; i < height; i ++) {
      for (int j = 0; j < width; j ++) {
          int index = i * width + j;
          if (!visited[index]) {         
              visited[index] = true;
              color[cnt] = pixels[index];
              size[cnt] = 1;
              id[index] = cnt;
              queue<int> qx, qy ;
              qx.push(i);
              qy.push(j);
              while (!qx.empty()) {
                  int x = qx.front();
                  qx.pop();
                  int y = qy.front();
                  qy.pop();
                  for (int u = 0; u < 4; u ++) {
                      int x_t = x + dx[u];
                      int y_t = y + dy[u];
                      if (x_t >= height || x_t < 0 || y_t >= width || y_t < 0)
                          continue;
                      int index_t = x_t * width + y_t;
                      if (same(pixels[index_t], pixels[index]) && !visited[index_t]) {
                          visited[index_t] = true;
                          qx.push(x_t);
                          qy.push(y_t);
                          size[cnt] += 1;
                          id[index_t] = cnt;
                      }
                  }
              }
              cnt += 1;
          }
      }
  }
}
  
	for (int i = 0; i < height; i ++)
      for (int j = 0; j < width; j ++) 
	        visited[i*width+j] = false;
 
  //recolor
#pragma omp parallel num_threads(NUM_THREADS)
{
  bool * visited = (bool *) calloc(n, sizeof(bool));
  int p = omp_get_thread_num();
  int cnt = p*n;
  #pragma omp for
  
  for (int i = 0; i < height; i ++) {
      for (int j = 0; j < width; j ++) {
          int index = i * width + j;
          if (!visited[index]) {
              map<KStruct, int> neighbor_cnt;
              queue<int> qx, qy ;
              qx.push(i);
              qy.push(j);
              while (!qx.empty()) {
                  int x = qx.front();
                  qx.pop();
                  int y = qy.front();
                  qy.pop();
                  for (int u = 0; u < 4; u ++) {
                      int x_t = x + dx[u];
                      int y_t = y + dy[u];
                      if (x_t >= height || x_t < 0 || y_t >= width || y_t < 0)
                          continue;
                      int index_t = x_t * width + y_t;
                      if (!same(pixels[index_t], pixels[index])) {
                          KStruct key(pixels[index_t][0], pixels[index_t][1], pixels[index_t][2]);
            							if(neighbor_cnt.find(key) != neighbor_cnt.end())
            								neighbor_cnt[key] += sqrt(size[id[index_t]]);
            							else
            								neighbor_cnt[key] = sqrt(size[id[index_t]]);
                      } else if (!visited[index_t]) {
                          visited[index_t] = true;
                          qx.push(x_t);
                          qy.push(y_t);
                      }
                  }    
              }  
              #ifdef RAND_MODE
              if (rand()/(double)RAND_MAX < sqrt(sqrt( size[cnt]/(double)n )) ) {
                  cnt += 1;
                  continue;
              }
              #endif
              if (size[cnt] < threshold) {
                  vector<pair<KStruct, int> > score_vec(neighbor_cnt.begin(), neighbor_cnt.end());  
        					sort(score_vec.begin(), score_vec.end(), cmp);
        					KStruct* key = & score_vec[0].first;
        					char * tempcolor = (char *) malloc(3*sizeof(char));
        					tempcolor[0] = key->a;
        					tempcolor[1] = key->b;
        					tempcolor[2] = key->c;
      					  color[cnt] = tempcolor;
              }
              cnt += 1;
          }
      }
  }
}
}

int main(int argc, char ** argv)
{
	int width, height;
  char ** pixels, **pixels2;
  char * pixs;
  char * all_pixs, *all_pixs2;
  bool * visited;
  int * id;
  int * size;
  char ** color;
  char * head;
  int n;
  int rank = 0, p = 0;
  int mysize, avgSize;
  
  int blockSize = 15;
	int iter = 5;
  double ctime = 0;
	string inputFile = "kmeans.bmp";
	string outputFile = "colorBlock_omp.bmp";
	
  if (argc>=2) {
    inputFile = argv[1];
  }
  if (argc>=3) {
    outputFile = argv[2];
  }
  if (argc>=4) {
    blockSize=atoi(argv[3]);
  }
  if (argc>=5) {
    iter=atoi(argv[4]);
  }
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	//read parameters
  if (rank == MASTER) {
    cout << "-----------------colorBlock part------------------" << endl;
  }
  get_pic(inputFile, width, height, n, mysize, avgSize, pixels, head, visited, id, size, color);
  double cstart = currentSeconds();
	for (int i = 1; i <= iter; i ++) {
		for (int j = 0; j < mysize ; j ++) {
			visited[j] = false;
			id[j] = 0;
			size[j] = 0;
			color[j] = NULL;
		}
		 
		bfs(blockSize, width, height, mysize, pixels, visited, id, size, color);
    update_pixels(mysize, pixels, color, id);
	}
 
#ifdef ADVANCE
  get_pic(inputFile, width, height, n, mysize, avgSize, pixels2, head, visited, id, size, color, 1);
  for (int i = 1; i <= iter; i ++) {
		for (int j = 0; j < mysize ; j ++) {
			visited[j] = false;
			id[j] = 0;
			size[j] = 0;
			color[j] = NULL;
		}
		 
		bfs(blockSize, width, height, mysize, pixels2, visited, id, size, color);
    update_pixels(mysize, pixels2, color, id);
	}
#endif
  pixs = (char *)malloc((avgSize * 3) * sizeof(char));
  for (int i = 0; i < avgSize; i ++)
    for (int j = 0; j < 3; j ++)
      pixs[i*3+j] = pixels[i][j];

  //gather the color for each pixel
  if (rank != MASTER)
    MPI_Gather(pixs, avgSize*3, MPI_CHAR, all_pixs, avgSize*3, MPI_CHAR, MASTER, MPI_COMM_WORLD);
  else {	
    all_pixs = (char *)malloc((avgSize * p * 3) * sizeof(char));
	  MPI_Gather(pixs, avgSize*3, MPI_CHAR, all_pixs, avgSize*3, MPI_CHAR, MASTER, MPI_COMM_WORLD);
  }

#ifdef ADVANCE  
  for (int i = 0; i < avgSize; i ++)
    for (int j = 0; j < 3; j ++)
      pixs[i*3+j] = pixels2[i][j];  
  //gather the color for each pixel
  if (rank != MASTER)
    MPI_Gather(pixs, avgSize*3, MPI_CHAR, all_pixs2, avgSize*3, MPI_CHAR, MASTER, MPI_COMM_WORLD);
  else {	
    all_pixs2 = (char *)malloc((avgSize * p * 3) * sizeof(char));
	  MPI_Gather(pixs, avgSize*3, MPI_CHAR, all_pixs2, avgSize*3, MPI_CHAR, MASTER, MPI_COMM_WORLD);
  }

  if (rank == MASTER) {
    int diff = height/p/2*width*3;
    for (int x = height/p/2; x < height - height/p/2; x ++) {
      if ( x%(height/p) < height/p/4 || x%(height/p) > height/p/4 * 3) {
        for (int i = 3*x*width; i < 3*x*width + 3*width; i ++)
          all_pixs[i] = all_pixs2[i-diff];
      }
    }
  }
#endif  
  
  if (rank == MASTER) {
    ctime = currentSeconds() - cstart; 
    put_pic(outputFile, n, head, all_pixs);
		free(all_pixs);
#ifdef ADVANCE
    free(all_pixs2);
#endif    
	}
   	
  for (int i = 0; i < avgSize; i++) {
		free(pixels[i]);
#ifdef ADVANCE
    free(pixels2[i]);
#endif
  }
  free(head);
  free(pixs);
	free(pixels);
#ifdef ADVANCE
  free(pixels2);
#endif
  free(color);
  free(visited);
  free(id);
  free(size);
  MPI_Finalize();
	if (rank == MASTER) {
    cout << "color block pruning time: " << ctime << " s" <<endl;  
    cout << "-----------------colorBlock done------------------" << endl;
  }
}
