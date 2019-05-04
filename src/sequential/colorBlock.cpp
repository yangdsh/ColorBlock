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
#include "../util/cycletimer.h"

//#define RAND_MODE

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

void get_pic(string strFile, int &width, int &height, int &n, char* &head, char** &pixels, bool* &visited, int *&id, int* &size, char** &color){
	FILE *fpi ;
	fpi=fopen(strFile.c_str(),"rb");
	
	if(fpi != NULL){
		//check file type
		short bfType;
		fread(&bfType,1,sizeof(short),fpi);
		if(0x4d42!=bfType)
		{
			cout<<"the file is not a bmp file!"<<endl;
		}
		//get the n of BMP 
		fseek(fpi,18,SEEK_SET) ;
		fread(&width,1,sizeof(int),fpi);
		fread(&height,1,sizeof(int),fpi);
		n = width * height ;
		//check the number of colors
		fseek(fpi,28,SEEK_SET) ;
		unsigned short colors ;
		fread(&colors,1,sizeof(unsigned short),fpi);
		if( colors != 24 )
		{
			cout << "The length of colors must be 24 bits" << endl ;
		}
		//get the file header
		fseek(fpi,0,SEEK_SET);
    head = (char *)malloc(54 * sizeof(char));
		fread(head,54,sizeof(char),fpi);
		fseek(fpi,54,SEEK_SET);
		
		pixels = (char **)malloc(n * sizeof(char *));
		visited = (bool *) malloc(n * sizeof(bool));
    id = (int *) malloc(n * sizeof(int));
    size = (int *) malloc(n * sizeof(int));
    color = (char **)malloc(n * sizeof(char *));

		for( int i = 0 ; i < n ; i ++ )
		{
			visited[i] = false;
      pixels[i] = (char *)malloc(3 * sizeof(char));
			fread(pixels[i], 3, sizeof(char), fpi) ;
		}
		fclose(fpi);
	}
	else
	{
		cout<<"file open error!"<<endl;
	}
}

void update_pixels(int n, char** pixels, char** color, int *id)
{	
	for (int i = 0 ; i < n ; ++i) {
	    for(int k = 0 ; k < 3; k ++ ) {
            pixels[i][k] = color[id[i]][k];
	    }
	}
}
	
void put_pic(string file_out, int n, char* &head, char** &pixels)
{	
	FILE *fpw ;
	if((fpw=fopen(file_out.c_str(),"wb"))==NULL) {
		cout<<"create the bmp file error!"<<endl;
	}
	fwrite(head,sizeof(char),54,fpw);

	for (int i = 0 ; i < n ; ++i) {
	    for(int k = 0 ; k < 3; k ++ ) {
            fwrite(&pixels[i][k], sizeof(char), 1, fpw);
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
  int cnt = 0;
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
  
	for (int i = 0; i < height; i ++)
      for (int j = 0; j < width; j ++) 
	        visited[i*width+j] = false;
 
  //recolor
  cnt = 0;
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

int main(int argc, char ** argv)
{
	int width, height;
  char ** pixels;
  bool * visited;
  int * id;
  int * size;
  char ** color;
  char * head;
  static int n;
  
  int blockSize = 1500;
	int iter = 10;
  double ctime = 0;
	string inputFile = "kmeans.bmp";
	string outputFile = "colorBlock.bmp";
	
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
		
		cout << "blockSize: ";
		cin.getline(line, 100);
		if (strlen(line) > 0)
			blockSize = atoi(line);
		else
      cout << blockSize << endl;
		
		cout << "iterations: ";
		cin.getline(line, 100);
		if (strlen(line) > 0)
			iter = atoi(line);
    else
		  cout << iter << endl;
	}
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
	
  get_pic(inputFile, width, height, n, head, pixels, visited, id, size, color);
	for (int i = 1; i <= iter; i ++) {
		for (int j = 0; j < n ; j ++) {
			visited[j] = false;
			id[j] = 0;
			size[j] = 0;
			color[j] = NULL;
		}
		
    double cstart = currentSeconds(); 
		bfs(blockSize, width, height, n, pixels, visited, id, size, color);
    double cend = currentSeconds(); 
    ctime += cend - cstart;
    update_pixels(n, pixels, color, id);
	}
  put_pic(outputFile, n, head, pixels);
  	
  for (int i = 0; i < n; i++) {;
		free(pixels[i]);
  }
  free(head);
	free(pixels);
  free(color);
  free(visited);
  free(id);
  free(size);
  cout << "total color block pruning time: " << ctime << " s" <<endl;
}
