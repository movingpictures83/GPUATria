#ifndef GPUATRIAPLUGIN_H
#define GPUATRIAPLUGIN_H

#include "Plugin.h"
#include "PluginProxy.h"
#include <string>
#include <vector>
#include <set>
//#define GSIZE 126
#define BLOCK_SIZE 256
#define _DTH cudaMemcpyDeviceToHost
#define _HTD cudaMemcpyHostToDevice


class GPUATriaPlugin : public Plugin 
{
  public:    
  std::string toString(){return "GPUATria";}
  ~GPUATriaPlugin();
  void input(std::string file);
  void run();
  void output(std::string file);

  private:
     float* OrigGraph;
     std::string* bacteria;//[GSIZE];
     float* dG;
     float* dTable;
     float* dPay;
     int* dSpots;
     char* dMark;
     std::vector<float> U; 
     std::set<int> ks;  
     float* H_pay;
     int GSIZE;
     //void _GPU_Floyd(float *H_G, /*int *H_Gpath,*/ float *H_pay, const int N);
     void getIndices(int bac1, int bac2, unsigned long* pos_pos, unsigned long* pos_neg, unsigned long* neg_pos, unsigned long* neg_neg);
     void cudaErrorCheck(cudaError_t err);
};

//__global__ void _Wake_GPU(int reps);
__global__ void _GPU_Floyd_kernel(int k, float *G, const int* __restrict__ spots,/*int *P,*/ char* m, const int N, const int a, const int _2GSIZE);
__global__ void _GPU_Pay_kernel(float* D, float* P, const int* __restrict__ spots, int N, int GSIZE);
__global__ void _GPU_Copy_kernel(float* dst, float* src, unsigned long N);
__global__ void _GPU_Triad_kernel(float* G, int maxnode, int GSIZE);
__global__ void _GPU_Sweep_kernel(float* G, int N);
__global__ void _GPU_Spots_kernel(int* S, int GSIZE);
//__constant__ int _GPU_GSIZE;
//__constant__ int _GPU_2GSIZE;

#endif

