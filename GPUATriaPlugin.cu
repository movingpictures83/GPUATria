#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include "PluMA.h"
#include "PluginManager.h"
#include "GPUATriaPlugin.h"
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "csv_parser/csv_parser.cpp"

bool floatComp(float A, float B) {
   float EPSILON=1e-7;
   float diff = A - B;
   return (diff < EPSILON) && (-diff < EPSILON);
}


GPUATriaPlugin::~GPUATriaPlugin() {
   if (OrigGraph) free(OrigGraph);
   if (bacteria) delete bacteria;
   if (H_pay) free(H_pay);
   
}


void GPUATriaPlugin::getIndices(int bac1, int bac2, unsigned long* pos_pos, unsigned long* pos_neg, unsigned long* neg_pos, unsigned long* neg_neg) {
         if (bac1 > bac2) {
		int tmp = bac1;
		bac1 = bac2;
		bac2 = tmp;
         }
 
				 unsigned long row_pos = 2*(2*GSIZE*bac1 - bac1*(bac1-1));
				 unsigned long row_neg = row_pos + 2*(GSIZE-bac1);;
                                 
                                 *pos_pos = row_pos + 2*(bac2-bac1);
                                 *pos_neg = *pos_pos + 1;
                                 *neg_pos = row_neg + 2*(bac2-bac1);
                                 *neg_neg = *neg_pos + 1;
}

void GPUATriaPlugin::cudaErrorCheck(cudaError_t err) {
	 if(err!=cudaSuccess){printf("%s in %s at line %d\n",cudaGetErrorString(err),__FILE__,__LINE__);exit(1);}
}

void GPUATriaPlugin::input(std::string file) {
//const int NumBytes=(GSIZE*2)*(GSIZE*2)*sizeof(float);
                //host allocations to create Adjancency matrix and result matrices with path matrices
                //OrigGraph=(float *)malloc(NumBytes);//will be original Adjancency matrix, will NOT be changed
 
                const char field_terminator = ',';
                const char line_terminator  = '\n';
                const char enclosure_char   = '"';
   // File is in CSV format
  csv_parser file_parser;
                file_parser.set_skip_lines(1);
                file_parser.init(file.c_str());
                file_parser.set_enclosed_char(enclosure_char, ENCLOSURE_OPTIONAL);
                file_parser.set_field_term_char(field_terminator);
                file_parser.set_line_term_char(line_terminator);

                GSIZE = 0;
                while (file_parser.has_more_rows()) {
                   file_parser.get_row();
                   GSIZE++;
                }
                const int NumBytes=(GSIZE*2)*(GSIZE*2)*sizeof(float);
                OrigGraph=(float *)malloc(NumBytes);//will be original Adjancency matrix, will NOT be changed
                bacteria = new std::string[GSIZE];

                file_parser.init(file.c_str());
                //set<int> ks;
                unsigned int row_count = 0;
                unsigned long pos_pos, pos_neg, neg_pos, neg_neg;
        	while(file_parser.has_more_rows())
        	{
                	unsigned int i = 0;

                	csv_row row = file_parser.get_row();

			bacteria[row_count] = row[0];
			int count = 0;
			//cout << row_count << " " << row.size() << endl;
                	for (i = 1; i < row.size(); i++) {
 			      int bac1 = row_count;
			      int bac2 = i-1;
			      float weight = atof(row[i].c_str());
			      if (bac1 == bac2)
				weight = 1;
                              //if (row_count == 637) cout << "I: " << i << endl;
                              if (bac1 <= bac2) {
				 getIndices(bac1, bac2, &pos_pos, &pos_neg, &neg_pos, &neg_neg);
                                 if (weight > 0) {
				    OrigGraph[pos_pos] = weight;
				    OrigGraph[neg_neg] = weight;
				    OrigGraph[pos_neg] = 0;
				    OrigGraph[neg_pos] = 0;
				    count++;
                                 }
                                 else if (weight < 0) {
				    OrigGraph[pos_neg] = weight;
				    OrigGraph[neg_pos] = weight;
				    OrigGraph[pos_pos] = 0;
				    OrigGraph[neg_neg] = 0;
                                    count++;
                                 }
                              	else {
				    OrigGraph[pos_pos] = 0;
				    OrigGraph[pos_neg] = 0;
				    OrigGraph[neg_pos] = 0;
				    OrigGraph[neg_neg] = 0;
				}
                              }
                        }
			if (count != 0) {ks.insert(2*row_count); ks.insert(2*row_count+1);}
                	row_count++;
        	}
                /*while(file_parser.has_more_rows())
                {
                        unsigned int i = 0;

                        csv_row row = file_parser.get_row();
                        bacteria[row_count] = row[0];

                        for (i = 1; i < row.size(); i++) {
                              int bac1 = row_count;
                              int bac2 = i-1;
                              float weight = atof(row[i].c_str());
                              if (bac1 != bac2) {
                                 if (weight > 0) {
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2)] = weight;
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2+1)] = weight;
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2)] = 0;
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2+1)] = 0;
                                 }
                                 else if (weight < 0) {
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2)] = weight;
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2+1)] = weight;
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2)] = 0;
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2+1)] = 0;
                                 }
                                else {
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2)] = 0;
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2+1)] = 0;
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2)] = 0;
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2+1)] = 0;
                                }
                              }
                              else {
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2)] = 1; // Start these at 1, because they are starting verts.
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2+1)] = 1;
                                    OrigGraph[(bac1*2+1)*(2*GSIZE)+(bac2*2)] = 0;
                                    OrigGraph[(bac1*2)*(2*GSIZE)+(bac2*2+1)] = 0;
                              }
                        }
                        row_count++;
                }*/
  

}

void GPUATriaPlugin::run() {
         const int NumBytes=(GSIZE*2)*(GSIZE*2)*sizeof(float);
         const int NumBytesPay=(GSIZE)*sizeof(float);   // N by 2N.  2 paths from every i to every j
         const int _2GSIZE = 2*GSIZE;
                H_pay = (float *)malloc(NumBytesPay);
         U.resize(GSIZE, 0.0);
	 PluginManager::log("GPU Global Memory Allocation...");
         /////////////////////////////////////////////////////////////////////
         // GPU Memory Allocation
	 cudaErrorCheck(cudaMalloc((float **)&dG,NumBytes));
	 cudaErrorCheck(cudaMalloc((float **)&dTable, NumBytes));
         cudaErrorCheck(cudaMalloc((float**)&dPay,NumBytesPay));
         cudaErrorCheck(cudaMalloc((int**)&dSpots, GSIZE*2*GSIZE*2*sizeof(int)));
	 dim3 pathThreads((GSIZE*2+BLOCK_SIZE-1)/BLOCK_SIZE,GSIZE*2);
	 //dim3 pathThreads((GSIZE*GSIZE*2*2+BLOCK_SIZE-1)/BLOCK_SIZE,GSIZE*2);
         dim3 payThreads((GSIZE+BLOCK_SIZE-1)/BLOCK_SIZE);
         dim3 spotThreads((GSIZE*2*GSIZE*2+BLOCK_SIZE-1)/BLOCK_SIZE, GSIZE*2);
         //int numGrids = (GSIZE*2*(GSIZE+1)+BLOCK_SIZE-1)/BLOCK_SIZE;
	 //dim3 copyThreads((numGrids+65535+1)/65535, 65535);
	 unsigned long copyThreads = (GSIZE*2*(GSIZE+1)+BLOCK_SIZE-1)/BLOCK_SIZE;
         unsigned long triadThreads = (GSIZE*GSIZE*8+BLOCK_SIZE-1)/BLOCK_SIZE;
	 unsigned long sweepThreads = (GSIZE*2*(GSIZE+1)+BLOCK_SIZE-1)/BLOCK_SIZE;
         ////////////////////////////////////////////////////////////////////
	 cudaErrorCheck(cudaMemcpy(dG,OrigGraph,NumBytes,_HTD));

	 free(OrigGraph);
         int currentrank=1;
         PluginManager::log("Done.");
         PluginManager::log("Computing Spots...");
         _GPU_Spots_kernel<<<spotThreads, BLOCK_SIZE>>>(dSpots, GSIZE);
         cudaErrorCheck(cudaGetLastError()); 
         cudaErrorCheck(cudaThreadSynchronize());
         PluginManager::log("Done.");
         for (int a = 0; a < GSIZE; a++) {

                ///////////////////////////////////////////////////////////////////////////////////////////
                // GPU Floyd Algorithm
		_GPU_Copy_kernel<<<copyThreads, BLOCK_SIZE>>>(dTable, dG, GSIZE*2*(GSIZE+1)); 
                cudaErrorCheck(cudaGetLastError());
		cudaErrorCheck(cudaThreadSynchronize());
                // Note: k is the node you are going *through*.  Not the start node.
	        for(set<int>::iterator k=ks.begin();k!=ks.end();k++){//main loop
		   _GPU_Floyd_kernel<<<pathThreads,BLOCK_SIZE>>>(*k,dTable,dSpots,/*dW,*//*dP*/dMark,_2GSIZE,a, _2GSIZE);
                   cudaErrorCheck(cudaGetLastError());
		   cudaErrorCheck(cudaThreadSynchronize());
	        }

                _GPU_Pay_kernel<<<payThreads, BLOCK_SIZE>>>(dTable, dPay, dSpots, GSIZE, GSIZE);
                cudaErrorCheck(cudaGetLastError());
	        cudaErrorCheck(cudaMemcpy(H_pay,dPay,NumBytesPay,_DTH));
                ////////////////////////////////////////////////////////////////////////////////////////
                vector<int> maxnodes;
		int mnode = -1;
                float maxpay = -1;
                for (int i = 0; i < GSIZE; i++) {
                   if (fabs(H_pay[i]) > maxpay) {
                      mnode = i;
                      maxpay = fabs(H_pay[i]);
                   }
                }
		//if (maxpay == 0.5) break;
                //U[mnode] = H_pay[mnode];
		maxnodes.push_back(mnode);
		for (int i = 0; i < GSIZE; i++) {
                   //if ((i != mnode) && fabs(H_pay[i]) == maxpay) {
                   if ((i != mnode) && floatComp(fabs(H_pay[i]), maxpay)) {
                      maxnodes.push_back(i);
	           }
		}
                if (maxpay == 0)
                   break;
		for (int w = 0; w < maxnodes.size(); w++) {
                        int maxnode = maxnodes[w];
                PluginManager::log(std::string("Node with highest pay: "+bacteria[maxnode]+": "+std::to_string(H_pay[maxnode])));
                U[maxnode] = currentrank;//H_pay[maxnode];
		ks.erase(maxnode*2);
		ks.erase(maxnode*2+1);
		//cudaErrorCheck(cudaThreadSynchronize());
	        //cudaErrorCheck(cudaMemcpy(dG,OrigGraph,NumBytes,_HTD));
                _GPU_Triad_kernel<<<triadThreads, BLOCK_SIZE>>>(dG, maxnode, GSIZE);
                cudaErrorCheck(cudaGetLastError());
		cudaErrorCheck(cudaThreadSynchronize());

		cudaErrorCheck(cudaThreadSynchronize());
		_GPU_Sweep_kernel<<<sweepThreads, BLOCK_SIZE>>>(dG, GSIZE*2*(GSIZE+1));
                cudaErrorCheck(cudaGetLastError());
		cudaErrorCheck(cudaThreadSynchronize());
		}
		currentrank += maxnodes.size();

           }
        free(H_pay);

        /////////////////////////////////////////////////////////////
        // FREE GPU MEM
	cudaErrorCheck(cudaFree(dG));
	cudaErrorCheck(cudaFree(dTable));
	cudaErrorCheck(cudaFree(dPay));
       
        /////////////////////////////////////////////////////////////

}




void GPUATriaPlugin::output(std::string file) {
for (int i = GSIZE-1; i >= 0; i--)
           for (int j = 0; j < i; j++) {
              if (fabs(U[j]) > fabs(U[j+1])) {
                 float tmp = U[j];
                 U[j] = U[j+1];
                 U[j+1] = tmp;
                 string tmp2 = bacteria[j];
                 bacteria[j] = bacteria[j+1];
                 bacteria[j+1] = tmp2;
              }
           }
/*for (int i = GSIZE-1; i >= 0; i--)
           for (int j = 0; j < i; j++) {
              if (fabs(U[j]) < fabs(U[j+1])) {
                 float tmp = U[j];
                 U[j] = U[j+1];
                 U[j+1] = tmp;
                 string tmp2 = bacteria[j];
                 bacteria[j] = bacteria[j+1];
                 bacteria[j+1] = tmp2;
              }
           }
*/
        std::ofstream noafile(file.c_str(), std::ios::out);
        //noafile << "Name\tCentrality\tRank" << endl;
        noafile << "Name\tCentrality\tRank" << endl;
        float min = 0;
        float max = 0;
        for (int i = 0; i < GSIZE; i++) {
           /*U[i] = fabs(U[i]);
           if (fabs(U[i]) > max)
              max = fabs(U[i]);
           if (fabs(U[i]) < min)
              min = fabs(U[i]);*/
           //noafile << bacteria[i] << "\t" << U[i] << "\t\t" << GSIZE-i << endl;
	   if (U[i] != 0)
           noafile << bacteria[i] << "\t" << "#" << U[i] << " " << bacteria[i] << "\t" << U[i] << endl;
	   else
           noafile << bacteria[i] << "\t" << bacteria[i] << "\t" << "NR" << endl;
        }
        /*std::ofstream noafile(file.c_str(), std::ios::out);
        noafile << "Name\tCentrality\tRank" << endl;
        float min = 0;
        float max = 0;
        for (int i = 0; i < GSIZE; i++) {
           U[i] = fabs(U[i]);
           if (fabs(U[i]) > max)
              max = fabs(U[i]);
           if (fabs(U[i]) < min)
              min = fabs(U[i]);
           noafile << bacteria[i] << "\t" << U[i] << "\t\t" << GSIZE-i << endl;
        }*/

}




__global__ void _GPU_Copy_kernel(float* dst, float* src, unsigned long N) {
   //unsigned long threadNum = blockIdx.x*gridDim.y*blockDim.x + blockIdx.y*blockDim.x + threadIdx.x;
   unsigned long threadNum = blockIdx.x*blockDim.x + threadIdx.x;
   if (threadNum >= N) return;
   dst[threadNum] = src[threadNum];
}

//int row_pos = 2*(2*GSIZE*bac1 - bac1*(bac1-1));
//                                 int row_neg = row_pos + 2*(GSIZE-bac1);;
//
//                                 *pos_pos = row_pos + 2*(bac2-bac1);
//                                 *pos_neg = *pos_pos + 1;
//                                 *neg_pos = row_neg + 2*(bac2-bac1);
//                                 *neg_neg = *neg_pos + 1;

// If +/+ = 0, +/-=1, -/+=2 and -/- = 3
//
// k	max/i	max/j	i/j
// 0	0	0	0
// 1	0	1	1
// 2	1	0	2
// 3	1	1	3
// 4	2	2	0
// 5	2	3	1
// 6	3	2	2
// 7	3	3	3

__global__ void _GPU_Triad_kernel(float* G, int maxnode, int GSIZE) {
   int _8GSIZE = GSIZE*8;
   int N = _8GSIZE*GSIZE;
   int threadNum = blockIdx.x*blockDim.x + threadIdx.x;
   if (threadNum > N) return;
   int i = threadNum / _8GSIZE;
   int offset = threadNum - _8GSIZE*i;
   int j = offset / 8;
   if (i >= j || maxnode == i || maxnode == j) return;
   int k = offset % 8;
 
   int bac1 = maxnode;
   int bac2 = i;
   unsigned long maxnode_i, maxnode_j, i_j;
   if (maxnode > i) {
      bac1 = i;
      bac2 = maxnode;
   }
   maxnode_i = 2*(2*GSIZE*bac1 - bac1*(bac1-1)) + 2*(bac2-bac1)   + (k/4)*(2*(GSIZE-bac1)) + (k/2)%2;
   
   bac1 = maxnode;
   bac2 = j;
   if (maxnode > j) {
      bac1 = j;
      bac2 = maxnode;
   }

   maxnode_j = 2*(2*GSIZE*bac1 - bac1*(bac1-1)) + 2*(bac2-bac1) +  (k/4)*(2*(GSIZE-bac1)) + (k%2);

   // i < j here
   bac1 = i;
   bac2 = j;
   i_j = 2*(2*GSIZE*bac1 - bac1*(bac1-1)) + 2*(bac2-bac1)  + ((k/2)%2)*(2*(GSIZE-bac1)) + (k%2);
   
   float G_mi = G[maxnode_i];
   float G_mj = G[maxnode_j];
   float G_ij = G[i_j];

   bool edgeMaxI = (G_mi != 0);
   bool edgeMaxJ = (G_mj != 0);


   if (edgeMaxI && edgeMaxJ  &&
       G_ij != 0) {
      G[maxnode_i] = 2;
      G[maxnode_j] = 2; 
      G[i_j] = 2;
   }
   else {
      if (edgeMaxI) G[maxnode_i] = 2;
      if (edgeMaxJ) G[maxnode_j] = 2;
   }
}

__global__ void _GPU_Sweep_kernel(float* G, int N) {
   unsigned long threadNum = blockIdx.x*blockDim.x + threadIdx.x;
   if (threadNum > N) return;
   
   if (G[threadNum] == 2) G[threadNum] = 0;
}

__global__ void _GPU_Pay_kernel(float* D, float* P, const int* __restrict__ spots, int N, int GSIZE) {
    int node = blockIdx.x*blockDim.x+threadIdx.x;
    int _2node = node<<1;
    int _2GSIZE = GSIZE<<1;
    if (node >= N) return;
    float sum = -1; // Better
    for (int i = 0; i <= _2node; i++)
       sum += D[spots[i*_2GSIZE+_2node]];
    int loc = _2node*_2GSIZE + _2node+1;
    for (int i = _2node+1; i < _2GSIZE; i++) {
       sum += D[spots[loc++]];
    }
    P[node] = sum;
}

__global__ void _GPU_Spots_kernel(int* spots, int GSIZE) {
   int bac1 = blockIdx.y;
   int bac2 = blockIdx.x*blockDim.x + threadIdx.x;
   if (bac1 >= GSIZE*2 || bac2 >= GSIZE*2) return;
 
    int b1 = bac1 >> 1;
    int b2 = bac2 >> 1;
    if (b1 > b2) {
        int tmp = b1;
        b1 = b2;
        b2 = tmp;
    }
    spots[bac1*GSIZE*2+bac2] = (((GSIZE<<1)*b1 - b1*(b1-1))<<1) + (bac1%2)*((GSIZE-b1)<<1) + ((b2-b1)<<1) + (bac2%2);
  
}

// Way this works:
// Every thread BLOCK does one value of i (=blockIdx.y), and every THREAD does one value of j.
// Every function CALL does one value of k (the node going through).
__global__ void _GPU_Floyd_kernel(const int k, float *G, const int* __restrict__ spots,/*float* W,*/char* mark,const int N,const int first, int _2GSIZE){
        int i = blockIdx.y;
        int j = blockIdx.x*blockDim.x + threadIdx.x;
        int _2iGSIZE = i*_2GSIZE;
        int GSIZE = _2iGSIZE >> 1;
	unsigned long idx=spots[_2iGSIZE+j];/*getSpot(i, j,GSIZE);*/ // Spot for (i, j)
        
	if(i >= j || j>=N || j == k /*|| (first != 0 && mark[markid] == 0)*/)return; // If j is out of bounds, return

        float best = G[spots[_2iGSIZE+k]];
	if(best==0)return; // Because, we use best.
        float tmp_b=G[spots[k*_2GSIZE+j]];//G[getSpot(k, j)];  // Path from k to j
	if(tmp_b==0)return; // If this is zero, no path, return	
	float cur=best*tmp_b; // New path from i to j is the best path from i to k times the path from k to j
        // If the shortest path from i to j THROUGH this k is better
        // than all k's so far, replace
        if (fabs(cur) > fabs(G[idx]))
		G[idx]=cur;
}


PluginProxy<GPUATriaPlugin> ATriaPluginProxy = PluginProxy<GPUATriaPlugin>("GPUATria", PluginManager::getInstance());

