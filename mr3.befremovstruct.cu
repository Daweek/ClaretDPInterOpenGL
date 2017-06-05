/*Archivo mr3.cu que contiene el codigo para CUDA.
    Renderizacion por OpenGL-CUDA interoperability
    Nucleo del codigo para calcular la fuerza entre particulas
    Creado por: Martinez Noriega Edgar Josafat
*/
#define GL_ON
#define KER
//#define DP
#define INTEROP
//#define TIME_MEMORY
/////////////
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <GL/glew.h>
#include <GL/freeglut.h>
// ***** CUDA includes
#include <cuda.h>
#include <nvcuvid.h>
#include <cudaGL.h>
#include <cuda_gl_interop.h>
#include <cuda_runtime_api.h>
#include <helper_cuda.h>

#define NMAX      8192
#define NTHRE       512
#define ATYPE        8
#define ATYPE2    (ATYPE * ATYPE)
#define ThreadsPB 512
//////For NaCl Optminized if_kernel
#define NTHREOPT      512
#define NDIVBIT      4
#define NDIV      (1<<NDIVBIT)
#define NTHREOPT2    (NTHREOPT/NDIV)

typedef struct {
  float r[3];
  int atype;
} VG_XVEC;

typedef struct {
  float pol;
  float sigm;
  float ipotro;
  float pc;
  float pd;
  float zz;
} VG_MATRIX;


/////////GLOBAL Variables/////////////////////////////////////////
int   	*d_atypemat;
VG_XVEC *d_x=NULL;
VG_XVEC	*vec=NULL;
float 	*d_force=NULL;
float   *d_side,*d_sideh;
float   *d_amass,*d_vl;
float 	*d_ekin1;
float 	*d_ekin,*d_xs,*d_mtemp,*d_mpres;
float		*d_poss,*d_colr;

int mem_flg=0;
int mem_flg2=0;
int mem_sp=5;
int mem_cpu=0;
int flg1=0,flg2=0,flg3=0;

extern GLuint g_possVBO, g_colorVBO;
extern cudaDeviceProp g_devprop;
extern struct cudaGraphicsResource* g_strucPossVBOCUDA;
extern struct cudaGraphicsResource* g_strucColorVBOCUDA;

__constant__
VG_MATRIX c_matrix[4]={[0].pol=1.250000,[0].sigm=2.340000,[0].ipotro=3.154574,[0].pc=0.072868,[0].pd=0.034699,[0].zz=1.000000,
	[1].pol=1.000000,[1].sigm=2.755000,[1].ipotro=3.154574,[1].pc=0.485784,[1].pd=0.602893,[1].zz=-1.000000,
	[2].pol=1.000000,[2].sigm=2.755000,[2].ipotro=3.154574,[2].pc=0.485784,[2].pd=0.602893,[2].zz=-1.000000,
	[3].pol=0.750000,[3].sigm=3.170000,[3].ipotro=3.154574,[3].pc=5.031334,[3].pd=10.106042,[3].zz=1.000000};

__constant__
float d_color_table[5][4]={ {0.35	,0.19	,0.19	,1.0},
														{0.19	,0.275,0.19	,1.0},
														{1.0	,0.4	,1.0	,1.0},
														{0.0	,0.8	,1.0	,1.0},
														{1.0	,1.0	,1.0	,1.0} };

//////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////FORCE CALCULATION WITH GPU/////////////////////////////////////
//////////////////////////////////////////////////////////////////////
__global__
void update_coor_kernel(int n3, float *vl,VG_XVEC *cd,float *xs,
                        float *fc,float *side){
#ifdef KER
	int tid  = threadIdx.x + blockIdx.x * blockDim.x;

	if (tid < n3){
		vl[tid]   =  (vl[tid]*(1-(*xs))+fc[tid])/(1+(*xs));
    cd[tid/3].r[tid % 3]   +=   vl[tid];
		if (cd[tid/3].r[tid % 3] < 0 || cd[tid/3].r[tid % 3] > side[tid % 3]) vl[tid] *= -1;
	}
#endif
}
//////////////////////////////////////////////////////////////////////////
__device__ __inline__
void inter_if(float xj[3], float xi[3], float fi[3], int t, float xmax,
		float xmax1) {
#ifdef KER

	int k;
	float dn2, r, inr, inr2, inr4, inr8, d3, dr[3];
	float pb = (float) (0.338e-19 / (14.39 * 1.60219e-19)), dphir;

	dn2 = 0.0f;
	for (k = 0; k < 3; k++) {
		dr[k] = xi[k] - xj[k];
		dr[k] -= rintf(dr[k] * xmax1) * xmax;
		dn2 += dr[k] * dr[k];
	}
	r = sqrtf(dn2);
#if 1
	inr = 1.0f / r;
#elif 0
	if(dn2 != 0.0f) inr = 1.0f / r;
	else inr = 0.0f;
#elif 0
	if(dn2 == 0.0f) inr = 0.0f;
	else inr = 1.0f / r;
#else
	inr = 1.0f / r;
	if(dn2 == 0.0f) inr = 0.0f;
#endif
	inr2 = inr * inr;
	inr4 = inr2 * inr2;
	inr8 = inr4 * inr4;
	d3 = pb * c_matrix[t].pol
			* expf((c_matrix[t].sigm - r) * c_matrix[t].ipotro);
	dphir =
			(d3 * c_matrix[t].ipotro * inr - 6.0f * c_matrix[t].pc * inr8
					- 8.0f * c_matrix[t].pd * inr8 * inr2
					+ inr2 * inr * c_matrix[t].zz);
#if 1
	if (dn2 == 0.0f)
		dphir = 0.0f;
#endif
	for (k = 0; k < 3; k++)
		fi[k] += dphir * dr[k];
#endif
}

__global__
void nacl_kernel_if2(VG_XVEC *x, int n, int nat, float xmax, float *fvec) {
#ifdef KER
	int tid = threadIdx.x;
	int jdiv = tid / NTHREOPT2;
	int i = blockIdx.x * NTHREOPT2 + (tid & (NTHREOPT2 - 1)); // Same + (tid %16)
	int j, k;
	float xmax1 = 1.0f / xmax;
	int atypei;
	float xi[3];
	__shared__ VG_XVEC s_xj[NTHREOPT];
	__shared__ float s_fi[NTHREOPT][3];

	for (k = 0; k < 3; k++)
		s_fi[tid][k] = 0.0f;
	for (k = 0; k < 3; k++)
		xi[k] = x[i].r[k];
	atypei = x[i].atype * nat;
	int na;
	na = n / NTHREOPT;
	na = na * NTHREOPT;
	for (j = 0; j < na; j += NTHREOPT) {
		__syncthreads();
		s_xj[tid] = x[j + tid];
		__syncthreads();
#pragma unroll 16
		for (int js = jdiv; js < NTHREOPT; js += NDIV)
			inter_if(s_xj[js].r, xi, s_fi[tid], atypei + s_xj[js].atype, xmax,
					xmax1);
	}
	for (j = na + jdiv; j < n; j += NDIV) {
		inter_if(x[j].r, xi, s_fi[tid], atypei + x[j].atype, xmax, xmax1);
	}
#if NTHREOPT>=512 && NTHREOPT2<=256
	__syncthreads();
	if(tid<256) for(k=0;k<3;k++) s_fi[tid][k]+=s_fi[tid+256][k];
#endif
#if NTHREOPT>=256 && NTHREOPT2<=128
	__syncthreads();
	if (tid < 128)
		for (k = 0; k < 3; k++)
			s_fi[tid][k] += s_fi[tid + 128][k];
#endif
#if NTHREOPT>=128 && NTHREOPT2<=64
	__syncthreads();
	if (tid < 64)
		for (k = 0; k < 3; k++)
			s_fi[tid][k] += s_fi[tid + 64][k];
#endif
#if NTHREOPT>=64 && NTHREOPT2<=32
	__syncthreads();
	if (tid < 32)
		for (k = 0; k < 3; k++)
			s_fi[tid][k] += s_fi[tid + 32][k];
#endif
#if NTHREOPT2<=16
	if (tid < 16)
		for (k = 0; k < 3; k++)
			s_fi[tid][k] += s_fi[tid + 16][k];
#endif
#if NTHREOPT2<=8
	if(tid<8) for(k=0;k<3;k++) s_fi[tid][k]+=s_fi[tid+8][k];
#endif
#if NTHREOPT2<=4
	if(tid<4) for(k=0;k<3;k++) s_fi[tid][k]+=s_fi[tid+4][k];
#endif
#if NTHREOPT2<=2
	if(tid<2) for(k=0;k<3;k++) s_fi[tid][k]+=s_fi[tid+2][k];
#endif
#if NTHREOPT2<=1
	if(tid<1) for(k=0;k<3;k++) s_fi[tid][k]+=s_fi[tid+1][k];
#endif
	if (jdiv == 0)
		for (k = 0; k < 3; k++)
			fvec[i * 3 + k] = s_fi[tid][k];
#endif
}

__global__
void velforce_kernel(int n3, float *fc, float *a_mass, float *vl,
                     VG_XVEC *atype, int *atype_mat, float hsq,float *ekin1,
                     float *poss, float *sideh){
#ifdef KER
	__shared__ float cache [ThreadsPB];
    int indx = threadIdx.x;
	int tid  = threadIdx.x + blockIdx.x * blockDim.x;

	cache [indx] = 0;

	if (tid < n3 ){
		fc[tid]-= fc[tid]/(n3/3);
		fc[tid] *= hsq/a_mass[atype_mat[atype[tid/3].atype]];
		cache [indx] = vl[tid]*vl[tid]*a_mass[atype_mat[atype[tid/3].atype]];
#ifdef INTEROP
		poss[tid] = atype[tid / 3].r[tid % 3]-sideh[tid % 3]; // for graphics VBO -- Position
#endif
	}
	__syncthreads();

	for (unsigned int s=blockDim.x/2; s>0; s>>=1)
	{
		if (indx < s)
		{
			cache[indx] += cache[indx + s];
		}
		__syncthreads();
	}
	if (indx == 0) ekin1[blockIdx.x] = cache [0];

#endif
}

__global__
void reduction (float *ekin,float *mtemp,float *mpres,float *xs,float tscale,
                    float nden, float vir,int s_num,int w_num,float rtemp,
					float lq,float hsq,float *ekin1, int limi){

#ifdef KER
	__shared__ float cache [NTHREOPT];

  int indx = threadIdx.x;

	cache [indx] = (indx < limi) ? ekin1[indx]:0.0f;

	__syncthreads();

	for (unsigned int s=NTHREOPT/2; s>0; s>>=1){
		if (indx < s)
		{
			cache[indx] += cache[indx + s];
		}
			__syncthreads();
	  }

	if (indx == 0){
		*ekin = cache [0];
		*ekin /= hsq;
		*mtemp = tscale * (*ekin);
		*mpres  = nden / 3.f * ((*ekin) - (vir)) / (s_num + w_num);
		*xs += (*mtemp - rtemp) /  lq * hsq *.5f;
	}

#endif
}

#ifdef INTEROP
__global__
void colorn4(int n4,float *vl,VG_XVEC *atype, int *atype_mat, float *colorvbo){
#ifdef KER
	int tid  = threadIdx.x + blockIdx.x * blockDim.x;
	float d0;
	float d0aux[4];

	d0 = (vl[tid/4]*vl[tid/4]+vl[tid/4+1]*vl[tid/4+1]+vl[tid/4+2]*vl[tid/4+2])*500;
	d0aux[0] 	= d0;
	d0aux[1] 	= d0/3;
	d0aux[2]	= d0/3;
	d0aux[3]	= 0;

	if (tid < n4){
		colorvbo[tid] = d_color_table[atype_mat[atype[tid/4].atype]][tid%4] + d0aux[tid%4];
	}
#endif
}
#endif


#ifdef DP
__global__
void md_loop_cuda (	int n3, float *vl,VG_XVEC *cd,float *xs,float *fc,float *side,
					int n, int nat, float xmax,
					float *a_mass, int *atype_mat, float hsq,float *ekin1,
					float *ekin,float *mtemp,float *mpres,float tscale,
					  float nden, float vir,int s_num,int w_num,float rtemp,
					  float lq,int limi,
					  int md_step, float *poss, float *sideh, float *colorvbo)
{
#if 1
	int  blocksPGrid = (n3 + ThreadsPB - 1)/(ThreadsPB);
	dim3 THREADS(NTHRE);
	dim3 BLOCKS((n3 + ThreadsPB - 1)/(ThreadsPB));
	dim3 threads(NTHREOPT);
	dim3 grid((n * NDIV + NTHREOPT - 1) / NTHREOPT);
	dim3 colorgridn4(((n*4) + ThreadsPB - 1)/(ThreadsPB));

	for(int md_loop = 0; md_loop < md_step; md_loop++){
		update_coor_kernel<<<BLOCKS,THREADS>>>(n3,vl,cd,xs,fc,side);
		nacl_kernel_if2<<<grid, threads>>>(cd, n, nat, xmax, fc);
		velforce_kernel<<<BLOCKS,THREADS>>>(n3,fc,a_mass,vl,cd,atype_mat,hsq,ekin1,poss,sideh);
		reduction<<<1,NTHRE>>>(ekin,mtemp,mpres,xs,tscale,nden,vir,s_num,w_num,rtemp,lq,hsq,ekin1,blocksPGrid);
	}
#ifdef INTEROP
	colorn4<<<colorgridn4,THREADS>>>(n*4,vl,cd,atype_mat,colorvbo);
#endif

#endif
}

#endif
//////////////////NaCl Optmized
///////////////////////////////

extern "C"
void mdlop(int n3,int grape_flg,double phi [3],double *phir,double *iphi, double *vir,int s_num3,
			timeval time_v,double *md_time0,double *md_time,int *m_clock,int md_step,double *mtemp,
			double tscale,double *mpres,double nden,int s_num,int w_num,double rtemp,double lq,
			double x[], int n, int atype[], int nat,
			double pol[], double sigm[], double ipotro[],
		 	double pc[], double pd[],double zz[],
		 	int tblno, double xmax, int periodicflag,
		 	double force[],
			double hsq,double a_mass [], int atype_mat [], double *ekin,double *vl,
			double *xs,double side [],int *firstmalloc, double sideh[]){

//////////////VARIABLES FROM THE BEGINING/////////////////
	//int md_loop;
	//float *forcef=NULL;
	int i,j;
	float xmaxf;
  if((periodicflag & 1)==0) xmax*=2.0;
	xmaxf=xmax;
  int n4 = n*4;
/////////////////////////////////////////////////////////
	int  blocksPGrid = (n3 + ThreadsPB - 1)/(ThreadsPB);
	dim3 THREADS(NTHRE);
	dim3 BLOCKS((n3 + ThreadsPB - 1)/(ThreadsPB));
	dim3 threads(NTHREOPT);
	dim3 grid((n * NDIV + NTHREOPT - 1) / NTHREOPT);
	dim3 colorgridn4((n4 + ThreadsPB - 1)/(ThreadsPB));

	float   fxs = *xs;
	float   fside[3],*ffc, fsideh[3];
	float   *vla;
	VG_XVEC	*veca;

	int     p = 0;
	float   hsqf = hsq;
	float   *fvl,fa_mass[4];

	float ftscale = tscale,fnden = nden,frtemp = rtemp,flq = lq,fvir = 0;
	float fmtemp = *mtemp,fmpres = *mpres;

	vla		= (float*)	malloc(n3*sizeof(float));
	veca  = (VG_XVEC*)malloc((n+NTHREOPT2)*sizeof(VG_XVEC));


	if(*firstmalloc == 0){

		printf("CUDA malloc time...\n");

		// Allocating memory for float conversion.
		ffc = (float*)		malloc(n3*sizeof(float));
		fvl = (float*)		malloc(n3*sizeof(float));
		vec = (VG_XVEC*) 	malloc((NMAX+NTHREOPT2)*sizeof(VG_XVEC));

		// Conversion from Double to Float
		for (p=0;p<4;p++) fa_mass[p] = (float) a_mass[p];
		for (p=0;p<3;p++) fside[p] 	 = (float) side[p];
		for (p=0;p<3;p++) fsideh[p]  = (float) sideh[p];
		for (p=0;p<n3;p++){
			fvl     [p] =  (float) *(vl +p);
			ffc     [p] =  (float) *(force +p);
		}

		for (i = 0; i < (n + NTHREOPT2 - 1) / NTHREOPT2 * NTHREOPT2; i++) {
			if (i < n) {
				for (j = 0; j < 3; j++) {
					vec[i].r[j] = x[i * 3 + j];
				}
				vec[i].atype = atype[i];
			}
			else {
				for (j = 0; j < 3; j++) {
					vec[i].r[j] = 0.0f;
				}
				vec[i].atype = 0;
			}
		}

		// Free CUDA memory. In case we already allocate
		checkCudaErrors(cudaFree(d_x));
		checkCudaErrors(cudaFree(d_force));
		checkCudaErrors(cudaFree(d_side));
		checkCudaErrors(cudaFree(d_sideh));
		checkCudaErrors(cudaFree(d_amass));
		checkCudaErrors(cudaFree(d_vl));
		checkCudaErrors(cudaFree(d_atypemat));
		checkCudaErrors(cudaFree(d_ekin));
		checkCudaErrors(cudaFree(d_xs));
		checkCudaErrors(cudaFree(d_mtemp));
		checkCudaErrors(cudaFree(d_mpres));
		checkCudaErrors(cudaFree(d_ekin1));


		// Allocate global memory to GPU
		checkCudaErrors(cudaMalloc((void**)&d_x,sizeof(VG_XVEC)* (NMAX + NTHREOPT2)));
		checkCudaErrors(cudaMalloc((void**)&d_force,sizeof(float)*(NMAX + NTHREOPT2)*3));
		checkCudaErrors(cudaMalloc((void**)&d_side,3*sizeof(float)));
		checkCudaErrors(cudaMalloc((void**)&d_sideh,3*sizeof(float)));
		checkCudaErrors(cudaMalloc((void**)&d_amass,4*sizeof(float)));
		checkCudaErrors(cudaMalloc((void**)&d_vl,n3*sizeof(float)));
		checkCudaErrors(cudaMalloc((void**)&d_atypemat,20*sizeof(int)));
		checkCudaErrors(cudaMalloc((void**)&d_ekin,sizeof(float)));
		checkCudaErrors(cudaMalloc((void**)&d_xs,sizeof(float)));
		checkCudaErrors(cudaMalloc((void**)&d_mtemp,sizeof(float)));
		checkCudaErrors(cudaMalloc((void**)&d_mpres,sizeof(float)));
		checkCudaErrors(cudaMalloc((void**)&d_ekin1,blocksPGrid*sizeof(float)));
		//checkCudaErrors(cudaMalloc((void**)&d_poss,n3*sizeof(float)));

		// Copy memory from CPU to GPU
		checkCudaErrors(cudaMemcpy(d_x,vec,sizeof(VG_XVEC)*((n + NTHREOPT2 - 1) / NTHREOPT2 * NTHREOPT2),cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_side,fside,sizeof(float)*3,cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_sideh,fsideh,sizeof(float)*3,cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_mtemp,&fmtemp,sizeof(float),cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_mpres,&fmpres,sizeof(float),cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_xs,&fxs,sizeof(float),cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_vl,fvl,sizeof(float)*n3,cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_amass,fa_mass,sizeof(float)*4,cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_atypemat,atype_mat,sizeof(int)*20,cudaMemcpyHostToDevice));
		checkCudaErrors(cudaMemcpy(d_force,ffc,sizeof(float)*n*3,cudaMemcpyHostToDevice));

		// Free the memory used to convert to Float
		free(ffc);
		free(fvl);
		free(vec);
	}

#ifdef INTEROP
	//Interoperability
	// Position
	size_t vbosizepos;
	checkCudaErrors(cudaGraphicsMapResources(1,&g_strucPossVBOCUDA,0));
	checkCudaErrors(cudaGraphicsResourceGetMappedPointer((void**)&d_poss,
																												&vbosizepos,
																												g_strucPossVBOCUDA));
	// Color
	size_t vbosizecol;
	checkCudaErrors(cudaGraphicsMapResources(1,&g_strucColorVBOCUDA,0));
	checkCudaErrors(cudaGraphicsResourceGetMappedPointer((void**)&d_colr,
																												&vbosizecol,
																												g_strucColorVBOCUDA));
#endif

	///////Md_loop///////////////////////////////////////////////

#ifdef DP
#ifndef TIME_MEMORY
	gettimeofday(&time_v,NULL);
	*md_time0 = (time_v.tv_sec + time_v.tv_usec / 1000000.0);
#endif
//	for (int m=0;m<1000;m++){
	md_loop_cuda<<<1,1>>>(n3,d_vl,d_x,d_xs,d_force,d_side,
						n,nat,xmaxf,
						d_amass,d_atypemat,hsqf,d_ekin1,
						d_ekin,d_mtemp,d_mpres,ftscale,fnden,fvir,s_num,w_num,frtemp,flq,blocksPGrid,
						md_step,d_poss,d_sideh,d_colr);
	//}
	*m_clock+=md_step;
	cudaDeviceSynchronize();

#ifndef TIME_MEMORY
	gettimeofday(&time_v,NULL);
	*md_time = (time_v.tv_sec + time_v.tv_usec / 1000000.0);
#endif

#else
	gettimeofday(&time_v,NULL);
	*md_time0 = (time_v.tv_sec + time_v.tv_usec / 1000000.0);
//	for (int m=0;m<10000;m++){
	for(int md_loop = 0; md_loop < md_step; md_loop++){
		update_coor_kernel<<<BLOCKS,THREADS>>>(n3,d_vl,d_x,d_xs,d_force,d_side);
		nacl_kernel_if2<<<grid, threads>>>(d_x, n, nat, xmaxf, d_force);
		velforce_kernel<<<BLOCKS,THREADS>>>(n3,d_force,d_amass,d_vl,d_x,d_atypemat,hsqf,d_ekin1,d_poss,d_sideh);
		reduction<<<1,threads>>>(d_ekin,d_mtemp,d_mpres,d_xs,ftscale,fnden,fvir,s_num,w_num,frtemp,flq,hsqf,d_ekin1,blocksPGrid);
	}
#ifdef INTEROP
	colorn4<<<colorgridn4,THREADS>>>(n4,d_vl,d_x,d_atypemat,d_colr); // Just update after the cycle. For color output.
#endif
//	}
	*m_clock+=md_step;
	cudaDeviceSynchronize();
	gettimeofday(&time_v,NULL);
	*md_time = (time_v.tv_sec + time_v.tv_usec / 1000000.0);
#endif

/////////////////Copy back to the CPU
	//CUDA_SAFE_CALL(cudaMemcpy(forcef,d_force,sizeof(float)*n*3,cudaMemcpyDeviceToHost));
	//CUDA_SAFE_CALL(cudaMemcpy(&fxs,d_xs,sizeof(float),cudaMemcpyDeviceToHost));
	//CUDA_SAFE_CALL(cudaMemcpy(&ekinaux,d_ekin,sizeof(float),cudaMemcpyDeviceToHost));
	//CUDA_SAFE_CALL(cudaMemcpy(&fmtemp,d_mtemp,sizeof(float),cudaMemcpyDeviceToHost));
	//CUDA_SAFE_CALL(cudaMemcpy(&fmpres,d_mpres,sizeof(float),cudaMemcpyDeviceToHost));

#ifdef TIME_MEMORY
	gettimeofday(&time_v,NULL);
	*md_time0 = (time_v.tv_sec + time_v.tv_usec / 1000000.0);
#endif

#ifdef INTEROP
	checkCudaErrors(cudaGraphicsUnmapResources(1,&g_strucPossVBOCUDA,0));
	checkCudaErrors(cudaGraphicsUnmapResources(1,&g_strucColorVBOCUDA,0));
#endif

#ifndef INTEROP
	checkCudaErrors(cudaMemcpy(vla,d_vl,n3*sizeof(float),cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy(veca,d_x,n*sizeof(VG_XVEC),cudaMemcpyDeviceToHost));
	//checkCudaErrors(cudaMemcpy(cord,d_poss,n3*sizeof(float),cudaMemcpyDeviceToHost));
#endif

#ifdef TIME_MEMORY
	gettimeofday(&time_v,NULL);
	*md_time = (time_v.tv_sec + time_v.tv_usec / 1000000.0);
#endif


	//for(i=0;i<n;i++) for(j=0;j<3;j++) force[i*3+j]=(double) forcef[i*3+j];
	for(p=0;p<n3;p++) *(vl+p) = (double) vla[p];
	for(i=0;i<n;i++)for(j=0;j<3;j++) *(x+i*3+j) = (double)veca[i].r[j];
	//for(i=0;i<n;i++)for(j=0;j<3;j++) *(x+i*3+j) = (double)cord[j+i*3];


	free(veca);
	free(vla);
	//free(cord);
	*firstmalloc = 1;

}



