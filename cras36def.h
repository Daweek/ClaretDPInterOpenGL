#ifndef CRAS36DEF_H_
#define CRAS36DEF_H_
#pragma once
/*
		Definitions and Global variables.

    Edg@r J.
*/
// App related
#define VER 1.00
#define LAP_TIME
#define C_MASS
#define SUBWIN
#define CROSS
#define INFO

// Accelerator related
#define VTGRAPE // use Virtualized GRAPE library
//#define CPUONLY

// CUDA related
#define KER
//#define INTEROP
//#define DP
//#define TIME_MEMORY

// OpenGL Graphics related
#define GL_ON
#define STEREO 	0
#define TEXTURE 1   //Defined to load textures in case of using them.

#if defined(MDGRAPE3) || defined(VTGRAPE)
#define MDM 2      /* 0:host 2:m2 */
#else
#define MDM 0      /* 0:host 2:m2 */
#endif
#define SPC 0
#define ST2 0
#define TIP5P 1
#define SYS 0 /* 0:NaCl 1:water(fcc) 2:water(ice) 3:water(ice2) 4:NaCl-water */

// Memory related
#define S_NUM_MAX 10*10*10*8
#define W_NUM_MAX 10*10*10*8
#define ZERO_P 1
#define V_SCALE 0
#define T_CONST 1
#define P_CONST 0
#define KNUM 5                    /* number of particle type */
#define VMAX 462 /*1535*/        /* max value of wave nubmer vector */
#define EFT 12000
#define my_min(x,y) ((x)<(y) ? (x):(y))
#define my_max(x,y) ((x)>(y) ? (x):(y))

// Constant related
#define PI 		M_PI              /* pi */
#define PIT 	M_PI*2.0         /* 2 * pi */
#define PI2 	M_PI*M_PI        /* pi*pi */
#define IPI 	M_1_PI           /* 1/pi */
#define ISPI 	M_2_SQRTPI*0.5  /* 1 / sqrt(pi) */

#endif // Header
