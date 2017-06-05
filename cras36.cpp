#include "cras36.h" //Include fisrt this one


void initGL(int argc, char **argv)
{

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(winWidth, winHeight);
  glutInitWindowPosition(900,500);
  WindowID = glutCreateWindow("Claret: MD simulator");

  glutDisplayFunc(display);
  glutReshapeFunc(reshape);
  glutMouseFunc(mouse);
  glutMotionFunc(motion);
  glutKeyboardFunc(keyboard);

	//	Glew initialization to have access to the extensions
	if (glewInit() != GLEW_OK)
		assert(!"Failed to initialize GLEW.\n");
	if (!glewIsSupported("GL_EXT_framebuffer_object"))
		assert(!"The GL_EXT_framebuffer_object extension is required.\n");

}

void initCUDA(){
	int 	devCount = 0;
	char 	devName [256];
	int 	n = n3/3; //Number of particles
	int 	errorcode;

#if 0
	checkCudaErrors(cuInit(0));
	checkCudaErrors(cuDeviceGetCount(&devCount));

	for(int i = 0; i < devCount; i++)
	{
		checkCudaErrors(cuDeviceGet(&g_device, i));
		checkCudaErrors(cuDeviceGetName(devName, 256, g_device));
	}
#endif

#ifndef CPUONLY
	checkCudaErrors(cudaSetDevice(0));
	checkCudaErrors(cudaGetDevice(&g_devgpu));
	checkCudaErrors(cudaGetDeviceProperties(&g_devprop,g_devgpu));

#ifdef INTEROP
	checkCudaErrors(cudaGLSetGLDevice(g_devgpu));
#endif
#endif

	printf("Device %d: %s is used!\n", g_devgpu, g_devprop.name);

}


int main(int argc, char **argv)
{
  int i,j,k;
  int i0,i1;
  double d0;
  char sbuf[50];
  char tt_name[256];
  int md_loop;

  ////Default Configuration (No arguments Passed to Claret)///
  np 					= 11;   // Number of particles from 1 - 12
  temp 				= 300;	// Initial Temperature
  md_step			=	10;		// Inner loop which calls to compute force
  md_loop			=	1; 		// For no OpenGL loop
  grape_flg 	= 1;		// CPU =0 , GPU = 1
  rendermode	= 3;		// 1=points, 2=Textures, 3=Spheres

  printf("MD_STEP=%d\tMD_LOOP=%d\tPARTICLE=%d\tACCEL=%d\tRENDERMODE=%d\n"
  		  				 ,md_step,md_loop,np,grape_flg,rendermode);
  ////////////////////////////

  /////////////////////////////////////////////
  //Reading Arguments
  if ( 1 <= argc){

  	if(argc == 2){
		  md_step = atoi(argv[1]);
		  printf("MD_STEP=%d\tMD_LOOP=%d\tPARTICLE=%d\tACCEL=%d\tRENDERMODE=%d\n"
		  				 ,md_step,md_loop,np,grape_flg,rendermode);
	  }

	  if(argc == 3){
		  md_step = atoi(argv[1]);
		  md_loop = atoi(argv[2]);
		  printf("MD_STEP=%d\tMD_LOOP=%d\tPARTICLE=%d\tACCEL=%d\tRENDERMODE=%d\n"
		  				 ,md_step,md_loop,np,grape_flg,rendermode);
	  }

	  if(argc == 4){
		  md_step = atoi(argv[1]);
		  md_loop = atoi(argv[2]);
		  np	  	= atoi(argv[3]);
		  printf("MD_STEP=%d\tMD_LOOP=%d\tPARTICLE=%d\tACCEL=%d\tRENDERMODE=%d\n"
		  				 ,md_step,md_loop,np,grape_flg,rendermode);
	  }

	  if(argc == 5){
		  md_step = atoi(argv[1]);
		  md_loop = atoi(argv[2]);
		  np	  	= atoi(argv[3]);
		  grape_flg	  = atoi(argv[4]);
		  printf("MD_STEP=%d\tMD_LOOP=%d\tPARTICLE=%d\tACCEL=%d\tRENDERMODE=%d\n"
		  				 ,md_step,md_loop,np,grape_flg,rendermode);
	  }

	  if(argc == 6){
		  md_step = atoi(argv[1]);
		  md_loop = atoi(argv[2]);
		  np	  	= atoi(argv[3]);
		  grape_flg	  = atoi(argv[4]);
		  rendermode  = atoi(argv[5]);
		  printf("MD_STEP=%d\tMD_LOOP=%d\tPARTICLE=%d\tACCEL=%d\tRENDERMODE=%d\n"
		  				 ,md_step,md_loop,np,grape_flg,rendermode);
	  }
  }

/////////////////////////////////////////////

#ifdef GL_ON
	sub_x = 1.5;
	sub_y = 1.5;
	temp_ymax = 2000;

 	initGL(argc, argv);
  init();
  initCUDA();

#endif

  init_MD();
  keep_mem(S_NUM_MAX,W_NUM_MAX*w_site);
  set_cd(1);

////Drawing With OpenGL///////////////////////
#ifdef GL_ON

  while(g_glutLoopContinue){

  	if(run_flg == 1)
  			md_run();

  	glutMainLoopEvent();

  }
#else
  //////////Using only results from Console///////////////
  char *acc = "CPU";


  if (grape_flg == 0) acc ="CPU";
  else acc = "GPU";
  printf("Starting MD for NaCl\n");
  printf("\nAccelerator type:%s\n",acc);
  printf("Number of Particles:%d\n",n1);
  printf("MD_LOOP:%d cicles\n",md_loop);
  printf("MD_step:%d steps\n\n\n",md_step);



  for (i=0;i<md_loop;i++){
	  gettimeofday(&time_v,NULL);
	  disp_time0 = (time_v.tv_sec + time_v.tv_usec / 1000000.0);

	  md_run();

	  gettimeofday(&time_v,NULL);
	  disp_time = (time_v.tv_sec + time_v.tv_usec / 1000000.0);
	  printf("Time to complete one full cicle of MD_LOOP: %3fsec\n\n",disp_time-disp_time0);
  }

#endif

  return 0;
}

