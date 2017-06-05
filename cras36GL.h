///////////////////////////For Test in ClaretEncDec
//////////////////OpenGL/////////////////////
//////////////////Functions//////////////////
//////////////////////////////////////////////
#ifdef NOMAINCLARET
extern bool	g_glutLoopContinue;
#endif
//#ifdef GL_ON
void CircleTable (float **sint,float **cost,const int n){
    int i;
    /* Table size, the sign of n flips the circle direction */
    const int size = abs(n);
    /* Determine the angle between samples */
    const float angle = 2*M_PI/(float)( ( n == 0 ) ? 1 : n );
    /* Allocate memory for n samples, plus duplicate of first entry at the end */
    *sint = (float *) calloc(sizeof(float), size+1);
    *cost = (float *) calloc(sizeof(float), size+1);
    /* Bail out if memory allocation fails, fgError never returns */
    if (!(*sint) || !(*cost))
    {
        free(*sint);
        free(*cost);
        printf("Failed to allocate memory in fghCircleTable");
		exit(0);
    }
    /* Compute cos and sin around the circle */
    (*sint)[0] = 0.0;
    (*cost)[0] = 1.0;
    for (i=1; i<size; i++)
    {
        (*sint)[i] = sin(angle*i);
        (*cost)[i] = cos(angle*i);
    }
    /* Last sample is duplicate of the first */
    (*sint)[size] = (*sint)[0];
    (*cost)[size] = (*cost)[0];
}

void mat_inv(double a[4][4])
{
  int i,j,k;
  double t, u, det;
  int n = 3;

  det = 1;
  for(k = 0; k < n; k++){
    t = a[k][k]; det *= t;
    for(i = 0; i < n; i++) a[k][i] /= t;
    a[k][k] = 1 / t;
    for(j = 0; j < n; j++)
      if(j != k){
        u = a[j][k];
        for(i = 0; i < n; i++)
          if(i != k) a[j][i] -= a[k][i] * u;
          else       a[j][i] = -u/t;
      }
  }
}
//#if TEXTURE == 1
void TexSphere(double r, int num)
{
  int i;

  glEnable(GL_TEXTURE_2D);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
#ifdef GL_VERSION_1_1
   glBindTexture(GL_TEXTURE_2D, sp_tex[num]);
#endif

  glBegin(GL_POLYGON);
  glNormal3d(1,0,0);
  for(i = 0; i < CIRCLE; i++){
    glTexCoord2f(circle_cd[i][1]*.5+.5, circle_cd[i][2]*.5+.5);
    glVertex3d(circle_cd[i][0],r*circle_cd[i][1],r*circle_cd[i][2]);
  }
  glEnd();
  glDisable(GL_TEXTURE_2D);
}
void TexCirSphere(double r, int num)
{

  glEnable(GL_TEXTURE_2D);
  glEnable(GL_ALPHA_TEST);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  //glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
#ifdef GL_VERSION_1_1
   glBindTexture(GL_TEXTURE_2D, sp_tex[num]);
#endif

  glBegin(GL_QUADS);
  glNormal3f(1.0, 0.0, 0.0);
  glTexCoord2f(0.0, 0.0); glVertex3f(0.0, -r, -r);
  glTexCoord2f(0.0, 1.0); glVertex3f(0.0,  r, -r);
  glTexCoord2f(1.0, 1.0); glVertex3f(0.0,  r,  r);
  glTexCoord2f(1.0, 0.0); glVertex3f(0.0, -r,  r);

  glEnd();
  glDisable(GL_ALPHA_TEST);
  glDisable(GL_TEXTURE_2D);
}
//#endif

void readtexture(char *filename)
{
  int i,j,k;
  char buf[256];
  char *buf2;
  FILE *fp;

  if((fp = fopen(filename,"r")) == NULL){
/*
    printf("texture file open error %s\n",filename);
    exit(1);
*/
    texf_flg = 0;
    tex_flg = 0;
  }

  if(texf_flg == 1){
   buf2 = fgets(buf,256,fp);
   buf2 = fgets(buf,256,fp);
   buf2 = fgets(buf,256,fp);

    for(i = 0; i < X_PIXEL; i++){
      for(j = 0; j < Y_PIXEL; j++){
        for(k = 0; k < 3; k++){
          teximage[j][Y_PIXEL-1-i][k==2 ? 0:k+1] = fgetc(fp);
        }
        if(teximage[j][Y_PIXEL-1-i][0] == 0 &&
           teximage[j][Y_PIXEL-1-i][1] == 0 &&
           teximage[j][Y_PIXEL-1-i][2] == 0){
          teximage[j][Y_PIXEL-1-i][3] = (GLubyte)0;
        } else {
          teximage[j][Y_PIXEL-1-i][3] = (GLubyte)255;
        }
      }
    }
    fclose(fp);
  }
}


void init(void)
{
  int i,j;
  int i0;
  int a_num = 1;
#if 0 // Can not be used under X forwarding
  GLfloat mat_specular[] = {0.2, 0.2, 0.2, 1.0};
  GLfloat mat_ambient[] = {0.1, 0.1, 0.1, 1.0};
  GLfloat mat_shininess[] = {64.0};
  GLfloat light_position[] = {1.0, 1.1, 1.2, 0.0};

  glShadeModel(GL_SMOOTH);
/*  glShadeModel(GL_FLAT);*/
  glLightfv(GL_LIGHT0, GL_SPECULAR, mat_specular);
  glLightfv(GL_LIGHT0, GL_SHININESS, mat_shininess);
  glLightfv(GL_LIGHT0, GL_AMBIENT, mat_ambient);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
#endif
  glMatrixMode(GL_MODELVIEW);
  glGetDoublev(GL_MODELVIEW_MATRIX,m_matrix);
  glGetDoublev(GL_MODELVIEW_MATRIX,i_matrix);

  base = glGenLists(128);
  for(i = 0; i < 128; i++){
    glNewList(base+i, GL_COMPILE);
    glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, i);
    glEndList();
  }
  glListBase(base);

  color_table[0][0] = 0.7;
  color_table[0][1] = 0.38;
  color_table[0][2] = 0.38;
  color_table[0][3] = 1;

  color_table[1][0] = 0.38;
  color_table[1][1] = 0.55;
  color_table[1][2] = 0.38;
  color_table[1][3] = 1;

  for(i = 0; i < 3; i++){
    color_table[0][i] /= 2.0;
    color_table[1][i] /= 2.0;
  }
/*
  printf("%f %f %f\n",color_table[0][0],color_table[0][1],color_table[0][2]);
  printf("%f %f %f\n",color_table[1][0],color_table[1][1],color_table[1][2]);
*/
  color_table[2][0] = 1;
  color_table[2][1] = .4;
  color_table[2][2] = 1;
  color_table[2][3] = 1;

  color_table[3][0] = 0;
  color_table[3][1] = 0.8;
  color_table[3][2] = 1;
  color_table[3][3] = 1;

  color_table[4][0] = 1;
  color_table[4][1] = 1;
  color_table[4][2] = 1;
  color_table[4][3] = 1;

  r_table[0] = 2.443/2;
  r_table[1] = 3.487/2;
  r_table[2] = 3.156/2;
  r_table[3] = .7;
  r_table[4] = .7;

//#if TEXTURE == 1
  char file_name[4][30];
  strcpy(file_name[0],"SP_00.PPM");
  strcpy(file_name[1],"SP_01.PPM");

  glGenTextures(2, sp_tex);

  for(i = 0; i < 2; i++){
     readtexture(file_name[i]);
     glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

     glBindTexture(GL_TEXTURE_2D, sp_tex[i]);

     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

     glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, Y_PIXEL, X_PIXEL,
                  0, GL_RGBA, GL_UNSIGNED_BYTE, teximage);


   }
  glAlphaFunc(GL_GREATER,0.5);

//#endif


  if(kflg == 1){
    if( ( fp = fopen( k_file, "w" ) ) == NULL ){
      printf("DATA file open error\n");
      exit( 1 );
    }
  }
#ifdef INFO
  for(i = 0; i < 3; i++)
    trans0[i] = 0;
  for(i = 0; i < 4; i++){
    for(j = 0; j < 4; j++){
      if(i == j){
	matrix0[i*4+j] = 1;
      } else {
	matrix0[i*4+j] = 0;
      }
    }
  }
#endif
}
void hako(int flg)
{
  double d0;
  int i;
  static GLfloat kabe[]  = { 0.0, 0.0, 0.4, 1.0 };
  static GLfloat kabe2[]  = { 0.0, 0.0, 0.8, 1.0 };
  double side_s[3],side_e[3];

  for(i = 0; i < 3; i++){
    side_s[i] = -sideh[i];
    side_e[i] = side[i]-sideh[i];
  }

  if(flg == 0){
    for(i = 0; i < 3; i++){
      side_s[i] += -radius;
      side_e[i] +=  radius;
    }
    /*    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, kabe);*/
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, kabe);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, kabe2);

    glBegin(GL_POLYGON);
    glNormal3d(1,0,0);
    glVertex3d(side_s[0],side_s[1],side_s[2]);
    glVertex3d(side_s[0],side_e[1],side_s[2]);
    glVertex3d(side_s[0],side_e[1],side_e[2]);
    glVertex3d(side_s[0],side_s[1],side_e[2]);
    glEnd();
    glBegin(GL_POLYGON);
    glNormal3d(-1,0,0);
    glVertex3d(side_e[0],side_s[1],side_s[2]);
    glVertex3d(side_e[0],side_s[1],side_e[2]);
    glVertex3d(side_e[0],side_e[1],side_e[2]);
    glVertex3d(side_e[0],side_e[1],side_s[2]);
    glEnd();

    glBegin(GL_POLYGON);
    glNormal3d(0,-1,0);
    glVertex3d(side_s[0],side_e[1],side_s[2]);
    glVertex3d(side_e[0],side_e[1],side_s[2]);
    glVertex3d(side_e[0],side_e[1],side_e[2]);
    glVertex3d(side_s[0],side_e[1],side_e[2]);
    glEnd();
    glBegin(GL_POLYGON);
    glNormal3d(0,1,0);
    glVertex3d(side_s[0],side_s[1],side_s[2]);
    glVertex3d(side_s[0],side_s[1],side_e[2]);
    glVertex3d(side_e[0],side_s[1],side_e[2]);
    glVertex3d(side_e[0],side_s[1],side_s[2]);
    glEnd();

    glBegin(GL_POLYGON);
    glNormal3d(0,0,1);
    glVertex3d(side_s[0],side_s[1],side_s[2]);
    glVertex3d(side_e[0],side_s[1],side_s[2]);
    glVertex3d(side_e[0],side_e[1],side_s[2]);
    glVertex3d(side_s[0],side_e[1],side_s[2]);
    glEnd();

    glBegin(GL_POLYGON);
    glNormal3d(0,0,-1);
    glVertex3d(side_s[0],side_s[1],side_e[2]);
    glVertex3d(side_s[0],side_e[1],side_e[2]);
    glVertex3d(side_e[0],side_e[1],side_e[2]);
    glVertex3d(side_e[0],side_s[1],side_e[2]);
    glEnd();
  }
  if(flg == 1){
    glColor3d(1.0,1.0,1.0);
    glBegin(GL_LINE_LOOP);
    glVertex3d(side_s[0],side_s[1],side_s[2]);
    glVertex3d(side_s[0],side_e[1],side_s[2]);
    glVertex3d(side_s[0],side_e[1],side_e[2]);
    glVertex3d(side_s[0],side_s[1],side_e[2]);
    glEnd();
    glBegin(GL_LINE_LOOP);
    glVertex3d(side_e[0],side_s[1],side_s[2]);
    glVertex3d(side_e[0],side_e[1],side_s[2]);
    glVertex3d(side_e[0],side_e[1],side_e[2]);
    glVertex3d(side_e[0],side_s[1],side_e[2]);
    glEnd();
    glBegin(GL_LINES);
    glVertex3d(side_e[0],side_s[1],side_s[2]);
    glVertex3d(side_s[0],side_s[1],side_s[2]);
    glVertex3d(side_e[0],side_e[1],side_s[2]);
    glVertex3d(side_s[0],side_e[1],side_s[2]);
    glVertex3d(side_e[0],side_e[1],side_e[2]);
    glVertex3d(side_s[0],side_e[1],side_e[2]);
    glVertex3d(side_e[0],side_s[1],side_e[2]);
    glVertex3d(side_s[0],side_s[1],side_e[2]);
    glEnd();
  }
}
void bou2(double x0, double y0, double z0, double x1, double y1, double z1
   ,double wid,int dit)
{
  double d0,d2;
  GLUquadricObj *qobj;

  qobj = gluNewQuadric();
  gluQuadricDrawStyle(qobj,GLU_FILL);
  gluQuadricNormals(qobj,GLU_SMOOTH);

  d0 = sqrt((x0-x1)*(x0-x1) + (y0-y1)*(y0-y1) + (z0-z1)*(z0-z1));

  glPushMatrix();
  glTranslated(x0,y0,z0);
  d2 =-acos((z1-z0)/d0)/M_PI*180;
  /*
  printf("%f %f %f  %f %f %f  %f %f\n",x0,y0,z0,x1,y1,z1,d0,d2);
  */
  if(y0 == y1 && x0 == x1)
    glRotatef(d2,1,0,0);
  else
    glRotatef(d2,(y1-y0),-(x1-x0),0);
  gluCylinder(qobj,wid,wid,d0,dit,1);
  glPopMatrix();

  glPushMatrix();
  glTranslated(x1,y1,z1);
  glPopMatrix();

  gluDeleteQuadric(qobj);
}
void bond_drow(int s_num, int e_num, int c_num, double wid,int dit)
{
  glMaterialfv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,bond_color);
  bou2(cd[s_num*3]-sideh[0],cd[s_num*3+1]-sideh[1],cd[s_num*3+2]-sideh[2],
       cd[e_num*3]-sideh[0],cd[e_num*3+1]-sideh[1],cd[e_num*3+2]-sideh[2]
       ,wid,dit);
}
void line_drow(int s_num, int e_num)
{
  int i;
  GLfloat color[4];

  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE,bond_color);
  glLineWidth(1.0);
  glBegin(GL_LINES);
  glVertex3d(cd[s_num*3]-sideh[0],cd[s_num*3+1]-sideh[1],cd[s_num*3+2]-sideh[2]);
  glVertex3d(cd[e_num*3]-sideh[0],cd[e_num*3+1]-sideh[1],cd[e_num*3+2]-sideh[2]);
  glEnd();

}
void small_font(double px, double py, double pz,const char *moji)
{
  int i;
  int len;
  double wid,adj;

  wid = 0.1;
  len = strlen(moji);
  glColor4fv(moji_c[(int)(clear_color+.5)]);
  for(i = 0;i < len; i++){
    if(moji[i] == '1') adj = 0.55;
    else if(moji[i] >= '2' && moji[i] <= '9') adj = 0.7;
    else if(moji[i] == '0') adj = 0.7;
    else if(moji[i] == 'B') adj = 0.7;
    else adj = 1.0;
    glRasterPos3d(px,py,pz);
    px += wid*adj;
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_10, moji[i]);
  }
}
void medium_font(double px, double py, double pz,const char *moji)
{
	int i;
  int len;
  double wid,adj;

  wid = 0.1;
  len = strlen(moji);
  glColor4fv(moji_c[(int)(clear_color+.5)]);
  for(i = 0;i < len; i++){
    if(moji[i] == '1') adj = 0.55;
    else if(moji[i] >= '2' && moji[i] <= '9') adj = 0.7;
    else if(moji[i] == '0') adj = 0.7;
    else if(moji[i] == 'B') adj = 0.7;
    else adj = 1.0;
    glRasterPos3d(px,py,pz);
    px += wid*adj;
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, moji[i]);
  }
}

void reshape(int w, int h)
{

#if 0
  glViewport(0, 0, (GLsizei)w, (GLsizei)h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0, (double)w / (double)h, 0.1, 1600.0);
  glMatrixMode(GL_MODELVIEW);
#endif
}
void mouse(int button, int state, int x, int y)
{
  switch (button) {
  case GLUT_LEFT_BUTTON:
    if (state == GLUT_DOWN) {
      mpos[0] = x;
      mpos[1] = y;
      mouse_l = 1;
    }
    if (state == GLUT_UP) {
      mouse_l = 0;
    }
    break;
  case GLUT_MIDDLE_BUTTON:
    if (state == GLUT_DOWN) {
      mpos[0] = x;
      mpos[1] = y;
      mouse_m = 1;
    }
    if (state == GLUT_UP) {
      mouse_m = 0;
    }
    break;
  case GLUT_RIGHT_BUTTON:
    if (state == GLUT_DOWN) {
      mpos[0] = x;
      mpos[1] = y;
      mouse_r = 1;
    }
    if (state == GLUT_UP) {
      mouse_r = 0;
    }
    break;
  default:
    break;
  }
}

void motion(int x, int y)
{
  double d0;
  double len = 10;

  len = eye_len;

  if(mouse_l == 1 && mouse_m == 1){
    trans[0] += (double)(y-mpos[1])*len/150;
    angle[0] = -(double)(x-mpos[0])*0.2;
  } else  if(mouse_m == 1 || (mouse_l == 1 && mouse_r == 1)){
    trans[1] += (double)(x-mpos[0])*len*.001;
    trans[2] -= (double)(y-mpos[1])*len*.001;
  } else if(mouse_r == 1){
    trans[0] -= (double)(y-mpos[1])*len/150;
    angle[0] =  (double)(x-mpos[0])*0.2;
  } else if(mouse_l == 1){
    d0 = len/50;
    if(d0 > 1.0) d0 = 1.0;
    angle[1] = (double)(y-mpos[1])*d0;
    angle[2] = (double)(x-mpos[0])*d0;
  }
  if(mouse_l == 1 || mouse_m == 1 || mouse_r == 1){
    mpos[0] = x;
    mpos[1] = y;
    //glutPostRedisplay();
  }

}


void keyboard(unsigned char key, int x, int y)
{
  int i,j,k;
  int i0,i1,i2;
  double d0,d1,d2,d3,d4,d5;
  double ang0,ang1,ang2,ang3;

  int c;
  int cf[6];
  double cfr[6];
  double cp[18];

  double l = 2.285852;


  l /= 2;


  if(key == '?'){
    printf("!   : initilize\n");
    printf("q,Q : quit program\n");
    printf("i,I : (print information of postion and angle)\n");
    printf("p,P : (set cliping area)\n");
    printf("a,A : auto mode on/off\n");
    printf("k,K : simulatin cell display on/off\n");
    printf("W   : start output bmp file\n");
    printf("c,C : change backgrand color\n");
    printf("r   : make radius small\n");
    printf("R   : make radius large\n");
    printf("d   : make ditail donw\n");
    printf("D   : make ditail up\n");
    printf("v,V : varbose on/off\n");
    printf("s   : md_step--\n");
    printf("S   : md_step++\n");
    printf("t,T : temp += 100\n");
    printf("g,G : temp -= 100\n");
    printf("y,Y : temp += 10\n");
    printf("h,H : temp -= 10\n");
    printf("z,Z : stop/restart\n");
    printf("0-9 : chage particle number\n");
    printf("n   : create a new positive ion\n");
    printf("m   : create a new negative ion\n");
    printf("N   : create 4 new ions\n");
    printf("M   : create 27 new ions\n");
    printf("SP  : shoot new ion(s)\n");
  }

  if(key == ' ' && c_flg != C_STEP && start_vl <= 0){
      grape_flg = (grape_flg == 0 ? 1:0);
#ifdef CPUONLY
      grape_flg = 0;
#endif
  }
  if(key == '!'){
    if(sc_flg != 1){
#ifdef GL_ON
      glLoadIdentity();
      glGetDoublev(GL_MODELVIEW_MATRIX,m_matrix);
      glGetDoublev(GL_MODELVIEW_MATRIX,i_matrix);
#ifdef SUBWIN
      p_count = 0;
      temp_ymax = 2000;
      for(i = 0; i < DATA_NUM; i++)
	temp_data[i] = 0;
#endif
#endif
      trans[0] = 0;
      trans[1] = 0;
      trans[2] = 0;
    }
    c_flg = 0;
    c_num = 0;
    m_clock = 0;
    set_cd(0);
  }

  if((key >= '1' && key <= '9') && c_flg == 0){
    if(sc_flg != 1){
#ifdef GL_ON
      glLoadIdentity();
      glGetDoublev(GL_MODELVIEW_MATRIX,m_matrix);
      glGetDoublev(GL_MODELVIEW_MATRIX,i_matrix);
#ifdef SUBWIN
      p_count = 0;
      temp_ymax = 2000;
      for(i = 0; i < DATA_NUM; i++)
	temp_data[i] = 0;
#endif
#endif
      trans[0] = 0;
      trans[1] = 0;
      trans[2] = 0;
    }
    firstmalloc = 0;
    np = key-'0';
    npx = np;
    npy = np;
    npz = np;
    c_flg = 0;
    c_num = 0;
    m_clock = 0;
    set_cd(0);
  }

#ifdef NOMAINCLARET
  if(key == 27)	g_glutLoopContinue = false;
#endif

  if(key == 'q' || key == 'Q'){
    if(kflg == 1)
      fclose(fp);
    exit(0);
  }
#ifdef INFO
  if(key == 'i' || key == 'I'){
    printf("(%f,%f,%f)\n"
	   ,trans[0]-trans0[0],trans[1]-trans0[1],trans[2]-trans0[2]);
    printf("(");
    for(i = 0; i < 4; i++){
      for(j = 0; j < 4; j++){
	if(i == 0 && j == 0)
	  printf("%f",m_matrix[i*4+j]-matrix0[i*4+j]);
	else
	  printf(",%f",m_matrix[i*4+j]-matrix0[i*4+j]);
      }
    }
    printf(")\n");
    for(i = 0; i < 3; i++){
      trans0[i] = trans[i];
    }
    for(i = 0; i < 16; i++)
      matrix0[i] = m_matrix[i];
  }
#endif
  //////ORIGINAL
//  if(key == 'f' || key == 'F') bond_flg = bond_flg == 0 ? 1:0;
////////////END ORIGINAL
//
  if(key == 'f' || key == 'F') bond_flg = bond_flg == 0 ? 0:0;
#ifdef GL_ON
  if(key == 'p' || key == 'P'){
    if(clip_flg == 4) clip_flg = 0; else clip_flg++;
  }
  if(key == 'a' || key == 'A') auto_flg = auto_flg == 0 ? 1:0;
  if(key == 'k' || key == 'K') kabe_flg = kabe_flg == 0 ? 1:0;
  if(key == 'W') save_flg = save_flg == 0 ? 1:0;
  if(key == 'c' || key == 'C') clear_color = clear_color == 0 ? 1:0;
  if(key == 'r' && (int)(radius*10+.5) > 1) {radius -= .1;}
  if(key == 'R') {radius += .1; }
  if(key == 'd' && ditail > 5)  {ditail -= 1; }
  if(key == 'D' && ditail < 20) {ditail += 1; }
  if(key == 'v' || key == 'V'){
    if(vflg != 0)
      vflg--;
    else
#ifdef LAP_TIME
      vflg = 3;
#else
      vflg = 1;
#endif
  }
#endif
#if defined(MDGRAPE3) || defined(VTGRAPE)
  if(key == 's' && md_step > 10){ md_step -= 10; md_stepf = 10;}
  if(key == 'S'){ md_step += 10; md_stepf = 10;}
#else
  if(key == 's' && md_step > 1){ md_step -= 1; md_stepf = 10;}
  if(key == 'S'){ md_step += 1; md_stepf = 10;}
#endif
  if(key == 't' || key == 'T'){
    temp += 100;
    rtemp = temp / epsv * kb;
  }
  if(key == 'g' || key == 'G'){
    if(temp > 100){
      temp -= 100;
      rtemp = temp / epsv * kb;
    }
  }
  if(key == 'h' || key == 'H'){
    if(temp > 10){
      temp -= 10;
      rtemp = temp / epsv * kb;
    }
  }
  if(key == 'y' || key == 'Y'){
    temp += 10;
    rtemp = temp / epsv * kb;
  }
  if(key == 'z' || key == 'Z'){
    run_flg *= -1;
#if 0
    if(sc_flg == 0){
      if(run_flg == 1)
        glutIdleFunc(md_run);
      else
        glutIdleFunc(NULL);
    } else {
      if(run_flg == 1)
        glutIdleFunc(md_run);
      else
	glutIdleFunc(NULL);
    }
#endif
  }
  /*
  if(key == '0')
    eye_pos = -1;
  if(key == '-')
    eye_pos = 0;
  if(key == '=')
    eye_pos = 1;
  */
  /*
#ifdef GL_ON
  if((key >= '1' && key <= '5') && c_flg == 0){
    drow_flg[key-'1'] = (drow_flg[key-'1'] == 1 ? 0:1);
  }
#endif
  */

///////////////////////////////////////////////////////
/*  if(key >= '0' && key <= '9' && c_flg == C_STEP){
    start_vl = .3*(key-'0'+1)/10*delt/2e-15;
    velp_flg = 1;
  }*/
////////////////////////////////////////////////////////

  //if((key == 'N' || key == 'M' ||/*key == 'B' ||*/
   //   key == 'n' || key == 'm'/* || key == 'b'*/) && c_flg == 0){
  if(key == 'M' ){

    c_flg = 1;
    w_add = s_add = 0;

    if(key == 'b' || key == 'B'){
      /*c_num = w_site;
      w_add = 1;*/
	printf("This option is disable,kick out the comments-");
	printf("line 1954,1946 and 1945\n");
    } else if(key == 'N'){
      c_num = 4;
      s_add = 4;
      r = 3;
    } else if(key == 'M'){
      c_num = 27;
      s_add = 27;
      r = 9;
    } else {
      c_num = 1;
      s_add = 1;
      r = 1;
    }

    d0 = (i_matrix[0]*(-trans[0])+
          i_matrix[4]*(-trans[1])+
          i_matrix[8]*(-trans[2]));
    d1 = (i_matrix[1]*(-trans[0])+
          i_matrix[5]*(-trans[1])+
          i_matrix[9]*(-trans[2]));
    d2 = (i_matrix[2]*(-trans[0])+
          i_matrix[6]*(-trans[1])+
          i_matrix[10]*(-trans[2]));
    d0 += sideh[0];
    d1 += sideh[1];
    d2 += sideh[2];

    d3 = (i_matrix[0]*(-trans[0]+eye_len-10)+
          i_matrix[4]*(-trans[1])+
          i_matrix[8]*(-trans[2]));
    d4 = (i_matrix[1]*(-trans[0]+eye_len-10)+
          i_matrix[5]*(-trans[1])+
          i_matrix[9]*(-trans[2]));
    d5 = (i_matrix[2]*(-trans[0]+eye_len-10)+
          i_matrix[6]*(-trans[1])+
          i_matrix[10]*(-trans[2]));
    d3 += sideh[0];
    d4 += sideh[1];
    d5 += sideh[2];
    t_cd[0] = d3-d0;
    t_cd[1] = d4-d1;
    t_cd[2] = d5-d2;
    /*
    for(i = 0; i < 4; i++){
      for(j = 0; j < 4; j++){
        printf("%f ",i_matrix[i*4+j]);
      }
      printf("\n");
    }
    */
    /*
    printf("%f %f %f %f  %f %f %f\n",trans[0],trans[1],trans[2],eye_len
           ,d3,d4,d5);
    */

    if(d3 < side[0]-r && d3 > r &&
       d4 < side[1]-r && d4 > r &&
       d5 < side[2]-r && d5 > r){
    } else {

      for(i = 0; i < 6; i++){
        cf[i] = 0;
        cfr[i] = -1;
      }

      i0 = 0;      /* x > */
      cp[i0]   = side[0]-r;
      cp[i0+1] = (cp[i0]-d3)/(d3-d0)*(d4-d1) + d4;
      cp[i0+2] = (cp[i0]-d3)/(d3-d0)*(d5-d2) + d5;
      if(cp[i0+1] > r && cp[i0+1] < side[1]-r &&
         cp[i0+2] > r && cp[i0+2] < side[2]-r) cf[i0/3] = 1;
      i0 += 3;    /* < x */
      cp[i0]   = r;
      cp[i0+1] = (cp[i0]-d3)/(d3-d0)*(d4-d1) + d4;
      cp[i0+2] = (cp[i0]-d3)/(d3-d0)*(d5-d2) + d5;
      if(cp[i0+1] > r && cp[i0+1] < side[1]-r &&
         cp[i0+2] > r && cp[i0+2] < side[1]-r) cf[i0/3] = 1;

      i0 += 3;    /* y > */
      cp[i0+1]   = side[1]-r;
      cp[i0]   = (cp[i0+1]-d4)/(d4-d1)*(d3-d0) + d3;
      cp[i0+2] = (cp[i0+1]-d4)/(d4-d1)*(d5-d2) + d5;
      if(cp[i0]   > r && cp[i0]   < side[0]-r &&
         cp[i0+2] > r && cp[i0+2] < side[2]-r) cf[i0/3] = 1;
      i0 += 3;    /* < y */
      cp[i0+1]   = r;
      cp[i0]   = (cp[i0+1]-d4)/(d4-d1)*(d3-d0) + d3;
      cp[i0+2] = (cp[i0+1]-d4)/(d4-d1)*(d5-d2) + d5;
      if(cp[i0]   > r && cp[i0]   < side[0]-r &&
         cp[i0+2] > r && cp[i0+2] < side[2]-r) cf[i0/3] = 1;

      i0 += 3;   /* z > */
      cp[i0+2]   = side[2]-r;
      cp[i0]   = (cp[i0+2]-d5)/(d5-d2)*(d3-d0) + d3;
      cp[i0+1] = (cp[i0+2]-d5)/(d5-d2)*(d4-d1) + d4;
      if(cp[i0]   > r && cp[i0]   < side[0]-r &&
         cp[i0+1] > r && cp[i0+1] < side[1]-r) cf[i0/3] = 1;
      i0 += 3;
      cp[i0+2]   = r;
      cp[i0]   = (cp[i0+2]-d5)/(d5-d2)*(d3-d0) + d3;
      cp[i0+1] = (cp[i0+2]-d5)/(d5-d2)*(d4-d1) + d4;
      if(cp[i0]   > r && cp[i0]   < side[0]-r &&
         cp[i0+1] > r && cp[i0+1] < side[1]-r) cf[i0/3] = 1;

      for(i = 0; i < 6; i++){
        if(cf[i] == 1){
          cfr[i] = sqrt((cp[i*3]  -d3)*(cp[i*3]  -d3)+
                        (cp[i*3+1]-d4)*(cp[i*3+1]-d4)+
                        (cp[i*3+2]-d5)*(cp[i*3+2]-d5));
        }
      }
      d0 = 10000;
      c = -1;
      for(i = 0; i < 6; i++){
        if(cf[i] == 1 && d0 > cfr[i]){
          d0 = cfr[i];
          c = i;
        }
      }
      if(c == -1) c_num = 0;
      c *= 3;
      d3 = cp[c];
      d4 = cp[c+1];
      d5 = cp[c+2];
    }

    if(key == 'b' || key == 'B'){

    } else {
      if(c_num == 1){
        cd[n3]   = d3;
        cd[n3+1] = d4;
        cd[n3+2] = d5;
        atype[n1] = ((key == 'n') ? 0:1);
      } else if(c_num == 4){
        i = 0;
        l *= 1.1;
        for(i0 = 0; i0 < 2; i0++){
          for(i1 = 0; i1 < 2; i1++){
            d0 = i_matrix[4]*(l*(i0*2-1))
                +i_matrix[8]*(l*(i1*2-1));
            d1 = i_matrix[5]*(l*(i0*2-1))
                +i_matrix[9]*(l*(i1*2-1));
            d2 = i_matrix[6]*(l*(i0*2-1))
                +i_matrix[10]*(l*(i1*2-1));
            cd[n3+i]   = d0+d3;
            cd[n3+i+1] = d1+d4;
            cd[n3+i+2] = d2+d5;
            atype[n1+i/3] = (i0+i1) %2;
            i += 3;
          }
        }
      } else if(c_num == 27){
        i = 0;
        l *= 1.2;
        for(i0 = 0; i0 < 3; i0++){
          for(i1 = 0; i1 < 3; i1++){
            for(i2 = 0; i2 < 3; i2++){
              d0 = i_matrix[0]*(l*2*(i0-1))
                  +i_matrix[4]*(l*2*(i1-1))
                  +i_matrix[8]*(l*2*(i2-1));
              d1 = i_matrix[1]*(l*2*(i0-1))
                  +i_matrix[5]*(l*2*(i1-1))
                  +i_matrix[9]*(l*2*(i2-1));
              d2 = i_matrix[2]*(l*2*(i0-1))
                  +i_matrix[6]*(l*2*(i1-1))
                  +i_matrix[10]*(l*2*(i2-1));
              cd[n3+i]   = d0+d3;
              cd[n3+i+1] = d1+d4;
              cd[n3+i+2] = d2+d5;
              atype[n1+i/3] = (i0+i1+i2) %2;
              i += 3;
            }
          }
        }
      }
    }
  ////////////////////Start Velocity
    start_vl = .3*(9+1)/10*delt/2e-15;
       velp_flg = 1;
  //////////////////////////////////////Shoot
       w_num += w_add;
        w_num3 = w_num*3;
        s_num += s_add;
        s_num3 = s_num*3;
        ws_num = w_num + s_num;
        ws_num3= ws_num*3;

        n1 = s_num + w_num*w_site;
        n3 = n1*3;
        tscale = 1. / 3. /((double)(s_num + w_num*2) - 1);

        r = sqrt(t_cd[0]*t_cd[0] + t_cd[1]*t_cd[1] + t_cd[2]*t_cd[2]);
        d3 = 0;
        for(i = 0; i < c_num*3; i+=3){
          vl[n3-3-i]   = -t_cd[0]/r *start_vl;
          vl[n3-3-i+1] = -t_cd[1]/r *start_vl;
          vl[n3-3-i+2] = -t_cd[2]/r *start_vl;
          d3 += (vl[n3-3-i]*vl[n3-3-i]+vl[n3-2-i]*vl[n3-2-i]+vl[n3-1-i]*vl[n3-1-i])
            *a_mass[atype_mat[atype[(n3-3-i)/3]]];
        }
        d3 *= tscale/hsq;
        rtemp += d3;
        temp = rtemp * epsv /kb;

        c_flg = 0;
        c_num = 0;
        velp_flg = 0;
        start_vl = -1;

  }


/////////////////////////////////////////////////////////////////
  /*
  if(key == ' ' && c_flg == C_STEP && start_vl > 0){

    w_num += w_add;
    w_num3 = w_num*3;
    s_num += s_add;
    s_num3 = s_num*3;
    ws_num = w_num + s_num;
    ws_num3= ws_num*3;

    n1 = s_num + w_num*w_site;
    n3 = n1*3;
    tscale = 1. / 3. /((double)(s_num + w_num*2) - 1);

    r = sqrt(t_cd[0]*t_cd[0] + t_cd[1]*t_cd[1] + t_cd[2]*t_cd[2]);
    d3 = 0;
    for(i = 0; i < c_num*3; i+=3){
      vl[n3-3-i]   = -t_cd[0]/r *start_vl;
      vl[n3-3-i+1] = -t_cd[1]/r *start_vl;
      vl[n3-3-i+2] = -t_cd[2]/r *start_vl;
      d3 += (vl[n3-3-i]*vl[n3-3-i]+vl[n3-2-i]*vl[n3-2-i]+vl[n3-1-i]*vl[n3-1-i])
        *a_mass[atype_mat[atype[(n3-3-i)/3]]];
    }
    d3 *= tscale/hsq;
    rtemp += d3;
    temp = rtemp * epsv /kb;

    c_flg = 0;
    c_num = 0;
    velp_flg = 0;
    start_vl = -1;

}*/
///////////////////////////////////////////////////////////////////////
#ifdef GL_ON
 // if(sc_flg != 1)
  //  glutPostRedisplay();
#endif
}

void single_display(int which)
{
  double d0,d1,d2,d3,d4,d5;
  double mag;
  int i,j;
  int i0;
  char str_buf[256];
  char str_buf2[256];
  GLfloat particle_color[4];
  GLfloat color[4];

  // Prepare Projection matrix without calling reshape
#ifndef NOMAINCLARET
  glViewport(0, 0, fboWidth, fboHeight);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(60.0, fboWidth / fboHeight, 0.1, 1600.0);
  glMatrixMode(GL_MODELVIEW);
#endif

  glClearColor(clear_color, 0.0,0.0 , 0.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glCullFace(GL_BACK);

  glLoadIdentity();
  glPushMatrix();

  d3 = atan((eye_width*which)/eye_len);
  d1 = sin(d3)*eye_len;
  d0 = cos(d3)*eye_len;

  //gluLookAt(d0, d1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);
  gluLookAt(d0, d1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0);

  glTranslated(trans[0], trans[1], trans[2]);

  glPushMatrix();
  glLoadIdentity();
  glRotatef( angle[0],1.0,0.0,0.0);
  glRotatef( angle[1],0.0,1.0,0.0);
  glRotatef( angle[2],0.0,0.0,1.0);

  glMultMatrixd(m_matrix);
  glGetDoublev(GL_MODELVIEW_MATRIX, m_matrix);
  glPopMatrix();

  for(i = 0; i < 16; i++)
    i_matrix[i] = m_matrix[i];
  mat_inv((double(*)[4])i_matrix);

  glMultMatrixd(m_matrix);


  if(kabe_flg == 1){
  	d2 = (i_matrix[0]*(1.0)+i_matrix[4]*(1.0)+i_matrix[8]*(1.0));
    d3 = (i_matrix[1]*(1.0)+i_matrix[5]*(1.0)+i_matrix[9]*(1.0));
    d4 = (i_matrix[2]*(1.0)+i_matrix[6]*(1.0)+i_matrix[10]*(1.0));

    d0 = side0/2;
    d1 = 0.3;
    color[0] = d1; color[1] = d1; color[2] = d1; color[3] = 1.0;
    glMaterialfv(GL_FRONT,GL_AMBIENT_AND_DIFFUSE,color);
    glLineWidth(2.0);
    glBegin(GL_LINES);
    	glNormal3d(d2,d3,d4);
    	glVertex3d(0,0,-d0);
    	glVertex3d(0,0, d0);
    glEnd();
    glBegin(GL_LINES);
    	glNormal3d(d2,d3,d4);
    	glVertex3d(0,-d0,0);
    	glVertex3d(0, d0,0);
    glEnd();
    glBegin(GL_LINES);
    	glNormal3d(d2,d3,d4);
    	glVertex3d(-d0,0,0);
    	glVertex3d( d0,0,0);
    glEnd();
  }



  if(clip_flg == 0){
    glDisable(GL_CLIP_PLANE0);
  } else if(clip_flg == 1){
    glClipPlane(GL_CLIP_PLANE0, clip[0]);
    glClipPlane(GL_CLIP_PLANE1, clip[1]);
    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
  } else if(clip_flg == 2){
    glClipPlane(GL_CLIP_PLANE0, clip[2]);
    glClipPlane(GL_CLIP_PLANE1, clip[3]);
    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
  } else if(clip_flg == 3){
    glClipPlane(GL_CLIP_PLANE0, clip[4]);
    glClipPlane(GL_CLIP_PLANE1, clip[5]);
    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
  } else if(clip_flg == 4){
    glClipPlane(GL_CLIP_PLANE0, clip[0]);
    glClipPlane(GL_CLIP_PLANE1, clip[1]);
    glClipPlane(GL_CLIP_PLANE2, clip[2]);
    glClipPlane(GL_CLIP_PLANE3, clip[3]);
    glClipPlane(GL_CLIP_PLANE4, clip[4]);
    glClipPlane(GL_CLIP_PLANE5, clip[5]);
    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
    glEnable(GL_CLIP_PLANE2);
    glEnable(GL_CLIP_PLANE3);
    glEnable(GL_CLIP_PLANE4);
    glEnable(GL_CLIP_PLANE5);
  }

  angle[0] = 0;
  if(mouse_l == 1 || mouse_m == 1 || mouse_r == 1){
    angle[1] = 0;
    angle[2] = 0;
  }
  if(ini_flg == 1){
    mouse_l = 0;
    ini_flg = 0;
  }

  if(kabe_flg == 1)
    hako(1);


#if defined(MDGRAPE3) || defined(VTGRAPE)
  if(bond_flg == 1 && grape_flg==0){
#else
  if(bond_flg == 1){
#endif
    d2 = (i_matrix[0]*(1.0)+i_matrix[4]*(1.0)+i_matrix[8]*(1.0));
    d3 = (i_matrix[1]*(1.0)+i_matrix[5]*(1.0)+i_matrix[9]*(1.0));
    d4 = (i_matrix[2]*(1.0)+i_matrix[6]*(1.0)+i_matrix[10]*(1.0));
    glNormal3d(d2,d3,d4);
    for(i = 0; i < n1; i++){
      for(j = 0; j < nig_num[i]; j++){
	/*	  bond_drow(i,nig_data[i][j],0, 0.03, 5);*/
    	line_drow(i,nig_data[i*6+j]);
      }
    }
  }


#if 1
///////////////////////////Start Particle Drawing/////////////////////////

  // This routine uses Inter-operability with CUDA. We do not call any memcpy to from GPU to CPU
  // We use velocity to compute color
  if (rendermode == 1 ){

#ifdef INTEROP
  	int n = n3/3;

		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);
		glEnable(GL_COLOR_MATERIAL);

		///////////////////////////
		glBindBuffer(GL_ARRAY_BUFFER, g_possVBO);
		glVertexPointer(3, GL_FLOAT, 0, 0);

		glBindBuffer(GL_ARRAY_BUFFER, g_colorVBO);
		glColorPointer(4,GL_FLOAT,0,0);

		glPointSize(8.0);
		glDrawArrays(GL_POINTS,0,n);

		///////////////////////////////
		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);
		glDisable(GL_COLOR_MATERIAL);

#else
		// This routine needs Position and Velocity from GPU. memcpy to CPU is performed every md_step.
		glEnableClientState(GL_VERTEX_ARRAY);
		glEnableClientState(GL_COLOR_ARRAY);
		glEnable(GL_COLOR_MATERIAL);

		GLuint buf[2];
		glGenBuffers(2, buf);

		int n = n3/3;
		int q=0;

		float *f_pointA,*f_clr;
		unsigned int size 			= (n*3)*sizeof(float);
		unsigned int size_color = (n*4)*sizeof(float);

		f_pointA    = (float*)malloc(n*3*sizeof(float));
		f_clr	  		= (float*)malloc(n*4*sizeof(float));

		///////////////////////////
		for(i=0; i<n3;i+=3){
			if(drow_flg[atype_mat[atype[i/3]]] == 1){
				// Compute coordinates
				f_pointA[i] 	= cd[i]-sideh[0];
				f_pointA[i+1]	= cd[i+1]-sideh[1];
				f_pointA[i+2]	= cd[i+2]-sideh[2];
				// Compute Color
				d0 = (vl[i]*vl[i]+vl[i+1]*vl[i+1]+vl[i+2]*vl[i+2])*500;
				*(f_clr+0+q) = color_table[atype_mat[atype[i/3]]][0]+d0;
				*(f_clr+1+q) = color_table[atype_mat[atype[i/3]]][1]+d0/3;
				*(f_clr+2+q) = color_table[atype_mat[atype[i/3]]][2]+d0/3;
				*(f_clr+3+q) = color_table[atype_mat[atype[i/3]]][3];
				q+=4;
				}
		}
		///////////////////////////
		glBindBuffer(GL_ARRAY_BUFFER, buf[0]);
		glBufferData(GL_ARRAY_BUFFER,size,f_pointA, GL_DYNAMIC_DRAW);
		glVertexPointer(3, GL_FLOAT, 0, 0);


		glBindBuffer(GL_ARRAY_BUFFER, buf[1]);
		glBufferData(GL_ARRAY_BUFFER,size_color,f_clr, GL_DYNAMIC_DRAW);
		glColorPointer(4,GL_FLOAT,0,0);

		glPointSize(8.0);
		glDrawArrays(GL_POINTS,0,n);
		///////////////////////////////
		glDisableClientState(GL_VERTEX_ARRAY);
		glDisableClientState(GL_COLOR_ARRAY);
		glDisable(GL_COLOR_MATERIAL);

		glDeleteBuffers(2, buf);

		free(f_pointA);
		free(f_clr);



#endif

 }

  // This also needs Position and velocity. Uses texture by CPU
  else if(rendermode == 2){
		for(i = 0; i < n3; i += 3){
			if(drow_flg[atype_mat[atype[i/3]]] == 1){
				glPushMatrix();

				glTranslated(cd[i]-sideh[0], cd[i+1]-sideh[1], cd[i+2]-sideh[2]);
				glMultMatrixd(i_matrix);
				glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, black);
				glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, tran);
				TexCirSphere(0.7,atype[i/3]);

				//if(atype[i/3] == 8){
				//	glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE,red);
				//exit(0);
				//}  else{

				d0 = (vl[i]*vl[i]+vl[i+1]*vl[i+1]+vl[i+2]*vl[i+2])*500;
				particle_color[0] = color_table[atype_mat[atype[i/3]]][0]+d0;
				particle_color[1] = color_table[atype_mat[atype[i/3]]][1]+d0/3;
				particle_color[2] = color_table[atype_mat[atype[i/3]]][2]+d0/3;
				particle_color[3] = color_table[atype_mat[atype[i/3]]][3];
				glMaterialfv(GL_FRONT, GL_AMBIENT,particle_color);

				particle_color[0] = color_table[atype_mat[atype[i/3]]][0]+d0/4;
				particle_color[1] = color_table[atype_mat[atype[i/3]]][1]+d0/12;
				particle_color[2] = color_table[atype_mat[atype[i/3]]][2]+d0/12;
				particle_color[3] = color_table[atype_mat[atype[i/3]]][3];
				glMaterialfv(GL_FRONT, GL_DIFFUSE,particle_color);
  //    }
				glPopMatrix();
    	}
  	}
  }

  else if (rendermode == 3){
  	for(i = 0; i < n3; i += 3){
  		if(drow_flg[atype_mat[atype[i/3]]] == 1){
  			glPushMatrix();
  			if(atype[i/3] == 8){
  				glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE,red);
  				exit(0);
  			}  else{

				d0 = (vl[i]*vl[i]+vl[i+1]*vl[i+1]+vl[i+2]*vl[i+2])*500;
				particle_color[0] = color_table[atype_mat[atype[i/3]]][0]+d0;
				particle_color[1] = color_table[atype_mat[atype[i/3]]][1]+d0/3;
				particle_color[2] = color_table[atype_mat[atype[i/3]]][2]+d0/3;
				particle_color[3] = color_table[atype_mat[atype[i/3]]][3];
				glMaterialfv(GL_FRONT, GL_AMBIENT,particle_color);

				particle_color[0] = color_table[atype_mat[atype[i/3]]][0]+d0/4;
				particle_color[1] = color_table[atype_mat[atype[i/3]]][1]+d0/12;
				particle_color[2] = color_table[atype_mat[atype[i/3]]][2]+d0/12;
				particle_color[3] = color_table[atype_mat[atype[i/3]]][3];
				glMaterialfv(GL_FRONT, GL_DIFFUSE,particle_color);

      }

    glTranslated(cd[i]-sideh[0], cd[i+1]-sideh[1], cd[i+2]-sideh[2]);
    glutSolidSphere(radius*r_table[atype_mat[atype[i/3]]], ditail, ditail/2);
    glPopMatrix();
    }
  	}
  }


#endif
///////////////////////////////////////////////////////////////////
//////////////////////END Particle Drawing////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
  glPopMatrix();

  glDisable(GL_DEPTH_TEST);
  if(clip_flg != 0){
    glDisable(GL_CLIP_PLANE0);
    glDisable(GL_CLIP_PLANE1);
    glDisable(GL_CLIP_PLANE2);
    glDisable(GL_CLIP_PLANE3);
    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);
  }

#ifdef LAP_TIME
#if defined(_WIN32) && !defined(__CYGWIN__)
  disp_time0 = disp_time;
  disp_time = (double)timeGetTime()/1000.;
#elif defined(MAC)
  disp_time0 = disp_time;
  disp_time = (double)clock()/60.;
#else
  gettimeofday(&time_v,NULL);
  disp_time0 = disp_time;
  disp_time = (time_v.tv_sec + time_v.tv_usec / 1000000.0);
#endif
#endif

  glDisable(GL_LIGHTING);
/////////////////////////////////////////////////////////
//////////////////////////////////////End for Android part
  if(vflg >= 1){

		d0 = -6.4;
		d1 =  6.2;
		d2 = -12;

		/*    sprintf(str_buf,"T=%.0fK N=%d (W%d N%d)",temp,n1,w_num,s_num);*/
		if(temp_unit_type == 1)
			sprintf(str_buf,"T=%.0fC N=%d",temp-273,n1);
		else
			sprintf(str_buf,"T=%.0fK N=%d",temp,n1);

		if(vflg >= 3 && auto_flg == 0){
			if(grape_flg == 1)
				strcat(str_buf,"  GPU:ON");
      else
      	strcat(str_buf,"  GPU:OFF");
    }

    glColor4fv(moji_c[(int)(clear_color+.5)]);
    glRasterPos3d(d0, d1, d2);
    glCallLists(strlen(str_buf), GL_BYTE, str_buf);

    d1 -= .3;
    if(temp_unit_type == 1)
      sprintf(str_buf,"temp:%4.0fC time:%.3es"
	      ,mtemp*epsv/kb-273,delt*m_clock);
    else
      sprintf(str_buf,"temp:%4.0fK time:%.3es"
	      ,mtemp*epsv/kb,delt*m_clock);

    glRasterPos3d(d0, d1, d2);
    glCallLists(strlen(str_buf), GL_BYTE, str_buf);

    d1 -= .3;
    sprintf(str_buf,"pressure:%.4ePa"
            ,mpres*epsj/(sigma*sigma*sigma));
    glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE
                 ,moji_c[(int)(clear_color+.5)]);
    glRasterPos3d(d0, d1, d2);
    glCallLists(strlen(str_buf), GL_BYTE, str_buf);


    if(velp_flg > 0){
      d1 -= .3;
      d4 = 0;
      d5 = 1. / 3. /((double)(s_num+s_add + (w_num+w_add)*2) - 1);
      for(i = 0; i < c_num; i++)
        d4 += pow(start_vl,2)*a_mass[atype_mat[atype[(n3+i*3)/3]]]*d5/hsq*epsv/kb;
      if(auto_flg == 1)
      	sprintf(str_buf,"Vc = %.0f(m/s) = %.0f(km/h)"
      			,start_vl/h*sigma/tmrdp
      			,start_vl/h*sigma/tmrdp*3.6);
      else
      	sprintf(str_buf,"Vc = %.0f(m/s) = %.0f(km/h) (%4.0fK)"
      			,start_vl/h*sigma/tmrdp
      			,start_vl/h*sigma/tmrdp*3.6,d4);
      glRasterPos3d(d0, d1, d2);
      glCallLists(strlen(str_buf), GL_BYTE, str_buf);
    }

    if(md_stepf > 0){
      d1 -= .3;
      sprintf(str_buf,"md_step=%d",md_step);
      md_stepf--;
      glRasterPos3d(d0, d1, d2);
      glCallLists(strlen(str_buf), GL_BYTE, str_buf);
    }
    if(c_flg == C_STEP  && start_vl <= 0){
      d1 -= .3;
      sprintf(str_buf,"select velocity [0]-[9] keys");
      glRasterPos3d(d0, d1, d2);
      glCallLists(strlen(str_buf), GL_BYTE, str_buf);
    }

    // More Information about Interop and DP
    d1 -= .6;
    sprintf(str_buf,"CUDA related info:");
    glRasterPos3d(d0, d1, d2);
    glCallLists(strlen(str_buf), GL_BYTE, str_buf);

    d1 -= .3;
#ifdef INTEROP
    sprintf(str_buf,"   OpenGL Interop: ON");
#else
    sprintf(str_buf,"   OpenGL Interop: OFF");
#endif
    glRasterPos3d(d0, d1, d2);
    glCallLists(strlen(str_buf), GL_BYTE, str_buf);

    d1 -= .3;
#ifdef DP
    sprintf(str_buf,"   DP: ON");
#else
    sprintf(str_buf,"   DP: OFF");
#endif
    glRasterPos3d(d0, d1, d2);
    glCallLists(strlen(str_buf), GL_BYTE, str_buf);


#ifdef LAP_TIME
    if(vflg >= 2){
      d1 = -6.2;
      if (grape_flg == 0)
    	  sprintf(str_buf,"%.8fs/step %.3fGflops",md_time-md_time0,
    			  (double)n1*(double)n1*78/(md_time-md_time0)*1e-9);
      if (grape_flg == 1)
    	  sprintf(str_buf,"%.8fs/step %.3fGflops",md_time-md_time0,
    			  (double)n1*(double)n1*78/(md_time-md_time0)*1e-9*md_step);//md_step from loop in mr3.cu

      glRasterPos3d(d0, d1, d2);
      glCallLists(strlen(str_buf), GL_BYTE, str_buf);
      d1 -= .3;
      sprintf(str_buf,"%.8fs/frm %.3ffrm/s",(disp_time-disp_time0)
	      ,1./(disp_time-disp_time0));
      glRasterPos3d(d0, d1, d2);
      glCallLists(strlen(str_buf), GL_BYTE, str_buf);
    }
#endif

    // For Drawing the sphere sample atom type
#if 0
		glEnable(GL_LIGHTING);
		d0 = 1.8;
		d1 = -2.4;
		d2 = -10;
		glPushMatrix();
		glTranslated(d0,d1,d2);
		d3 = temp*0.00016;
		particle_color[0] = color_table[atype_mat[0]][0]+d3;
		particle_color[1] = color_table[atype_mat[0]][1]+d3/3;
		particle_color[2] = color_table[atype_mat[0]][2]+d3/3;
		particle_color[3] = color_table[atype_mat[0]][3];
		glMaterialfv(GL_FRONT, GL_AMBIENT,particle_color);
		particle_color[0] = color_table[atype_mat[0]][0]+d3/4;
		particle_color[1] = color_table[atype_mat[0]][1]+d3/12;
		particle_color[2] = color_table[atype_mat[0]][2]+d3/12;
		particle_color[3] = color_table[atype_mat[0]][3];
		glMaterialfv(GL_FRONT, GL_DIFFUSE,particle_color);
		glutSolidSphere(radius*r_table[atype_mat[0]]/3.5, ditail, ditail/2);
		glPopMatrix();
		glDisable(GL_LIGHTING);
		medium_font(d0-0.08,d1-.04,d2,"Na");
		small_font(d0+0.08,d1+.04,d2,"+");

		glEnable(GL_LIGHTING);
		d0 += 0.5;
		glPushMatrix();
		glTranslated(d0,d1,d2);
		particle_color[0] = color_table[atype_mat[1]][0]+d3;
		particle_color[1] = color_table[atype_mat[1]][1]+d3/3;
		particle_color[2] = color_table[atype_mat[1]][2]+d3/3;
		particle_color[3] = color_table[atype_mat[1]][3];
		glMaterialfv(GL_FRONT, GL_AMBIENT,particle_color);
		particle_color[0] = color_table[atype_mat[1]][0]+d3/4;
		particle_color[1] = color_table[atype_mat[1]][1]+d3/12;
		particle_color[2] = color_table[atype_mat[1]][2]+d3/12;
		particle_color[3] = color_table[atype_mat[1]][3];
		glMaterialfv(GL_FRONT, GL_DIFFUSE,particle_color);
		glutSolidSphere(radius*r_table[atype_mat[1]]/3.5, ditail, ditail/2);
		glPopMatrix();
		glDisable(GL_LIGHTING);
		medium_font(d0-0.06,d1-.04,d2,"Cl");
		small_font(d0+0.07,d1+.04,d2,"-");
#endif
  }
  glDisable(GL_LIGHT0);
  glDisable(GL_CULL_FACE);

#ifndef NOMAINCLARET
 glutPostRedisplay();
 glutSwapBuffers();
#endif
}
  void display(void)
{
#if STEREO == 1
  glDrawBuffer(GL_BACK_LEFT);
  single_display(-1);

  glDrawBuffer(GL_BACK_RIGHT);
  single_display(1);
#else
  single_display(eye_pos);
#endif

 // glutSwapBuffers();

}

///////////////////////End OpenGL Functions//////////////////////////////////
