//////////////////////////////////////////////////////////////////////////
//
// Test program for the class RBFInterpolator
// 
// 2009 Karsten Noe
//
// Read the blog at cg.alexandra.dk for more information
//
//////////////////////////////////////////////////////////////////////////

#include <stdarg.h>
#include <iostream>
#include "glew/glew.h"												// Header File For The OpenGL32 Library
#include <GL/glut.h>												// Header File for Glut
#include "math.h"													// math functions

#pragma comment( lib, "opengl32.lib" )								// Search For openhl32.lib while linking
#pragma comment( lib, "glu32.lib" )									// Search For glu32.lib while linking
#pragma comment( lib, "glut32.lib" )								// Search For glut32.lib while linking
#pragma comment( lib, "glew32.lib" )								// Search For glew32.lib while linking

#include "ArcBall.h"												// ArcBall for navigation
#include "TriangleMesh.h"											// the structure in which the loaded mesh is stored
#include "RBFInterpolator.h"										// interpolation by radial basis functions
#include "NeHe_IPicture.h"											// Texture loader

#define M_PI 3.1415926535f

using namespace std;

/* windows size and position constants */
int GL_WIN_WIDTH = 600;
int GL_WIN_HEIGHT = 600;
int GL_WIN_INITIAL_X = 0;
int GL_WIN_INITIAL_Y = 0;

Matrix4fT   Transform   = {  1.0f,  0.0f,  0.0f,  0.0f,				// Arcball navigation: Final Transform
                             0.0f,  1.0f,  0.0f,  0.0f,
                             0.0f,  0.0f,  1.0f,  0.0f,
                             0.0f,  0.0f,  0.0f,  1.0f };

Matrix3fT   LastRot     = {  1.0f,  0.0f,  0.0f,					// Arcball navigation: Last Rotation
                             0.0f,  1.0f,  0.0f,
                             0.0f,  0.0f,  1.0f };

Matrix3fT   ThisRot     = {  1.0f,  0.0f,  0.0f,					// Arcball navigation: This Rotation
                             0.0f,  1.0f,  0.0f,
                             0.0f,  0.0f,  1.0f };

ArcBallT    ArcBall(GL_WIN_WIDTH, GL_WIN_HEIGHT );				    // NEW: ArcBall Instance
Point2fT    MousePt;												// Arcball navigation: Current Mouse Point
bool		bMouseRDown = false;
//bool        isClicked  = false;										// Arcball navigation: Clicking The Mouse?
//bool        isDragging = false;					                    // Arcball navigation: Dragging The Mouse?

RenderTriangleMesh* pighinMeshRenderer;									// Used for rendering a TriangleMesh.
//RenderTriangleMesh* pighinHairRenderer;
RenderTriangleMesh* pighinFaceRenderer;
RenderTriangleMesh* samsungMeshRenderer;
//RenderTriangleMesh* samsungHairRenderer;
RenderTriangleMesh* samsungFaceRenderer;

vector<float> controlPointPosX;										// X-coordinates of control points
vector<float> controlPointPosY;										// Y-coordinates of control points
vector<float> controlPointPosZ;										// Z-coordinates of control points

vector<float> controlPointDisplacementX;							// Displacement of control points in the X-direction
vector<float> controlPointDisplacementY;							// Displacement of control points in the Y-direction
vector<float> controlPointDisplacementZ;							// Displacement of control points in the Z-direction

vector<float> destinationPointPosX;									// Displacement of control points in the X-direction
vector<float> destinationPointPosY;									// Displacement of control points in the Y-direction
vector<float> destinationPointPosZ;									// Displacement of control points in the Z-direction

TriangleMesh *pighinMeshUndeformed;									// The original triangle mesh
TriangleMesh *pighinMeshDeformed;									// This mesh will be continuously deformed and rendered
//TriangleMesh *pighinHairUndeformed;
//TriangleMesh *pighinHairDeformed;
TriangleMesh *pighinFaceUndeformed;
TriangleMesh *pighinFaceDeformed;
TriangleMesh *samsungMeshUndeformed;
TriangleMesh *samsungMeshDeformed;
//TriangleMesh *samsungHairUndeformed;
//TriangleMesh *samsungHairDeformed;
TriangleMesh *samsungFaceUndeformed;
TriangleMesh *samsungFaceDeformed;
TriangleMesh *WWMeshUndeformed;
TriangleMesh *SHMeshUndeformed;
TriangleMesh *JQMeshUndeformed;

BoundingSphere boundingSphere;										// Bounding sphere to give information about the size of the object.

RBFInterpolator pighinRBFX, pighinRBFY, pighinRBFZ;					// Our wonderful radial basis function interpolation things :)
RBFInterpolator samsungRBFX, samsungRBFY, samsungRBFZ;

GLuint	imgfile[5];													// Number of Textures (image files)
bool normalUpdatesDesired = true;									// When the program is running press 'n' to disable normal updates
bool bShowCtrPts = false;											// Show control points for deforming the mesh. press 'b'
bool bTransparency = false;											// Make the mesh transparent. press 'm'
bool bShowTexture = false;											// Map Texture. press 't'

bool bPause = false; /* pause time */
static float morphtime = 0;	/* Mesh morphing time with cos(time) */
Vector3 MaFace;	/* Vertex finder */
unsigned int iMeshToDeform = 1;	/* mesh chooser 1:Pighin, 2:Samsung */
unsigned int iTexToRender = 0;	/* texture chooser 1:Samsung, 2:WW, 3:SH, 4:JQ */

// Texture fit for Pighin Pighin:{.022,.024,.503,.537},	Samsung:{0,0,0,0}, WW:{.018,.0245,.623,.397}, SH:{.0155,.0205,.548,.0402}
// Texture fit for Samsung Pighin:{0,0,0,0}	Samsung:{.022,.023,.502,.588}, WW:{.019,.025,.623,.392}, SH:{.018,.0215,.0548,.402}
GLfloat texscaleX=0.022;
GLfloat texscaleY=0.026;
GLfloat texlocX=0.50;
GLfloat texlocY=0.53;
GLfloat hairscaleX=-0.0235;
GLfloat hairscaleY=0.032;
GLfloat hairlocX=0.55;
GLfloat hairlocY=0.43;

// we want control points which we place at different vertex positions
const int numControlPoints = 24;	/*35*/	// can't be less than 7

const int pighinControlArray[numControlPoints] = {/*3974,*/3250,/*3072,3079,*/ /*536,*/226,/*0,55,*/
	3130,3977,3620,3270,  93,956,597,247,
	/*2842,*/332,2443,2609,2411,3097,76,  4954,2360,1931,1547,2309,1506,  1833,2207};
	//4475,2049,  3300,213
const int samsungControlArray[numControlPoints] = {/*201,*/265,/*270,199,*/ /*852,*/995,/*930,849,*/
	1247,1248,100,591,  792,671,1252,1250,
	/*616,*/574,610,785,578,32,812,  92,1228,800,1236,1242,620,  619,607};
	//140,705, 89,662
const int WWControlArray[numControlPoints] = {0,1,2,
	3,4,5,6,  7,8,9,10,
	11,12,13,14,15,16,  17,18,19,20,21,22,  23};
const int SHControlArray[numControlPoints] = {0,1,2,
	3,4,5,6,  7,8,9,10,
	11,12,13,14,15,16,  17,18,19,20,21,22,  23};	/* Sadly, there's no code to just copy an array to another */
const int JQControlArray[numControlPoints] = {0,1,2,
	3,4,5,6,  7,8,9,10,
	11,12,13,14,15,16,  17,18,19,20,21,22,  23};

void InitializeGlew()
{

	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
		/* Problem: glewInit failed, something is seriously wrong. */
		fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	}

	if (!(GLEW_ARB_vertex_buffer_object
		&& GLEW_EXT_framebuffer_object))
	{
		cerr << "Extensions not supported" << endl;
		exit(1);

	}

}

void DrawNet(GLfloat size, GLint LinesX, GLint LinesZ)
{
	/*Camera1*/
	glBegin(GL_LINES);
	for (int xc = 0; xc < LinesX; xc++)
	{
		glVertex3f(	-size / 2.0 + xc / (GLfloat)(LinesX-1)*size,
					0.0,
					size / 2.0);
		glVertex3f(	-size / 2.0 + xc / (GLfloat)(LinesX-1)*size,
					0.0,
					size / -2.0);
	}
	for (int zc = 0; zc < LinesX; zc++)
	{
		glVertex3f(	size / 2.0,
					0.0,
					-size / 2.0 + zc / (GLfloat)(LinesZ-1)*size);
		glVertex3f(	size / -2.0,
					0.0,
					-size / 2.0 + zc / (GLfloat)(LinesZ-1)*size);
	}
	glEnd();
}

// Draw grids using DrawNet
void DrawGrids ()
{
	GLfloat size = 2.0;
	GLint LinesX = 15;
	GLint LinesZ = 15;
	GLfloat halfsize = size / 2.0;

	glDisable(GL_LIGHTING);
	glPushMatrix();
		glScalef(60.0, 60.0, 60.0);
		glColor3f(0.5, 0.5, 0.5);
		glPushMatrix();	// top and bottom
			glTranslatef(0.0,-halfsize ,0.0);
			DrawNet(size,LinesX,LinesZ);
			glTranslatef(0.0,size,0.0);
			DrawNet(size,LinesX,LinesZ);
		glPopMatrix();
		glPushMatrix();	// left and right
			glTranslatef(-halfsize,0.0,0.0);	
			glRotatef(90.0,0.0,0.0,halfsize);
			DrawNet(size,LinesX,LinesZ);
			glTranslatef(0.0,-size,0.0);
			DrawNet(size,LinesX,LinesZ);
		glPopMatrix();
		glPushMatrix(); // front and back
			glTranslatef(0.0,0.0,-halfsize);	
			glRotatef(90.0,halfsize,0.0,0.0);
			DrawNet(size,LinesX,LinesZ);
			glTranslatef(0.0,size,0.0);
			DrawNet(size,LinesX,LinesZ);
		glPopMatrix();
	glPopMatrix();
	glEnable(GL_LIGHTING);
}

/* http://nehe.gamedev.net/data/articles/article.asp?article=13 */
/* http://lists.apple.com/archives/Mac-opengl/2001/Sep/msg00077.html */
Vector3 GetGLPos(int x, int y)
{
	/* GetGLPos */
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	GLfloat winX, winY, winZ;
	GLdouble posX, posY, posZ;
	//GLdouble nearX, nearY, nearZ;
	//GLdouble farX, farY, farZ;

	glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
	glGetDoublev( GL_PROJECTION_MATRIX, projection );
	glGetIntegerv( GL_VIEWPORT, viewport );
	
	winX = (float)x;
	winY = (float)viewport[3] - (float)y;
	glReadPixels( x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
	gluUnProject( winX, winY, winZ, modelview, projection, viewport, &posX, &posY, &posZ);
	
	return Vector3(posX, posY, posZ);
}

/* http://www.lighthouse3d.com/opengl/maths/index.php?raytriint */
//#define vector(a,b,c) \
//	(a)[0] = (b)[0] - (c)[0];	\
//	(a)[1] = (b)[1] - (c)[1];	\
//	(a)[2] = (b)[2] - (c)[2];
//
//int rayIntersectsTriangle(float *p, float *d, float *v0, float *v1, float *v2) {
//
//	float e1[3],e2[3],h[3],s[3],q[3];
//	float a,f,u,v;
//	
//	vector(e1,v1,v0);
//	vector(e2,v2,v0);
//	crossProduct(h,d,e2);
//	a = innerProduct(e1,h);
//	
//	if (a > -0.00001 && a < 0.00001)
//		return(false);
//	
//	f = 1/a;
//	vector(s,p,v0);
//	u = f * (innerProduct(s,h));
//	
//	if (u < 0.0 || u > 1.0)
//		return(false);
//	
//	crossProduct(q,s,e1);
//	v = f * innerProduct(d,q);
//	if (v < 0.0 || u + v > 1.0)
//		return(false);
//	// at this stage we can compute t to find out where 
//	// the intersection point is on the line
//	t = f * innerProduct(e2,q);
//	if (t > 0.00001) // ray intersection
//		return(true);
//	else // this means that there is a line intersection  
//		 // but not a ray intersection
//		 return (false);
//}

void InitializeOpenGL()					// Any GL Init Code & User Initialization Goes Here
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45,640/480,1,-100);
	glMatrixMode(GL_MODELVIEW);

	// Start Of User Initialization
	bMouseRDown	= false;
    //isClicked   = false;								            // NEW: Clicking The Mouse?
    //isDragging  = false;							                // NEW: Dragging The Mouse?

	glClearColor(0.2,0.2,0.2,1.0);									// clear color is gray
	//glClearDepth (1.0f);											// Depth Buffer Setup
	//glDepthFunc (GL_LEQUAL);										// The Type Of Depth Testing (Less Or Equal)
	glEnable (GL_DEPTH_TEST);										// Enable Depth Testing
	glShadeModel (GL_SMOOTH);										

	/*GLdouble zer_pl[]={0,0,1, -1};
	GLdouble one_pl[]={0,0,1, 1};
	glMatrixMode(GL_TEXTURE);
	glClipPlane(GL_CLIP_PLANE0, zer_pl);
	glEnable(GL_CLIP_PLANE0);
	glClipPlane(GL_CLIP_PLANE1, one_pl);
	glEnable(GL_CLIP_PLANE1);
	glMatrixMode(GL_MODELVIEW);
	glDisable(GL_CLIP_PLANE0);
	glDisable(GL_CLIP_PLANE1);*/
	/* Now the Texture Stuff */
	if (!BuildTexture("Data/head.jpg", imgfile[0]))
		perror("OMG");
	if (!BuildTexture("Data/WW_center.jpg", imgfile[1]))
		perror("OMG");
	if (!BuildTexture("Data/SH_center.jpg", imgfile[2]))
		perror("OMG");
	if (!BuildTexture("Data/JQ_center.jpg", imgfile[3]))
		perror("OMG");
	if (!BuildTexture("Data/hair.jpg", imgfile[4]))
		perror("OMG");

	glBindTexture(GL_TEXTURE_2D,imgfile[4]);	// not needed (i think)
	GLfloat hair_planeS[]={0,hairscaleX,0, hairlocX};
	GLfloat hair_planeT[]={hairscaleY,0,0, hairlocY};
	glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);	// GL_OBJECT_LINEAR, GL_EYE_LINEAR
	glTexGenfv(GL_S, GL_OBJECT_PLANE, hair_planeS);				// GL_OBJECT_PLANE, GL_EYE_PLANE
	glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);	
	glTexGenfv(GL_T, GL_OBJECT_PLANE, hair_planeT);
	//glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_BORDER);
	//glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_BORDER);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);	// GL_MODULATE, GL_ADD

	glBindTexture(GL_TEXTURE_2D,imgfile[0]);
	GLfloat ref_planeS[]={texscaleX,0,0, texlocX};
	GLfloat ref_planeT[]={0,texscaleY,0, texlocY};
	glTexGeni(GL_S, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);	// GL_OBJECT_LINEAR, GL_EYE_LINEAR
	glTexGenfv(GL_S, GL_OBJECT_PLANE, ref_planeS);				// GL_OBJECT_PLANE, GL_EYE_PLANE
	glTexGeni(GL_T, GL_TEXTURE_GEN_MODE, GL_OBJECT_LINEAR);	
	glTexGenfv(GL_T, GL_OBJECT_PLANE, ref_planeT);
	//glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_BORDER);
	//glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_BORDER);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);	// GL_MODULATE, GL_ADD
	glEnable(GL_TEXTURE_GEN_S);
	glEnable(GL_TEXTURE_GEN_T);
	//glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);

	glBlendFunc(GL_SRC_ALPHA,GL_ONE);		// Blending Function For Translucency Based On Source Alpha Value ( NEW )
	
	//GLfloat lightpos[] = {0,0,0,1};		/* lol fixed the ghost possession */
	//glLightfv(GL_LIGHT0, GL_POSITION, lightpos);					// set light position
	//GLfloat lightcolor[] = {1,1,1,1};
	//glLightfv(GL_LIGHT0, GL_DIFFUSE, lightcolor);					// set light color

	//GLfloat light_ambient[] = {0.2, 0.2, 0.2, 1.0};
	//GLfloat light_diffuse[] = {0.7, 0.7, 0.7, 1.0};
	//GLfloat light_specular[] = {0.7, 0.7, 0.7, 1.0};
	GLfloat light_ambient[] = {0.2, 0.2, 0.2, 1.0};
	GLfloat light_diffuse[] = {0.6, 0.6, 0.6, 1.0};
	GLfloat light_specular[] = {0.6, 0.6, 0.6, 1.0};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);

	glEnable(GL_LIGHT0);											// Enable Default Light
	glEnable(GL_LIGHTING);											// Enable Lighting

	glPointSize(3);	/**/
}

void glutResize (int w, int h)
{ 
	GL_WIN_WIDTH = w;
	GL_WIN_HEIGHT = h;
	ArcBall.setBounds(GL_WIN_WIDTH, GL_WIN_HEIGHT);
	glViewport (0, 0, GL_WIN_WIDTH, GL_WIN_HEIGHT);
}

void glutKeyboard (unsigned char key, int x, int y)
{ 
    switch (key)
    {
		case ' ':	// When the Space Bar is pressed
			bPause=!bPause;
		break;
		case 's':
			bPause=false;
			if (morphtime!=M_PI)
				morphtime = M_PI;
			else
				morphtime = 0;
		break;
		case 'r':
			bPause=false;
			morphtime=0;
			Matrix3fSetIdentity(&LastRot);	// three lines here reset all rotations
			Matrix3fSetIdentity(&ThisRot);								
			Matrix4fSetRotationScaleFromMatrix3f(&Transform, &ThisRot);	
		break;
		case 'b':
			bShowCtrPts=!bShowCtrPts;
		break;
		case 'm':
			bTransparency=!bTransparency;
		break;
		case 't':
			bShowTexture=!bShowTexture;
		break;
		case '1':
			for (unsigned int i = 0; i<controlPointPosX.size();  i++)
			{
				Vector3 pos = pighinMeshUndeformed->getParticles()[pighinControlArray[i]].getPos();
				destinationPointPosX[i]=pos[0];
				destinationPointPosY[i]=pos[1];
				destinationPointPosZ[i]=pos[2];
			}
			bPause=false;
			morphtime=0;
			iTexToRender=0;
			//glBindTexture(GL_TEXTURE_2D,imgfile[iTexToRender]);
		break;
		case '2':
			for (unsigned int i = 0; i<controlPointPosX.size();  i++)
			{
				Vector3 pos = samsungMeshUndeformed->getParticles()[samsungControlArray[i]].getPos();
				destinationPointPosX[i]=pos[0];
				destinationPointPosY[i]=pos[1];
				destinationPointPosZ[i]=pos[2];
			}
			bPause=false;
			morphtime=0;
			glBindTexture(GL_TEXTURE_2D,imgfile[0]);
			iTexToRender=0;
		break;
		case '3':
			for (unsigned int i = 0; i<controlPointPosX.size();  i++)
			{
				Vector3 pos = WWMeshUndeformed->getParticles()[WWControlArray[i]].getPos();
				destinationPointPosX[i]=pos[0];
				destinationPointPosY[i]=pos[1];
				destinationPointPosZ[i]=pos[2];
			}
			bPause=false;
			morphtime=0;	// morphtime=M_PI to make the texture show when key is pressed
			glBindTexture(GL_TEXTURE_2D,imgfile[1]);
			iTexToRender=1;
		break;
		case '4':
			for (unsigned int i = 0; i<controlPointPosX.size();  i++)
			{
				Vector3 pos = SHMeshUndeformed->getParticles()[SHControlArray[i]].getPos();
				destinationPointPosX[i]=pos[0];
				destinationPointPosY[i]=pos[1];
				destinationPointPosZ[i]=pos[2];
			}
			bPause=false;
			morphtime=0;
			glBindTexture(GL_TEXTURE_2D,imgfile[2]);
			iTexToRender=2;
		break;
		case '5':
			for (unsigned int i = 0; i<controlPointPosX.size();  i++)
			{
				Vector3 pos = JQMeshUndeformed->getParticles()[JQControlArray[i]].getPos();
				destinationPointPosX[i]=pos[0];
				destinationPointPosY[i]=pos[1];
				destinationPointPosZ[i]=pos[2];
			}
			bPause=false;
			morphtime=0;
			glBindTexture(GL_TEXTURE_2D,imgfile[3]);
			iTexToRender=3;
		break;

		//case 'w':
		//	//zpos-=1.0;
		//	GLfloat xrotrad, yrotrad;
		//	yrotrad = (yrot / 180 * 3.141592654f);
		//	xrotrad = (xrot / 180 * 3.141592654f); 
		//	xpos += GLfloat(sin(yrotrad));
		//	zpos -= GLfloat(cos(yrotrad));
		//	ypos -= GLfloat(sin(xrotrad));
		//break;
		//case 's':
		//	//zpos+=1.0;
		//	GLfloat xrotrad, yrotrad;
		//	yrotrad = (yrot / 180 * 3.141592654f);
		//	xrotrad = (xrot / 180 * 3.141592654f); 
		//	xpos -= GLfloat(sin(yrotrad)) * 0.5;
		//	zpos += GLfloat(cos(yrotrad)) * 0.5;
		//	ypos += GLfloat(sin(xrotrad)) * 0.5;
		//break;
		//case 'd':
		//	//xpos+=1.0;
		//	GLfloat yrotrad;
		//	yrotrad = (yrot / 180 * 3.141592654f);
		//	xpos += GLfloat(cos(yrotrad));
		//	zpos += GLfloat(sin(yrotrad));
		//break;
		//case 'a':
		//	//xpos-=1.0;
		//	GLfloat yrotrad;
		//	yrotrad = (yrot / 180 * 3.141592654f);
		//	xpos -= GLfloat(cos(yrotrad));
		//	zpos -= GLfloat(sin(yrotrad));
		//break;
		case 'z':
			Matrix4fMulRotationScale(&Transform, 1.2);
		break;
		case 'x':
			Matrix4fMulRotationScale(&Transform, 0.8);
		break;
		case '7':
			hairscaleX+=0.0005;
		break;
		case '8':
			hairscaleX-=0.0005;
		break;
		case '9':
			hairscaleY+=0.0005;
		break;
		case '0':
			hairscaleY-=0.0005;
		break;
		case 'u':
			hairlocX+=0.005;
		break;
		case 'i':
			hairlocX-=0.005;
		break;
		case 'o':
			hairlocY+=0.005;
		break;
		case 'p':
			hairlocY-=0.005;
		break;
		case 'c':
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); /* uncomment this function if you only want to draw wireframe model */
		break;
		case 'g':
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); /* uncomment this function if you only want to draw wireframe model */
		break;
		case 'v':
			glPolygonMode(GL_FRONT_AND_BACK, GL_POINT); /* uncomment this function if you only want to draw wireframe model */
		break;
		case 'n': 
			normalUpdatesDesired = !normalUpdatesDesired;
		break;
		case 'l':
			printf("\n Mesh:%d  Texture:%d - %f, %f, %f, %f", iMeshToDeform, iTexToRender, hairscaleX, hairscaleY, hairlocX, hairlocY);
		break;
		case 27: // ESC
			exit (1);
		break;
    }

	if (iMeshToDeform==1 && key=='1')
	{	texscaleX=0.022;	texscaleY=0.026;	texlocX=0.503;		texlocY=0.527;
		hairscaleX=-0.0235;	hairscaleY=0.032;	hairlocX=0.553;		hairlocY=0.427;}
	if (iMeshToDeform==1 && key=='2')
	{	texscaleX=0.022;	texscaleY=0.0225;	texlocX=0.502;		texlocY=0.593;
		hairscaleX=-0.0295;	hairscaleY=0.029;	hairlocX=0.537;		hairlocY=0.458;}
	if (iMeshToDeform==1 && key=='3')
	{	texscaleX=0.018;	texscaleY=0.0245;	texlocX=0.623;		texlocY=0.397;
		hairscaleX=-0.021;	hairscaleY=0.026;	hairlocX=0.417;		hairlocY=0.460;}
	if (iMeshToDeform==1 && key=='4')
	{	texscaleX=0.0155;	texscaleY=0.0205;	texlocX=0.548;		texlocY=0.402;
		hairscaleX=-0.020;	hairscaleY=0.0315;	hairlocX=0.417;		hairlocY=0.333;}
	if (iMeshToDeform==1 && key=='5')
	{	texscaleX=0.017;	texscaleY=0.020;	texlocX=0.553;		texlocY=0.432;
		hairscaleX=-0.0155;	hairscaleY=0.0275;	hairlocX=0.507;		hairlocY=0.428;}

	if (iMeshToDeform==2 && key=='1')
	{	texscaleX=0.022;	texscaleY=0.026;	texlocX=0.503;		texlocY=0.527;
		hairscaleX=-0.025;	hairscaleY=0.031;	hairlocX=0.552;		hairlocY=0.47;}
	if (iMeshToDeform==2 && key=='2')
	{	texscaleX=0.022;	texscaleY=0.0225;	texlocX=0.502;		texlocY=0.593;
		hairscaleX=-0.0305;	hairscaleY=0.035;	hairlocX=0.537;		hairlocY=0.40;}
	if (iMeshToDeform==2 && key=='3')
	{	texscaleX=0.019;	texscaleY=0.025;	texlocX=0.623;		texlocY=0.392;
		hairscaleX=-0.023;	hairscaleY=0.027;	hairlocX=0.432;		hairlocY=0.463;}
	if (iMeshToDeform==2 && key=='4')
	{	texscaleX=0.018;	texscaleY=0.0215;	texlocX=0.548;		texlocY=0.402;
		hairscaleX=-0.022;	hairscaleY=0.0285;	hairlocX=0.452;		hairlocY=0.368;}
	if (iMeshToDeform==2 && key=='5')
	{	texscaleX=0.017;	texscaleY=0.020;	texlocX=0.553;		texlocY=0.432;
		hairscaleX=-0.0145;	hairscaleY=0.0255;	hairlocX=0.487;		hairlocY=0.438;}

	glutPostRedisplay();
}

void glutMouse(int button, int state, int x, int y)
{
	MousePt.s.X=x;
	MousePt.s.Y=y;

	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		MaFace = GetGLPos(x, y);
		bMouseRDown = false;
	}
	if (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN)
	{
		LastRot = ThisRot;										// Set Last Static Rotation To Last Dynamic One
		ArcBall.click(&MousePt);								// Update Start Vector And Prepare For Dragging
		bMouseRDown = true;
	}
	else
	{
		bMouseRDown = false;
	}
}

void glutMotion(int x, int y)									// Perform Motion Updates Here
{
	MousePt.s.X=x;
	MousePt.s.Y=y;

	Quat4fT     ThisQuat;

	if(bMouseRDown == true)
	{
		ArcBall.drag(&MousePt, &ThisQuat);						// Update End Vector And Get Rotation As Quaternion
		Matrix3fSetRotationFromQuat4f(&ThisRot, &ThisQuat);		// Convert Quaternion Into Matrix3fT
		Matrix3fMulMatrix3f(&ThisRot, &LastRot);				// Accumulate Last Rotation Into This One
		Matrix4fSetRotationFromMatrix3f(&Transform, &ThisRot);	// Set Our Final Transform's Rotation From This One
	} 
}

void glutMenu(int entryID)
{
	if (entryID==1)
	{
		iMeshToDeform=1;
		for (unsigned int i = 0; i<controlPointPosX.size(); i++)
		{
			Vector3 pos = pighinMeshUndeformed->getParticles()[pighinControlArray[i]].getPos();
			controlPointPosX[i]=pos[0];
			controlPointPosY[i]=pos[1];
			controlPointPosZ[i]=pos[2];
		}
	}
	else if (entryID==2)
	{
		iMeshToDeform=2;
		for (unsigned int i = 0; i<controlPointPosX.size(); i++)
		{
			Vector3 pos = samsungMeshUndeformed->getParticles()[samsungControlArray[i]].getPos();
			controlPointPosX[i]=pos[0];
			controlPointPosY[i]=pos[1];
			controlPointPosZ[i]=pos[2];
		}
	}
	bPause=false;
	morphtime=0;
}

// Code for deforming the mesh 'initialObject' based on the current interpolation functions (global variables). 
// The deformed vertex positions will be stored in the mesh 'res'
// The triangle connectivity is assumed to be already correct in 'res'  
void deformObject(TriangleMesh* res, TriangleMesh* initialObject, RBFInterpolator rbfX, RBFInterpolator rbfY, RBFInterpolator rbfZ)
{
	for (unsigned int i = 0; i < res->getParticles().size(); i++)
	{
		Vector3 oldpos = initialObject->getParticles()[i].getPos();

		Vector3 newpos;
		newpos[0] = oldpos[0] + rbfX.interpolate(oldpos[0], oldpos[1], oldpos[2]);
		newpos[1] = oldpos[1] + rbfY.interpolate(oldpos[0], oldpos[1], oldpos[2]);
		newpos[2] = oldpos[2] + rbfZ.interpolate(oldpos[0], oldpos[1], oldpos[2]);

		res->getParticles()[i].setPos(newpos);
	}
}

void glutDisplay(void)
{
	if (bPause) {
		morphtime = morphtime + 0.05;
	}
	
	// move control points
	for (unsigned int i = 0; i<controlPointPosX.size(); i++ )
	{
		controlPointDisplacementX[i] = (destinationPointPosX[i]-controlPointPosX[i])*(-cosf(morphtime)/2+0.5);	// -cosine()/2+1/2 go from 0 to 1
		controlPointDisplacementY[i] = (destinationPointPosY[i]-controlPointPosY[i])*(-cosf(morphtime)/2+0.5);
		controlPointDisplacementZ[i] = (destinationPointPosZ[i]-controlPointPosZ[i])*(-cosf(morphtime)/2+0.5);
	}

	const GLfloat mat_ambient[] = {0.625, 0.395, 0.1, 1.0};
	const GLfloat mat_diffuse[] = {0.625, 0.395, 0.1, 1.0};
	const GLfloat mat_specular[] = {0.2, 0.2, 0.2, 1.0};
	const GLfloat mat_shininess[] = {100.0};
	glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
	glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);				// Clear Screen And Depth Buffer
 	glLoadIdentity();												// Reset The Current Modelview Matrix
    gluLookAt(0,0,30.0, 0,0,0 ,0.f,1.f,0.f);
	glMultMatrixf(Transform.M);

	glPushMatrix();
	glScalef(10.0/boundingSphere.radius,10.0/boundingSphere.radius,10.0/boundingSphere.radius);
	glTranslatef(-boundingSphere.center[0],-boundingSphere.center[1],-boundingSphere.center[2]);
	
	// set transparency
	if (bTransparency)
	{
		glEnable(GL_BLEND);			// Turn Blending On (glBlendFunc(GL_SRC_ALPHA,GL_ONE) should be called)
		glDisable(GL_DEPTH_TEST);	// Turn Depth Testing Off
	}
	else
	{
		glDisable(GL_BLEND);		// Turn Blending Off
		glEnable(GL_DEPTH_TEST);	// Turn Depth Testing On
	}

	// Render the triangle mesh
	if (bShowTexture)
	{
		glEnable(GL_TEXTURE_2D);
	}

	GLfloat ref_planeS[]={texscaleX,0,0, texlocX};
	GLfloat ref_planeT[]={0,texscaleY,0, texlocY};
	GLfloat hair_planeS[]={0,0,hairscaleX, hairlocX};
	GLfloat hair_planeT[]={0,hairscaleY,0, hairlocY};
	switch(iMeshToDeform)
	{
	case 1:
		// update the control points based on the new control point positions
		pighinRBFX.UpdateFunctionValues(controlPointDisplacementX);
		pighinRBFY.UpdateFunctionValues(controlPointDisplacementY);
		pighinRBFZ.UpdateFunctionValues(controlPointDisplacementZ);
		// deform the object to render
		deformObject(pighinMeshDeformed, pighinMeshUndeformed, pighinRBFX, pighinRBFY, pighinRBFZ);
		//deformObject(pighinHairDeformed, pighinHairUndeformed, pighinRBFX, pighinRBFY, pighinRBFZ);
		deformObject(pighinFaceDeformed, pighinFaceUndeformed, pighinRBFX, pighinRBFY, pighinRBFZ);
		// make sure normal are up to date.
		if (normalUpdatesDesired)
		{
			pighinMeshDeformed->updateNormals();
			//pighinHairDeformed->updateNormals();
			pighinFaceDeformed->updateNormals();
		}

		glTexGenfv(GL_S, GL_OBJECT_PLANE, ref_planeS);		// GL_OBJECT_PLANE, GL_EYE_PLANE
		glTexGenfv(GL_T, GL_OBJECT_PLANE, ref_planeT);
		glBindTexture(GL_TEXTURE_2D,imgfile[iTexToRender]);
		pighinFaceRenderer->draw();

		glTexGenfv(GL_S, GL_OBJECT_PLANE, hair_planeS);				// GL_OBJECT_PLANE, GL_EYE_PLANE
		glTexGenfv(GL_T, GL_OBJECT_PLANE, hair_planeT);
		glBindTexture(GL_TEXTURE_2D,imgfile[4]);	// Hair image
		pighinMeshRenderer->draw();
	break;
	case 2:
		samsungRBFX.UpdateFunctionValues(controlPointDisplacementX);
		samsungRBFY.UpdateFunctionValues(controlPointDisplacementY);
		samsungRBFZ.UpdateFunctionValues(controlPointDisplacementZ);
		deformObject(samsungMeshDeformed, samsungMeshUndeformed, samsungRBFX, samsungRBFY, samsungRBFZ);
		//deformObject(samsungHairDeformed, samsungHairUndeformed, samsungRBFX, samsungRBFY, samsungRBFZ);
		deformObject(samsungFaceDeformed, samsungFaceUndeformed, samsungRBFX, samsungRBFY, samsungRBFZ);
		if (normalUpdatesDesired)
		{
			samsungMeshDeformed->updateNormals();
			//samsungHairDeformed->updateNormals();
			samsungFaceDeformed->updateNormals();
		}
		glTexGenfv(GL_S, GL_OBJECT_PLANE, ref_planeS);		// GL_OBJECT_PLANE, GL_EYE_PLANE
		glTexGenfv(GL_T, GL_OBJECT_PLANE, ref_planeT);
		glBindTexture(GL_TEXTURE_2D,imgfile[iTexToRender]);
		samsungFaceRenderer->draw();

		glTexGenfv(GL_S, GL_OBJECT_PLANE, hair_planeS);				// GL_OBJECT_PLANE, GL_EYE_PLANE
		glTexGenfv(GL_T, GL_OBJECT_PLANE, hair_planeT);
		glBindTexture(GL_TEXTURE_2D,imgfile[4]);	// Hair image
		samsungMeshRenderer->draw();
	break;
	}	
	glDisable(GL_TEXTURE_2D);	// whether we're showing the texture or not, disable texture here
	
	// Render grids in the background
	DrawGrids();

	// Render control points
	glEnable(GL_COLOR_MATERIAL);
	if (bShowCtrPts)
	{

		glColor3f(0.33f,0.33f,1.0);
		for (unsigned int i = 0; i<controlPointPosX.size(); i++ )
		{
			glPushMatrix();
			glTranslatef(controlPointPosX[i]+controlPointDisplacementX[i],controlPointPosY[i]+controlPointDisplacementY[i],controlPointPosZ[i]+controlPointDisplacementZ[i]);
			glutSolidSphere(0.02*boundingSphere.radius, 15, 15);
			glPopMatrix();
		} 
		glColor3f(0.33,1.0,0.33);
		for (unsigned int i = 0; i<controlPointPosX.size(); i++ )
		{
			glPushMatrix();
			glTranslatef(destinationPointPosX[i],destinationPointPosY[i],destinationPointPosZ[i]);
			glutSolidSphere(0.03*boundingSphere.radius, 15, 15);
			glPopMatrix();
		} 
		//glColor3f(0.33,1.0,0.33);
		glBegin(GL_LINES);
		for (unsigned int i = 0; i<controlPointPosX.size(); i++ )
		{
			glPushMatrix();
			glVertex3f(controlPointPosX[i],controlPointPosY[i],controlPointPosZ[i]);
			glVertex3f(destinationPointPosX[i],destinationPointPosY[i],destinationPointPosZ[i]);
			glPopMatrix();
		}
		glEnd();
	}
	glPopMatrix();

	glPushMatrix();
		glColor3f(1.0,0.0,0.33);
		glTranslatef(MaFace.e[0],MaFace.e[1],MaFace.e[2]);
		glutSolidSphere(0.02*boundingSphere.radius, 15, 15);
	glPopMatrix();

	glDisable(GL_COLOR_MATERIAL);	/**/

	glutSwapBuffers();
	glutPostRedisplay();
}

void loadMeshAndSetupControlPoints()
{

	// open an OBJ file to deform
	string pighinOBJ = "Data/pighin_mesh.obj";
	//string pighinHAIR = "Data/pighin_mesh_hair.obj";
	string pighinFACE = "Data/pighin_mesh_face.obj";
	string samsungOBJ = "Data/samsung_mesh.obj";
	//string samsungHAIR = "Data/samsung_mesh_hair.obj";
	string samsungFACE = "Data/samsung_mesh_face.obj";
	string WWOBJ = "Data/WW.obj";
	string SHOBJ = "Data/SH.obj";
	string JQOBJ = "Data/JQ.obj";
	pighinMeshUndeformed = new TriangleMesh(pighinOBJ);
	pighinMeshUndeformed->updateNormals();
	pighinMeshDeformed = new TriangleMesh(pighinOBJ);
	pighinMeshDeformed->updateNormals();
	//pighinHairUndeformed = new TriangleMesh(pighinHAIR);
	//pighinHairUndeformed->updateNormals();
	//pighinHairDeformed = new TriangleMesh(pighinHAIR);
	//pighinHairDeformed->updateNormals();
	pighinFaceUndeformed = new TriangleMesh(pighinFACE);
	pighinFaceUndeformed->updateNormals();
	pighinFaceDeformed = new TriangleMesh(pighinFACE);
	pighinFaceDeformed->updateNormals();
	samsungMeshUndeformed = new TriangleMesh(samsungOBJ);
	samsungMeshUndeformed->updateNormals();
	samsungMeshDeformed = new TriangleMesh(samsungOBJ);
	samsungMeshDeformed->updateNormals();
	//samsungHairUndeformed = new TriangleMesh(samsungHAIR);
	//samsungHairUndeformed->updateNormals();
	//samsungHairDeformed = new TriangleMesh(samsungHAIR);
	//samsungHairDeformed->updateNormals();
	samsungFaceUndeformed = new TriangleMesh(samsungFACE);
	samsungFaceUndeformed->updateNormals();
	samsungFaceDeformed = new TriangleMesh(samsungFACE);
	samsungFaceDeformed->updateNormals();
	WWMeshUndeformed = new TriangleMesh(WWOBJ);
	SHMeshUndeformed = new TriangleMesh(SHOBJ);
	JQMeshUndeformed = new TriangleMesh(JQOBJ);	

	boundingSphere = pighinMeshUndeformed->getBoundingSphere();
	//samsungMeshUndeformed->updateToBoundingSphere(boundingSphere);
	//samsungMeshDeformed->updateToBoundingSphere(boundingSphere);
	WWMeshUndeformed->updateToBoundingSphere(boundingSphere);
	SHMeshUndeformed->updateToBoundingSphere(boundingSphere);
	JQMeshUndeformed->updateToBoundingSphere(boundingSphere);

	// set up a renderer for the mesh 
	pighinMeshRenderer = new RenderTriangleMesh(*pighinMeshDeformed);
	//pighinHairRenderer = new RenderTriangleMesh(*pighinHairDeformed);
	pighinFaceRenderer = new RenderTriangleMesh(*pighinFaceDeformed);
	samsungMeshRenderer = new RenderTriangleMesh(*samsungMeshDeformed);
	//samsungHairRenderer = new RenderTriangleMesh(*samsungHairDeformed);
	samsungFaceRenderer = new RenderTriangleMesh(*samsungFaceDeformed);

	for (int i = 0; i<numControlPoints; i++)
	{
		Vector3 pos = pighinMeshUndeformed->getParticles()[pighinControlArray[i]].getPos();
		controlPointPosX.push_back(pos[0]);
		controlPointPosY.push_back(pos[1]);
		controlPointPosZ.push_back(pos[2]);
	}

	// allocate vectors for storing displacements
	for (unsigned int i = 0; i<controlPointPosX.size();  i++)
	{
		controlPointDisplacementX.push_back(0.0f);
		controlPointDisplacementY.push_back(0.0f);
		controlPointDisplacementZ.push_back(0.0f);
	}

	// allocate vectors for storing destinations	/**/
	for (unsigned int i = 0; i<controlPointPosX.size();  i++)
	{
		Vector3 pos = samsungMeshUndeformed->getParticles()[samsungControlArray[i]].getPos();
		destinationPointPosX.push_back(pos[0]);
		destinationPointPosY.push_back(pos[1]);
		destinationPointPosZ.push_back(pos[2]);
	}

	// initialize interpolation functions
	pighinRBFX = RBFInterpolator(controlPointPosX, controlPointPosY, controlPointPosZ, controlPointDisplacementX );
	pighinRBFY = RBFInterpolator(controlPointPosX, controlPointPosY, controlPointPosZ, controlPointDisplacementY );
	pighinRBFZ = RBFInterpolator(controlPointPosX, controlPointPosY, controlPointPosZ, controlPointDisplacementZ );
	samsungRBFX = RBFInterpolator(destinationPointPosX, destinationPointPosY, destinationPointPosZ, controlPointDisplacementX );
	samsungRBFY = RBFInterpolator(destinationPointPosX, destinationPointPosY, destinationPointPosZ, controlPointDisplacementY );
	samsungRBFZ = RBFInterpolator(destinationPointPosX, destinationPointPosY, destinationPointPosZ, controlPointDisplacementZ );
}	


int main(int argc, char** argv)
{

	glutInitDisplayMode( GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);
    glutInitWindowPosition( GL_WIN_INITIAL_X, GL_WIN_INITIAL_Y );
    glutInitWindowSize( GL_WIN_WIDTH, GL_WIN_HEIGHT );
    glutInit( &argc, argv );

	glutCreateWindow("Face Modeling");

	InitializeGlew();

	loadMeshAndSetupControlPoints();
    /*
       The function below are called when the respective event
       is triggered.
    */
    glutReshapeFunc(glutResize);       // called every time the screen is resized
    glutDisplayFunc(glutDisplay);      // called when window needs to be redisplayed
    glutKeyboardFunc(glutKeyboard);    // called when the application receives a input from the keyboard
    glutMouseFunc(glutMouse);          // called when the application receives a input from the mouse
    glutMotionFunc(glutMotion);        // called when the mouse moves over the screen with one of this button pressed
	
	/* a user menu interface */
	glutCreateMenu(glutMenu);
	glutAddMenuEntry("Use Pighin Mesh", 1);
	glutAddMenuEntry("Use Samsung Mesh", 2);
	glutAttachMenu(GLUT_RIGHT_BUTTON);

	/*
        Do lighting, material, etc initialization or
        configuration here.
    */
    InitializeOpenGL();

	/*
       Application's main loop. All the above functions
	 are called whe the respective events are triggered
    */
	glutMainLoop();
}
  