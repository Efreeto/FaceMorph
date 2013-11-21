
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
bool		bMouseDown = false;

RenderTriangleMesh* pighinMeshRenderer;								// Used for rendering a TriangleMesh.
RenderTriangleMesh* samsungMeshRenderer;

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
TriangleMesh *samsungMeshUndeformed;
TriangleMesh *samsungMeshDeformed;

BoundingSphere boundingSphere;										// Bounding sphere to give information about the size of the object.

RBFInterpolator pighinRBFX, pighinRBFY, pighinRBFZ;					// Our wonderful radial basis function interpolation things :)
RBFInterpolator samsungRBFX, samsungRBFY, samsungRBFZ;

bool normalUpdatesDesired = true;									// When the program is running press 'n' to disable normal updates
bool bShowCtrPts = true;											// Show control points for deforming the mesh. press 'b'
bool bTransparency = false;											// Make the mesh transparent. press 'm'

bool bPause = false; /* pause time */
static float morphtime = 0;	/* Mesh morphing time with cos(time) */
unsigned int iMeshToDeform = 1;	/* mesh selector 1:Pighin, 2:Samsung */

// we want control points which we place at different vertex positions
const int numControlPoints = 24;	/* can't be less than 7 */

// corresponding point indices
const int pighinControlArray[numControlPoints] = {3250,226,	3130,3977,3620,3270,  93,956,597,247,
	332,2443,2609,2411,3097,76,  4954,2360,1931,1547,2309,1506,  1833,2207};
const int samsungControlArray[numControlPoints] = {265,995,	1247,1248,100,591,  792,671,1252,1250,
	574,610,785,578,32,812,  92,1228,800,1236,1242,620,  619,607};

void InitializeGlew()
{
	GLenum err = glewInit();
	if (GLEW_OK != err)
	{
		/* Problem: glewInit failed, something is seriously wrong. */
		fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
	}

	if (!(GLEW_ARB_vertex_buffer_object	&& GLEW_EXT_framebuffer_object))
	{
		cerr << "Extensions not supported" << endl;
		exit(1);
	}
}

void DrawNet(GLfloat size, GLint LinesX, GLint LinesZ)
{
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

void InitializeOpenGL()					// Any GL Init Code & User Initialization Goes Here
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(45,640/480,1,-100);
	glMatrixMode(GL_MODELVIEW);

	bMouseDown	= false;

	glClearColor(0.2,0.2,0.2,1.0);
	//glClearDepth (1.0f);			
	//glDepthFunc (GL_LEQUAL);
	glEnable (GL_DEPTH_TEST);		
	glShadeModel (GL_SMOOTH);										

	//glEnable(GL_AUTO_NORMAL);
	glEnable(GL_NORMALIZE);

	glBlendFunc(GL_SRC_ALPHA,GL_ONE);		// Blending Function For Translucency Based On Source Alpha Value ( NEW )

	GLfloat light_ambient[] = {0.2, 0.2, 0.2, 1.0};
	GLfloat light_diffuse[] = {0.6, 0.6, 0.6, 1.0};
	GLfloat light_specular[] = {0.6, 0.6, 0.6, 1.0};
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);

	glEnable(GL_LIGHT0);	
	glEnable(GL_LIGHTING);	

	glPointSize(3);
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
		case 'p':
			bPause=!bPause;
			break;
		case 's':
			bPause=true;
			if (morphtime!=M_PI)
				morphtime = M_PI;
			else
				morphtime = 0;
			break;
		case 'r':
			bPause=true;
			morphtime=0;
			Matrix3fSetIdentity(&LastRot);	// three lines here reset all rotations
			Matrix3fSetIdentity(&ThisRot);								
			Matrix4fSetRotationScaleFromMatrix3f(&Transform, &ThisRot);	
			break;
		case ' ':	// When the Space Bar is pressed
			bShowCtrPts=!bShowCtrPts;
			break;
		case 'm':
			bTransparency=!bTransparency;
			break;
		break;

		case 'z':
			Matrix4fMulRotationScale(&Transform, 1.2);
			break;
		case 'x':
			Matrix4fMulRotationScale(&Transform, 0.8);
			break;
		case 'c':
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			break;
		case 'g':
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			break;
		case 'v':
			glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
			break;
		case 'n': 
			normalUpdatesDesired = !normalUpdatesDesired;
			break;
		case 27: // ESC
			exit (1);
			break;
    }

	glutPostRedisplay();
}

void glutMouse(int button, int state, int x, int y)
{
	MousePt.s.X=x;
	MousePt.s.Y=y;

	if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN)
	{
		LastRot = ThisRot;										// Set Last Static Rotation To Last Dynamic One
		ArcBall.click(&MousePt);								// Update Start Vector And Prepare For Dragging
		bMouseDown = true;
	}
	else
	{
		bMouseDown = false;
	}
}

void glutMotion(int x, int y)									// Perform Motion Updates Here
{
	MousePt.s.X=x;
	MousePt.s.Y=y;

	Quat4fT     ThisQuat;

	if(bMouseDown == true)
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
		for (unsigned int i = 0; i<controlPointPosX.size(); i++)
		{
			iMeshToDeform=1;
			Vector3 src = pighinMeshUndeformed->getParticles()[pighinControlArray[i]].getPos();
			controlPointPosX[i]=src[0];
			controlPointPosY[i]=src[1];
			controlPointPosZ[i]=src[2];

			Vector3 dst = samsungMeshUndeformed->getParticles()[samsungControlArray[i]].getPos();
			destinationPointPosX[i]=dst[0];
			destinationPointPosY[i]=dst[1];
			destinationPointPosZ[i]=dst[2];
		}
	}
	else if (entryID==2)
	{
		for (unsigned int i = 0; i<controlPointPosX.size(); i++)
		{
			iMeshToDeform=2;
			Vector3 src = samsungMeshUndeformed->getParticles()[samsungControlArray[i]].getPos();
			controlPointPosX[i]=src[0];
			controlPointPosY[i]=src[1];
			controlPointPosZ[i]=src[2];
			
			Vector3 dst = pighinMeshUndeformed->getParticles()[pighinControlArray[i]].getPos();
			destinationPointPosX[i]=dst[0];
			destinationPointPosY[i]=dst[1];
			destinationPointPosZ[i]=dst[2];
		}
	}
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
	if (!bPause) {
		morphtime = morphtime + 0.05;
	}
	
	// move control points
	for (unsigned int i = 0; i<controlPointPosX.size(); i++ )
	{
		controlPointDisplacementX[i] = (destinationPointPosX[i]-controlPointPosX[i])*(-cosf(morphtime)/2+0.5);	// -cos()/2+1/2 go from 0 to 1
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
	
	switch(iMeshToDeform)
	{
	case 1:
		// update the control points based on the new control point positions
		pighinRBFX.UpdateFunctionValues(controlPointDisplacementX);
		pighinRBFY.UpdateFunctionValues(controlPointDisplacementY);
		pighinRBFZ.UpdateFunctionValues(controlPointDisplacementZ);
		// deform the object to render
		deformObject(pighinMeshDeformed, pighinMeshUndeformed, pighinRBFX, pighinRBFY, pighinRBFZ);
		// make sure normal are up to date.
		if (normalUpdatesDesired)
		{
			pighinMeshDeformed->updateNormals();
		}		
		pighinMeshRenderer->draw();
		break;
	case 2:
		// update the control points based on the new control point positions
		samsungRBFX.UpdateFunctionValues(controlPointDisplacementX);
		samsungRBFY.UpdateFunctionValues(controlPointDisplacementY);
		samsungRBFZ.UpdateFunctionValues(controlPointDisplacementZ);
		// deform the object to render
		deformObject(samsungMeshDeformed, samsungMeshUndeformed, samsungRBFX, samsungRBFY, samsungRBFZ);
		// make sure normal are up to date.
		if (normalUpdatesDesired)
		{
			samsungMeshDeformed->updateNormals();
		}
		samsungMeshRenderer->draw();
		break;
	}

	DrawGrids();	// Render grids in the background

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

	glDisable(GL_COLOR_MATERIAL);

	glutSwapBuffers();
	glutPostRedisplay();
}

void loadMeshAndSetupControlPoints()
{

	// open an OBJ file to deform
	string pighinOBJ = "Data/pighin_mesh.obj";
	string samsungOBJ = "Data/samsung_mesh.obj";
	pighinMeshUndeformed = new TriangleMesh(pighinOBJ);
	pighinMeshUndeformed->updateNormals();
	pighinMeshDeformed = new TriangleMesh(pighinOBJ);
	pighinMeshDeformed->updateNormals();
	samsungMeshUndeformed = new TriangleMesh(samsungOBJ);
	samsungMeshUndeformed->updateNormals();
	samsungMeshDeformed = new TriangleMesh(samsungOBJ);
	samsungMeshDeformed->updateNormals();

	boundingSphere = pighinMeshUndeformed->getBoundingSphere();

	// set up a renderer for the mesh 
	pighinMeshRenderer = new RenderTriangleMesh(*pighinMeshDeformed);
	samsungMeshRenderer = new RenderTriangleMesh(*samsungMeshDeformed);

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
  