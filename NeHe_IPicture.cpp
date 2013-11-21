/****************************************************************************
*                                     *                                     *
*  Jeff Molofee's IPicture Basecode   *    Huge Thanks To: Dave Richards    *
*       http://nehe.gamedev.net       *                    Bobby Ward &     *
*                2002                 *                    The MSDN         *
*                                     *                                     *
*****************************************************************************
*                                                                           *
*   Loads  : BMP, EMF, GIF, ICO, JPG, WMF                                   *
*   Source : Reads From Disk Or The Internet                                *
*   Extras : Images Can Be Any Width Or Height                              *
*                                                                           *
*****************************************************************************
*                                                                           *
*   ReshapeGL  : Set Your Aspect Ratio, How You Want                        *
*   WindowProc : Add Custom WM_ Events (Mouse, Etc)                         *
*   WinMain    : Set The Window Title                                       *
*                Set Resolution & Color Depth                               *
*                Remove 4 Lines Of Code To Force Fullscreen (Commented)     *
*                                                                           *
*****************************************************************************
*                                                                           *
*   Free To Use In Projects Of Your Own.  All I Ask For Is A Simple Greet   *
*   Or Mention Of The Site In Your Readme Or The Project Itself :)          *
*                                                                           *
****************************************************************************/

#include <windows.h>													// Header File For Windows
#include <gl\gl.h>														// Header File For The OpenGL32 Library
#include <gl\glu.h>														// Header File For The GLu32 Library
#include <olectl.h>														// Header File For The OLE Controls Library	(Used In BuildTexture)
#include <math.h>														// Header File For The Math Library			(Used In BuildTexture)

//#include "NeHeGL.h"														// Header File For NeHeGL

#pragma comment( lib, "opengl32.lib" )									// Search For OpenGL32.lib While Linking
#pragma comment( lib, "glu32.lib" )										// Search For GLu32.lib While Linking

int BuildTexture(char *szPathName, GLuint &texid)						// Load Image And Convert To A Texture
{
	HDC			hdcTemp;												// The DC To Hold Our Bitmap
	HBITMAP		hbmpTemp;												// Holds The Bitmap Temporarily
	IPicture	*pPicture;												// IPicture Interface
	OLECHAR		wszPath[MAX_PATH+1];									// Full Path To Picture (WCHAR)
	char		szPath[MAX_PATH+1];										// Full Path To Picture
	long		lWidth;													// Width In Logical Units
	long		lHeight;												// Height In Logical Units
	long		lWidthPixels;											// Width In Pixels
	long		lHeightPixels;											// Height In Pixels
	GLint		glMaxTexDim ;											// Holds Maximum Texture Size

	if (strstr(szPathName, "http://"))									// If PathName Contains http:// Then...
	{
		strcpy(szPath, szPathName);										// Append The PathName To szPath
	}
	else																// Otherwise... We Are Loading From A File
	{
		GetCurrentDirectoryA(MAX_PATH, (LPSTR)szPath);							// Get Our Working Directory
		strcat(szPath, "\\");											// Append "\" After The Working Directory
		strcat(szPath, szPathName);										// Append The PathName
	}

	MultiByteToWideChar(CP_ACP, 0, szPath, -1, wszPath, MAX_PATH);		// Convert From ASCII To Unicode
	HRESULT hr = OleLoadPicturePath(wszPath, 0, 0, 0, IID_IPicture, (void**)&pPicture);

	if(FAILED(hr))														// If Loading Failed
		return FALSE;													// Return False

	hdcTemp = CreateCompatibleDC(GetDC(0));								// Create The Windows Compatible Device Context
	if(!hdcTemp)														// Did Creation Fail?
	{
		pPicture->Release();											// Decrements IPicture Reference Count
		return FALSE;													// Return False (Failure)
	}

	glGetIntegerv(GL_MAX_TEXTURE_SIZE, &glMaxTexDim);					// Get Maximum Texture Size Supported
	
	pPicture->get_Width(&lWidth);										// Get IPicture Width (Convert To Pixels)
	lWidthPixels	= MulDiv(lWidth, GetDeviceCaps(hdcTemp, LOGPIXELSX), 2540);
	pPicture->get_Height(&lHeight);										// Get IPicture Height (Convert To Pixels)
	lHeightPixels	= MulDiv(lHeight, GetDeviceCaps(hdcTemp, LOGPIXELSY), 2540);

	// Resize Image To Closest Power Of Two
	if (lWidthPixels <= glMaxTexDim) // Is Image Width Less Than Or Equal To Cards Limit
		lWidthPixels = 1 << (int)floor((log((double)lWidthPixels)/log(2.0f)) + 0.5f); 
	else  // Otherwise  Set Width To "Max Power Of Two" That The Card Can Handle
		lWidthPixels = glMaxTexDim;
 
	if (lHeightPixels <= glMaxTexDim) // Is Image Height Greater Than Cards Limit
		lHeightPixels = 1 << (int)floor((log((double)lHeightPixels)/log(2.0f)) + 0.5f);
	else  // Otherwise  Set Height To "Max Power Of Two" That The Card Can Handle
		lHeightPixels = glMaxTexDim;
	
	//	Create A Temporary Bitmap
	BITMAPINFO	bi = {0};												// The Type Of Bitmap We Request
	DWORD		*pBits = 0;												// Pointer To The Bitmap Bits

	bi.bmiHeader.biSize			= sizeof(BITMAPINFOHEADER);				// Set Structure Size
	bi.bmiHeader.biBitCount		= 32;									// 32 Bit
	bi.bmiHeader.biWidth		= lWidthPixels;							// Power Of Two Width
	bi.bmiHeader.biHeight		= lHeightPixels;						// Make Image Top Up (Positive Y-Axis)
	bi.bmiHeader.biCompression	= BI_RGB;								// RGB Encoding
	bi.bmiHeader.biPlanes		= 1;									// 1 Bitplane

	//	Creating A Bitmap This Way Allows Us To Specify Color Depth And Gives Us Imediate Access To The Bits
	hbmpTemp = CreateDIBSection(hdcTemp, &bi, DIB_RGB_COLORS, (void**)&pBits, 0, 0);
	
	if(!hbmpTemp)														// Did Creation Fail?
	{
		DeleteDC(hdcTemp);												// Delete The Device Context
		pPicture->Release();											// Decrements IPicture Reference Count
		return FALSE;													// Return False (Failure)
	}

	SelectObject(hdcTemp, hbmpTemp);									// Select Handle To Our Temp DC And Our Temp Bitmap Object

	// Render The IPicture On To The Bitmap
	pPicture->Render(hdcTemp, 0, 0, lWidthPixels, lHeightPixels, 0, lHeight, lWidth, -lHeight, 0);

	// Convert From BGR To RGB Format And Add An Alpha Value Of 255
	for(long i = 0; i < lWidthPixels * lHeightPixels; i++)				// Loop Through All Of The Pixels
	{
		BYTE* pPixel	= (BYTE*)(&pBits[i]);							// Grab The Current Pixel
		BYTE  temp		= pPixel[0];									// Store 1st Color In Temp Variable (Blue)
		pPixel[0]		= pPixel[2];									// Move Red Value To Correct Position (1st)
		pPixel[2]		= temp;											// Move Temp Value To Correct Blue Position (3rd)

		// This Will Make Any Black Pixels, Completely Transparent		(You Can Hardcode The Value If You Wish)
		if ((pPixel[0]==0) && (pPixel[1]==0) && (pPixel[2]==0))			// Is Pixel Completely Black
			pPixel[3]	=   0;											// Set The Alpha Value To 0
		else															// Otherwise
			pPixel[3]	= 255;											// Set The Alpha Value To 255
	}

	glGenTextures(1, &texid);											// Create The Texture

	// Typical Texture Generation Using Data From The Bitmap
	glBindTexture(GL_TEXTURE_2D, texid);								// Bind To The Texture ID
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);		// (Modify This For The Type Of Filtering You Want)
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);     // (Modify This For The Type Of Filtering You Want)
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, lWidthPixels, lHeightPixels, 0, GL_RGBA, GL_UNSIGNED_BYTE, pBits);	// (Modify This If You Want Mipmaps)
	
	DeleteObject(hbmpTemp);												// Delete The Object
	DeleteDC(hdcTemp);													// Delete The Device Context

	pPicture->Release();												// Decrements IPicture Reference Count

	return TRUE;														// Return True (All Good)
}
//
//BOOL Initialize (GL_Window* window, Keys* keys)							// Any GL Init Code & User Initialiazation Goes Here
//{
//	g_window	= window;												// Window Values
//	g_keys		= keys;													// Key Values
//
//	// Start Of User Initialization
//
//	// Load .BMP From A File		(1st Texture)
//	if (!BuildTexture("Data/NeHe.bmp", texture[0]))																// (Example Code... Can Be Removed)
//		return FALSE;													// Return False If Loading Failed		// (Example Code... Can Be Removed)
//
//	// Load .JPG From A URL			(2nd Texture)																// (Example Code... Can Be Removed)
//	if (!BuildTexture("http://nehe.gamedev.net/data/downloads/o/opalis.jpg", texture[1]))						// (Example Code... Can Be Removed)
//		return FALSE;													// Return False If Loading Failed		// (Example Code... Can Be Removed)
//
//	// Load .GIF From A File		(3rd Texture)																// (Example Code... Can Be Removed)
//	if (!BuildTexture("Data/lady.gif", texture[2]))																// (Example Code... Can Be Removed)
//		return FALSE;													// Return False If Loading Failed		// (Example Code... Can Be Removed)
//
//	glEnable(GL_TEXTURE_2D);											// Enable Texture Mapping				// (Example Code... Can Be Removed)
//	glClearColor (0.0f, 0.0f, 0.0f, 0.5f);								// Black Background						// (Set To Any Color You Wish)
//	glClearDepth (1.0f);												// Depth Buffer Setup
//	glDepthFunc (GL_LEQUAL);											// The Type Of Depth Testing			// (Select The Depth Testing You Want)
//	glEnable (GL_DEPTH_TEST);											// Enable Depth Testing
//	glShadeModel (GL_SMOOTH);											// Select Smooth Shading				// (Set To Flat Shading If You Wish)
//	glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);					// Set Perspective Calculations To Most Accurate
//
//	glAlphaFunc(GL_GREATER,0.1f);										// Set Alpha Testing (To Make Black Transparent)
//	glEnable(GL_ALPHA_TEST);											// Enable Alpha Testing (To Make Black Transparent)
//	return TRUE;														// Return TRUE (Initialization Successful)
//}