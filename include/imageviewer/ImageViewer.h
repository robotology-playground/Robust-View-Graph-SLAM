/*
	A basic image viewer
	Copyright (C) 2012  www.scratchapixel.com

	This program is free software: you can redistribute it and/or modify
	it under the terms of the GNU General Public License as published by
	the Free Software Foundation, either version 3 of the License, or
	(at your option) any later version.

	This program is distributed in the hope that it will be useful,
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
	GNU General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <iostream>
#include <cmath>

#include <cstring>
#include <cassert>

#include <Cg/cgGL.h>
#include <GL/glfw.h>

/*
** 
** GLOBALS
**
*/

GLuint texID;
bool mouse[3] = { false };
bool update = false;
bool exposureKey = false;
int lastPosX = 0;

float exposure = 0; 	// Current exposure setting.  All pixels
						// are multiplied by pow(2,exposure) before
						// they appear on the screen

CGcontext cgContext;
CGprogram cgProgram;
CGprofile cgProfile;

/*
** 
** OpenGL related functions
**
*/

void handleCgErrors () {
	std::cerr << cgGetErrorString (cgGetError()) << std::endl;
	std::cerr << cgGetLastListing (cgContext) << std::endl;
	exit (1);
}

void checkGlErrors(const char where[]) {
	GLenum error = glGetError();
	if (error != GL_NO_ERROR) {
		//std::cerr << where << ": " << gluErrorString (error) << std::endl;
		exit (1);
	}
}

//
// Shader for RGB images
//

const char shaderRgbSource[] =
"                                                                          \n"
"    struct Out                                                            \n"
"    {                                                                     \n"
"        half3 pixel: COLOR;                                               \n"
"    };                                                                    \n"
"                                                                          \n"
"    Out                                                                   \n"
"    main (float2 tc: TEXCOORD0,                                           \n"
"          uniform sampler2D inputImage: TEXUNIT0,                         \n"
"          uniform float expMult)                                          \n"
"    {                                                                     \n"
"        //                                                                \n"
"        // Apply exposure                                                 \n"
"        //                                                                \n"
"                                                                          \n"
"        half3 color = tex2D(inputImage, tc).rgb * expMult;                \n"
"                                                                          \n"
"                                                                          \n"
"        //                                                                \n"
"        // Assign fragment color                                          \n"
"        //                                                                \n"
"                                                                          \n"
"        Out output;                                                       \n"
"        output.pixel = color;                                             \n"
"        return output;                                                    \n"
"    }                                                                     \n"
"                                                                          \n";

void initShaderRgb() {
	cgSetErrorCallback (handleCgErrors);

	cgContext = cgCreateContext();
	cgProfile = cgGLGetLatestProfile(CG_GL_FRAGMENT);
	cgGLSetOptimalOptions(cgProfile);

	cgProgram = cgCreateProgram (cgContext, CG_SOURCE, shaderRgbSource, cgProfile, "main", 0);

	cgGLLoadProgram(cgProgram);
	cgGLBindProgram(cgProgram);
	cgGLEnableProfile(cgProfile);

	CGparameter emParam = cgGetNamedParameter (cgProgram, "expMult");
	cgSetParameter1f(emParam, pow (2.0f, exposure));
}

void initRgbTexture(const uint32_t &w, const uint32_t &h, const float *rawPixels) {
	glBindTexture(GL_TEXTURE_2D, texID);
	glPixelStorei (GL_UNPACK_ALIGNMENT, 1);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_FLOAT, rawPixels);

	checkGlErrors("initRgbTexture");
}

void draw(const unsigned char *screenPixels) {
	int w, h;
	glfwGetWindowSize(&w, &h);
    
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	glViewport (0, 0, w, h); // optional
	glOrtho(0, w, 0, h, -1, 1);
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity();

	glClearColor (0.3, 0.3, 0.3, 1.0);
	glClear(GL_COLOR_BUFFER_BIT);
	glActiveTexture(GL_TEXTURE0);
	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, texID);
	cgGLEnableProfile (cgProfile);

	glBegin(GL_POLYGON);
	glTexCoord2f(0, 1); glVertex2f(0, 0);
	glTexCoord2f(1, 1); glVertex2f(w, 0);
	glTexCoord2f(1, 0); glVertex2f(w, h);
	glTexCoord2f(0, 0); glVertex2f(0, h);
	glEnd();

	glDisable(GL_TEXTURE_2D);

	checkGlErrors("draw");
}

void openImageFile( const char *filename, uint32_t &w, uint32_t &h, 
	unsigned char *&screenPixels, float *&rawPixels) {
	w = 640, h = 480; // default size
	std::ifstream ifs(filename);
	if (ifs.good()) {
		std::string header;
		ifs >> header; // P6
		assert(strcmp(header.c_str(), "P6") == 0);
		uint32_t maxcol;
		ifs >> w >> h >> maxcol;
	}
	else {
		fprintf(stderr, "Can't open file %s\n", filename);
		ifs.close();
	}

	screenPixels = new unsigned char[w * h * 3];
	rawPixels = new float[w * h * 3];

	if (ifs.good()) {
		ifs.ignore();
		ifs.read((char*)screenPixels, w * h * 3);
	}
	else {
		unsigned char *p = screenPixels;
		for (unsigned j = 0; j < h; ++j) {
			for (unsigned i = 0; i < w; ++i) {
				*p = *(p + 1) = *(p + 2) = ((i & 32) ^ (j & 32)) ? 80 : 150;
				p += 3;
			}
		}
	}

	for (uint32_t i = 0; i < w * h * 3; ++i)
		rawPixels[i] = screenPixels[i] / 255.f;
	ifs.close();	
}

/*
** 
** GLFW callback functions
**
*/

void keyboardCB(int key, int action) {
	switch (key) {
		case 82: // r(eset)
			if (action) exposure = 0, update = true;
			break;
		case 69: // e(xposure)
			exposureKey = (action == 1) ? true : false;
			break;
		default:
			break;
	}
}

void mousePosCB(int x, int y) {
	if (exposureKey && mouse[GLFW_MOUSE_BUTTON_LEFT]) {
		int dx = x - lastPosX;
		exposure +=  dx / 1000.f;
		// update image data and update texture
		CGparameter emParam = cgGetNamedParameter(cgProgram, "expMult");
		cgSetParameter1f(emParam, pow (2.0f, exposure));
	}
	lastPosX = x;
}

void mouseButtonCB(int id, int action) {
	switch (id) {
		case GLFW_MOUSE_BUTTON_LEFT:
			mouse[GLFW_MOUSE_BUTTON_LEFT] = (action == GLFW_PRESS) ? true : false;
			break;
		case GLFW_MOUSE_BUTTON_MIDDLE:
			mouse[GLFW_MOUSE_BUTTON_MIDDLE] = (action == GLFW_PRESS) ? true : false;
			break;
		case GLFW_MOUSE_BUTTON_RIGHT:
			mouse[GLFW_MOUSE_BUTTON_RIGHT] = (action == GLFW_PRESS) ? true : false;
			break;
		default: break;
	}	
}

void ImageViewer(const char *filename_1, const char *filename_2, uint32_t height, uint32_t width) {

	unsigned char *screenPixels_1 = NULL;
    unsigned char *screenPixels_2 = NULL;
	float *rawPixels_1 = NULL;
    float *rawPixels_2 = NULL;
	openImageFile (filename_1, width, height, screenPixels_1, rawPixels_1) ;
    openImageFile (filename_2, width, height, screenPixels_2, rawPixels_2) ;
    
    //GLFWwindow* window_1, window_2;
    //window_1 = glfwCreateWindow(width, height, "left", NULL, NULL);
    //window_2 = glfwCreateWindow(width, height, "right", NULL, NULL);

	uint32_t redBits = 8, greenBits = 8, blueBits = 8;
	uint32_t alphaBits = 8, depthBits = 32, stencilBits = 0;

	glfwInit () ;
	glfwOpenWindowHint (GLFW_WINDOW_NO_RESIZE, GL_TRUE) ;
	glfwOpenWindow (width, height, redBits, greenBits, blueBits, alphaBits, depthBits, stencilBits, GLFW_WINDOW) ;
    
    // Mouse and keyboard stuff
	glfwSetKeyCallback (&keyboardCB) ;
	glfwSetMousePosCallback (&mousePosCB) ;
	glfwSetMouseButtonCallback (&mouseButtonCB) ;
	glfwEnable (GLFW_KEY_REPEAT) ;
	glfwEnable (GLFW_STICKY_MOUSE_BUTTONS) ;
	
	// init Cg shader
	initShaderRgb () ;

	// init texture
	glGenTextures (1, &texID) ;

	// you need an OpenGL context before you can create a texture
	//initRgbTexture (width, height, rawPixels_1) ;
    initRgbTexture (width, height, rawPixels_2) ;

	while(!glfwGetKey(GLFW_KEY_ESC) && glfwGetWindowParam(GLFW_OPENED)) {
		glClear (GL_COLOR_BUFFER_BIT) ;
        
        //glfwInit ( );
        //glfwOpenWindowHint (GLFW_WINDOW_NO_RESIZE, GL_TRUE) ;
        //glfwOpenWindow (width, height, redBits, greenBits, blueBits, alphaBits, depthBits, stencilBits, GLFW_WINDOW) ;
        
        //glfwMakeContextCurrent (window_1) ;
        //glfwSetWindowTitle ("Simple Image Viewer 1") ;
		//draw (screenPixels_1) ;
        
        //glfwInit () ;
        //glfwOpenWindowHint (GLFW_WINDOW_NO_RESIZE, GL_TRUE) ;
        //glfwOpenWindow (width, height, redBits, greenBits, blueBits, alphaBits, depthBits, stencilBits, GLFW_WINDOW) ;
        
        //glfwMakeContextCurrent (window_2) ;
        glfwSetWindowTitle ("Simple Image Viewer 2") ;
        draw (screenPixels_2) ;
        
		glfwSwapBuffers () ;
		glfwWaitEvents () ;
	}

	glfwCloseWindow () ;
	glfwTerminate () ;
	delete [] rawPixels_1 ;
	delete [] screenPixels_1 ;
    delete [] rawPixels_2 ;
	delete [] screenPixels_2 ;
	//exit(EXIT_SUCCESS) ;
}
