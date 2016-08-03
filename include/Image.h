#ifndef IMAGE_H
#define IMAGE_H

/*This basic image class will only support floating point images
 * (a must these days) composed of three channels only. You can easily
 * extend this class as an exercise, for example you can use template to
 * define the image data type and the number of channels on the fly..
 *
 * The image resolution defines the dimension of this 2D array.
 * The width (horizontal dimension) and the height (vertical dimension)
 * define the overall resolution of the image. A pixel in the digital world
 * is generally defined by three digital values, one for each elementary or
 * primary color.
 *
 * The number of bits (not bytes) used per pixel is often abbreviated with
 * the acronym bpp; it stands for bits per pixel. Truecolor uses 24 bpp.
 * Sometimes the number of bits used is defined per channel (i.e. per color).
 * In this case the acronym is bpc (bits per channel) but it is not often used.
 * The number of bits used for each color component of a single pixel,
 * is known as the color depth or bit depth.
 *
 * Storing pixel color as floats was originally motivated by the need for
 * high dynamic range images (or HDR images), images capable of storing the
 * wider range of colors and light levels one can observe in the physical
 * world.
 *
 * From a programming point of view, all we need to do, is to create a
 * simple Image class, in which we will store the width and the height of
 * the image, as well as a 1D array of pixels (you can use a two dimensional
 * arrays if you want, but this is not necessarily more practical and is
 * not as efficient). We need two constructors, one to create an empty
 * image and one to create an image with a user specified resolution (line 20)
 * which can also be filled with a given constant color (this color default
 * to black).
 *
 * The pixels themselves are represented with a special structure named Rgb
 * (line 5). This is very convenient because we can write operators to
 * manipulate the three floats of a pixel at once. It is trivial to access
 * pixels from the image by overloading the bracket operator []
 * (line 25 and 26) and applying every mathematical operations we want on
 * these pixels: multiplying them (line 11), comparing them (line 10),
 * adding them (line 12), etc. The function (line 13) is interesting.
 * It can be used to compute the brightness of a pixel and accumulate
 * the result into a float. It will be used in the next lesson.
 *
 * Finally pixels are stored in the image as a 1D array of Rgb elements
 * (line 29). The size of this area is equal to the width of the image
 * multiplied by the image height. It is not necessary to multiply this
 * number by the number of channels since each pixel (of type Rgb) already
 * holds the memory to store three colors.
 *
 * References:
 * http://www.scratchapixel.com/
 * http://www.scratchapixel.com/old
 */

#include <string.h>
#include <fstream>
#include "GL/gl.h"

class Image
{
public:
    // Rgb structure, i.e. a pixel
    struct Rgb
    {
        Rgb() : r(0), g(0), b(0)  {}
        Rgb(float c) : r(c), g(c), b(c) {}
        Rgb(float _r, float _g, float _b) : r(_r), g(_g), b(_b) {}
        bool operator != (const Rgb &c) const { return c.r != r || c.g != g || c.b != b; }
        Rgb& operator *= (const Rgb &rgb) { r *= rgb.r, g *= rgb.g, b *= rgb.b; return *this; }
        Rgb& operator += (const Rgb &rgb) { r += rgb.r, g += rgb.g, b += rgb.b; return *this; }
        friend float& operator += (float &f, const Rgb rgb)
        { f += (rgb.r + rgb.g + rgb.b) / 3.f; return f; }
        float r, g, b;
    };
    
    Image() : w(0), h(0), pixels(NULL)
    { /* empty image */ }
    Image(const unsigned int &_w, const unsigned int &_h, const Rgb &c = kBlack) : w(_w), h(_h), pixels(NULL)
    {
        pixels = new Rgb[w * h];
        for (int i = 0; i < w * h; ++i) pixels[i] = c;
    }
    const Rgb& operator [] (const unsigned int &i) const { return pixels[i]; }
    Rgb& operator [] (const unsigned int &i) { return pixels[i]; }
    ~Image() { if (pixels != NULL) delete [] pixels; }
    unsigned int w, h; // image resolution
    Rgb *pixels; // 1D array of pixels
    static const Rgb kBlack, kWhite, kRed, kGreen, kBlue; // preset colors
};

const Image::Rgb Image::kBlack = Image::Rgb(0);
const Image::Rgb Image::kWhite = Image::Rgb(1);
const Image::Rgb Image::kRed = Image::Rgb(1,0,0);
const Image::Rgb Image::kGreen = Image::Rgb(0,1,0);
const Image::Rgb Image::kBlue = Image::Rgb(0,0,1);//?

// /*
// ** 
// ** GLOBALS
// **
// */
// 
// GLuint texID;
// bool mouse[3] = { false };
// bool update = false;
// bool exposureKey = false;
// int lastPosX = 0;
// 
// float exposure = 0; 	// Current exposure setting.  All pixels
// 						// are multiplied by pow(2,exposure) before
// 						// they appear on the screen
// 
// CGcontext cgContext;
// CGprogram cgProgram;
// CGprofile cgProfile;
// 
// /*
// ** 
// ** OpenGL related functions
// **
// */
// 
// void handleCgErrors ()
// {
//     std::cerr << cgGetErrorString (cgGetError()) << std::endl;
//     std::cerr << cgGetLastListing (cgContext) << std::endl;
//     exit (1);
// }
// 
// void checkGlErrors(const char where[])
// {
//     GLenum error = glGetError();
//     
//     if (error != GL_NO_ERROR) {
//         //std::cerr << where << ": " << gluErrorString (error) << std::endl;
//         exit (1);
//     }
// }
// 
// //
// // Shader for RGB images
// //
// /* The exposure of the displayed can be changed image by holding down the 'e'
//  * key and moving the mouse (left button) to the right to brighten the image
//  * or to the left to make it darker (effectively changing the exposure of
//  * the image). Rather than changing the pixels of the image itself and
//  * resetting the OpenGL texture every time the exposure changes, we can use
//  * a fragment shader instead. The fragment shader is attached to the rectangular
//  * polygon. Its arguments are the texture, the texture coordinates for the
//  * current fragment being rendered and the exposure value. The color from
//  * the texture for the current texture coordinates is retrieved (line 17)
//  * and multiplied by a number which is computed with the following formula:
//  */
// const char shaderRgbSource[] =
//         "                                                                          \n"
//         "    struct Out                                                            \n"
//         "    {                                                                     \n"
//         "        half3 pixel: COLOR;                                               \n"
//         "    };                                                                    \n"
//         "                                                                          \n"
//         "    Out                                                                   \n"
//         "    main (float2 tc: TEXCOORD0,                                           \n"
//         "          uniform sampler2D inputImage: TEXUNIT0,                         \n"
//         "          uniform float expMult)                                          \n"
//         "    {                                                                     \n"
//         "        //                                                                \n"
//         "        // Apply exposure                                                 \n"
//         "        //                                                                \n"
//         "                                                                          \n"
//         "        half3 color = tex2D(inputImage, tc).rgb * expMult;                \n"
//         "                                                                          \n"
//         "                                                                          \n"
//         "        //                                                                \n"
//         "        // Assign fragment color                                          \n"
//         "        //                                                                \n"
//         "                                                                          \n"
//         "        Out output;                                                       \n"
//         "        output.pixel = color;                                             \n"
//         "        return output;                                                    \n"
//         "    }                                                                     \n"
//         "                                                                          \n";
// 
// void initShaderRgb()
// {
//     cgSetErrorCallback (handleCgErrors);
//     
//     cgContext = cgCreateContext();
//     cgProfile = cgGLGetLatestProfile(CG_GL_FRAGMENT);
//     cgGLSetOptimalOptions(cgProfile);
//     
//     cgProgram = cgCreateProgram(cgContext, CG_SOURCE, shaderRgbSource, cgProfile, "main", 0);
//     
//     cgGLLoadProgram(cgProgram);
//     cgGLBindProgram(cgProgram);
//     cgGLEnableProfile(cgProfile);
//     
//     CGparameter emParam = cgGetNamedParameter (cgProgram, "expMult");
//     cgSetParameter1f(emParam, pow (2.0f, exposure));
// }

/* readPPM :
 *
 * The header (the few lines at the top of the file providing the necessary
 * information about the image we are about to read) is defined by three lines.
 * The first line contains a string which can either be "P1", "P2", ... "P6".
 * We are only concerned by "P6" PPM images which stands for binary color
 * images. The next line contains two integers representing the resolution
 * of the image. The third line contains an another integer representing
 * the number of levels each color is encoded with (i.e. the number of grey
 * values between black and white, generally 255). You can insert comments
 * in the header; they start with the letter '#', but we will ignore them
 * in this lesson.
 *
 * In the PPM format, colors are encoded using 24 bpp. In other worlds each
 * color is encoded as 1 bytes (in C++ an unsigned char). Each pixel is
 * encoded in the file as a sequence of three bytes, one for each
 * color.
 *
 * Pixels are encoded in scanline order, from the first row to the last row
 * in the image, for left to right. In other words the first three
 * bytes in the file are for the pixel at location (1,1) in the image
 * (where (1,1) denotes the upper left corner of the image), and the last
 * three bytes in the file represent the lower right pixel in the frame
 * (location (640, 480) assuming the image has resolution 640x480).
 */
Image readPPM(const char *filename)
{
    std::cout<<"Reading image: "<<filename<<std::endl ;
    std::ifstream ifs;
    ifs.open(filename, std::ios::binary); // need to spec. binary mode for Windows users
    Image img;
    try {
        if (ifs.fail()) { throw("Can't open input file"); }
        std::string header;
        int w, h, b;
        ifs >> header;
        if (strcmp(header.c_str(), "P6") != 0) throw("Can't read input file");
        ifs >> w >> h >> b;
        img.w = w; img.h = h;
        img.pixels = new Image::Rgb[w * h]; // this is throw an exception if bad_alloc
        ifs.ignore(256, '\n'); // skip empty lines in necessary until we get to the binary data
        unsigned char pix[3];
        // read each pixel one by one and convert bytes to floats
        for (int i = 0; i < w * h; ++i) {
            ifs.read(reinterpret_cast<char *>(pix), 3);
            img.pixels[i].r = pix[0] / 255.f;
            img.pixels[i].g = pix[1] / 255.f;
            img.pixels[i].b = pix[2] / 255.f;
        }
        ifs.close();
    }
    catch (const char *err) {
        fprintf(stderr, "%s\n", err);
        ifs.close();
    }
    return img;
}

// /* We try to avoid using third party libraries as much as possible however
//  * they are cases where this is not really possible. Displaying images to
//  * the screen requires to deal with the windows, key and mouse events.
//  * APIs to access these features are specific to the OS on which the program
//  * is compiled and are generally quite complex. GLFW  is cross-platform
//  * (Linux, Mac, Windows) and easy to compile. It can't really be used for
//  * creating complex GUIs however it perfectly suits the job for creating a
//  * simple image viewer.
//  *
//  * The most common technique to display an image using OpenGL is to convert
//  * the input image to a 2D texture and map this texture to a polygonal
//  * rectangle. The idea of the program is to read the PPM image first
//  * (PPM is an 8 bpp format and the image is internally converted to an 32
//  * bpp of float image), then create a window using the image's width and
//  * height. Once the window is created (with GLFW), a 2D texture is created
//  * using the following code:
//  */
// void initRgbTexture(const uint32_t &w, const uint32_t &h, const float *rawPixels)
// {
//     glBindTexture(GL_TEXTURE_2D, texID);
//     glPixelStorei (GL_UNPACK_ALIGNMENT, 1);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
//     glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
//     glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, w, h, 0, GL_RGB, GL_FLOAT, rawPixels);
// }
// 
// /* The draw function of the program, is responsible for setting up the
//  * viewport and the projection matrix (we can use an orthographic projection
//  * in that case) as well as rendering a rectangular polygon (which size is
//  * the same as the image). The previously created 2D texture is mapped to
//  * this rectangle effectively displaying the PPM file to the screen when
//  * the polygon is rendered.
//  */
// void draw(const unsigned char *screenPixels)
// {
//     int w, h;
//     glfwGetWindowSize(&w, &h);
//     glMatrixMode( GL_PROJECTION );
//     glLoadIdentity();
//     glViewport (0, 0, w, h); // optional
//     glOrtho(0, w, 0, h, -1, 1);
//     glMatrixMode( GL_MODELVIEW );
//     glLoadIdentity();
//     
//     glClearColor (0.3, 0.3, 0.3, 1.0);
//     glClear(GL_COLOR_BUFFER_BIT);
//     
//     glActiveTexture(GL_TEXTURE0);
//     glEnable(GL_TEXTURE_2D);
//     glBindTexture(GL_TEXTURE_2D, texID);
//     
//     cgGLEnableProfile (cgProfile);
//     
//     glBegin(GL_POLYGON);
//     glTexCoord2f(0, 1); glVertex2f(0, 0);
//     glTexCoord2f(1, 1); glVertex2f(w, 0);
//     glTexCoord2f(1, 0); glVertex2f(w, h);
//     glTexCoord2f(0, 0); glVertex2f(0, h);
//     glEnd();
//     
//     glDisable(GL_TEXTURE_2D);
//     
//     checkGlErrors("draw");
// }

#endif /* IMAGE_H */

/*
 * Best wishes;
 * Tariq Abuhashim, for iCub Facility
 * September, 2016
 * Thanks to Tim Bailey and Lorenzo Natale
 */