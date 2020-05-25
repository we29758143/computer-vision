#include <cstdio>
#include <cstring>
#include <cassert>
#include <cmath>
#include <iostream>
using namespace std;

#include "image.h"

using namespace std;

// HW0 #3
// const Image& im: input image
// return the corresponding grayscale image
Image rgb_to_grayscale(const Image& im)
  {
  assert(im.c == 3); // only accept RGB images
  Image gray(im.w,im.h,1); // create a new grayscale image (note: 1 channel)
  
  // TODO: calculate the pixels of 'gray'
  float para[] = {0.299, 0.587, 0.114};
  for (int i = 0; i < im.c; i++) {
    for (int j = 0; j < im.h * im.w; j++) {
      gray.data[j] += im.data[i * im.h * im.w + j] * para[i];
    }
  }
  
  // NOT_IMPLEMENTED();
  
  return gray;
  }



// Example function that changes the color of a grayscale image
Image grayscale_to_rgb(const Image& im, float r, float g, float b)
  {
  assert(im.c == 1);
  Image rgb(im.w,im.h,3);
  
  for(int q2=0;q2<im.h;q2++)for(int q1=0;q1<im.w;q1++)
    {
    rgb(q1,q2,0)=r*im(q1,q2);
    rgb(q1,q2,1)=g*im(q1,q2);
    rgb(q1,q2,2)=b*im(q1,q2);
    }
  
  return rgb;
  }




// HW0 #4
// Image& im: input image to be modified in-place
// int c: which channel to shift
// float v: how much to shift
void shift_image(Image& im, int c, float v)
  {
   
  assert(c>=0 && c<im.c); // needs to be a valid channel
  
  // TODO: shift all the pixels at the specified channel
 
  for (int y = 0; y < im.h; y++) {
    for (int x = 0; x < im.w; x++) {
      float value = im(x, y, c);
      float new_value = value + v;
      im.set_pixel(x, y, c, new_value);
    }
  }
  
  // clamp_image(im);
  // NOT_IMPLEMENTED();
  
  }

// HW0 #8
// Image& im: input image to be modified in-place
// int c: which channel to scale
// float v: how much to scale
void scale_image(Image& im, int c, float v)
  {
  assert(c>=0 && c<im.c); // needs to be a valid channel
  
  // TODO: scale all the pixels at the specified channel
  // for (int i = 0; i < im.h * im.w; i++) {
  //   im.data[im.h * im.w * c + i] *= v;
  // }

  for (int ch = 0; ch < im.c; ch++) {
    for (int y = 0; y < im.h; y++) {
      for (int x = 0; x < im.w; x++) {
        int origin = im.pixel_address(x, y, ch);
        int scale_value = origin * v;
        im.set_pixel(x, y, ch, scale_value);
      }
    }
  }
  // clamp_image(im);
  // NOT_IMPLEMENTED();
  }


// HW0 #5
// Image& im: input image to be modified in-place
void clamp_image(Image& im)
  {
  // TODO: clamp all the pixels in all channel to be between 0 and 1
  

  for (int ch = 0; ch < im.c; ch++) {
    for (int y = 0; y < im.h; y++) {
      for (int x = 0; x < im.w; x++) {

        float value = im.clamped_pixel(x, y, ch);
        if (value > 1) {
          value = 1;
        } else if (value < 0) {
          value = 0;
        }
        im.set_pixel(x, y, ch, value);
      }
    }
  }
  
  // NOT_IMPLEMENTED();
  
  }

// These might be handy
float max(float a, float b, float c)
  {
  return max({a,b,c});
  }

float min(float a, float b, float c)
  {
  return min({a,b,c});
  }


// HW0 #6
// Image& im: input image to be modified in-place
void rgb_to_hsv(Image& im)
  {
  assert(im.c==3 && "only works for 3-channels images");
  
  // TODO: Convert all pixels from RGB format to HSV format
  float r, g, b, min_value, max_value, c, h, s, v;
  int height = im.h;
  int weight = im.w;
  int k = 0;

  // cout << "rgb ----" << height << endl;
  // cout << "rgb ----" << weight << endl;


  for (int x = 0; x < weight; x++) {
    for (int y = 0; y < height; y++) {
      k++;
      r = im.data[x + y * weight];
      g = im.data[x + y * weight + weight * height];
      b = im.data[x + y * weight + weight * height * 2];

      min_value = min(r, g, b);
      max_value = max(r, g, b);

      c = max_value - min_value;
      v = max_value;

      if (v != 0) {
        s = c / v;
      } else {
        s = 0;
      }

      if (c == 0) {
        h = 0;
      } else {
        if (v == r) {
          h = (g - b) / c;
        } else if (v == g) {
          h = (b - r) / c + 2;
        } else if (v == b) {
          h = (r - g) / c +4;
        }
      }

      if (h < 0) {
        h = (h / 6) + 1;
      } else {
        h = (h / 6);
      }

      im.data[x + y * weight] = h;
      im.data[x + y * weight + weight * height] = s;
      im.data[x + y * weight + weight * height * 2] = v;


    }
  } 
  
  
  // NOT_IMPLEMENTED();
  }

// HW0 #7
// Image& im: input image to be modified in-place
void hsv_to_rgb(Image& im)
  {
  assert(im.c==3 && "only works for 3-channels images");
  
  // TODO: Convert all pixels from HSV format to RGB format
  

  float r, g, b, c, h, s, v, m, x, fm, mod;
  int height = im.h;
  int weight = im.w;
  // int k = 0;
  // cout << "hsv ----" << height << endl;
  // cout << "hsv ----" << weight << endl;


  for (int i = 0; i < weight; i++) {
    for (int j = 0; j < height; j++) {
      
      h = im.data[i + j * weight];
      s = im.data[i + j * weight + weight * height];
      v = im.data[i + j * weight + weight * height * 2];
     
      h *= 6;
      c = v * s;
      mod = fmod(h, 2);
      fm = abs(mod - 1);
      
      x = c * (1 - fm);
      m = v - c;

      if (h >= 0 && h < 1) {
        r = c;
        g = x;
        b = 0;
      } else if (h >= 1 && h < 2) {
        r = x;
        g = c;
        b = 0;
      } else if (h >= 2 && h < 3) {
        r = 0;
        g = c;
        b = x;
      } else if (h >= 3 && h < 4) {
        r = 0;
        g = x;
        b = c;
      } else if (h >= 4 && h < 5) {
        r = x;
        g = 0;
        b = c;
      } else {
        r = c;
        g = 0;
        b = x;
      } 

      im.data[i + j * weight] = r + m;
      im.data[i + j * weight + weight * height] = g + m;
      im.data[i + j * weight + weight * height * 2] = b + m;

    }
  }

  // NOT_IMPLEMENTED();
  
  }

// HW0 #9
// Image& im: input image to be modified in-place
void rgb_to_lch(Image& im)
  {
  assert(im.c==3 && "only works for 3-channels images");
  
  // TODO: Convert all pixels from RGB format to LCH format
  int weight = im.w;
  int height = im.h;
  int ch = im.c;
  // cout << "Im" << endl;
  // cout << im(0, 0, 0) << endl;
  // cout << im(0, 0, 1) << endl;
  // cout << im(0, 0, 2) << endl;

  Image rgb = Image(weight, height, ch);
  Image xyz = Image(weight, height, ch);
  Image luv = Image(weight, height, ch);

  // https://colorcalculations.wordpress.com/xyz-to-rgb/?fbclid=IwAR119uecDS-lmHP1AXJYYMO3T53Fyoi43llsBCcZl1qC7LTYqyrzAYTtgrc
  // sRGB to RGB
  // for (int i = 0; i < weight * height * ch; i++) {
  //   cout << rgb.data[i] << endl;
  // }
  for (int i = 0; i < weight * height * ch; i++) {
    rgb.data[i] = im.data[i] > 0.04045 ? powf((im.data[i] + 0.055) / 1.055, 2.4) : im.data[i] / 12.92;
  }
  
  // cout << "RGB" << endl;
  // cout << rgb(0, 0, 0) << endl;
  // cout << rgb(0, 0, 1) << endl;
  // cout << rgb(0, 0, 2) << endl;

  // https://www.image-engineering.de/library/technotes/958-how-to-convert-between-srgb-and-ciexyz
  // sRGB to XYZ
  float xyz_scale[3][3] = {{0.412453, 0.357580, 0.180423},
                          {0.212671, 0.715160, 0.072169}, 
                          {0.019334, 0.119193, 0.950277}};
  
  for (int i = 0; i < ch; ++i){
    for (int j = 0; j < ch; ++j){
      for (int k = 0; k < height * weight; ++k){
        xyz.data[i * height * weight + k] += xyz_scale[i][j] * rgb.data[j * height * weight + k];
      }
    }
  }

  // cout << "XYZ" << endl;
  // cout << xyz(0, 0, 0) << endl;
  // cout << xyz(0, 0, 1) << endl;
  // cout << xyz(0, 0, 2) << endl;

  // https://en.wikipedia.org/wiki/CIELUV
  // http://framewave.sourceforge.net/Manual/fw_function_020_0060_00330.html?fbclid=IwAR0jUlkGGhnkg8Zd_L38EI46oQZKYlcHKObENCf33FWnkybTnfIJVGX_aCQ
  // XYZ -> Luv
  // float xn = 0.312713;
  // float yn = 0.329016;
  float y_n = 1.0;

  // float un = 4 * xn / (-2 * xn + 12 * yn + 3);
  // float vn = 9 * yn / (-2 * xn + 12 * yn + 3);

  // openCV value
  float un = 0.19793943;
  float vn = 0.46831096;

  for (int i = 0; i < height * weight; ++i){
    float yr = xyz.data[height * weight + i] / y_n;
    if (yr > powf(6.0 / 29, 3)){
      luv.data[i] = 116 * powf(yr, 1.0/3) - 16;
    } else{
      luv.data[i] = powf(29.0/3, 3) * yr;
    }
    // calculate denominator
    float den = xyz.data[i] + 15*xyz.data[im.h*im.w+i] + 3*xyz.data[2*im.h*im.w+i];
    // u = x / den, v = y / den
    float u = 4 * xyz.data[i] / den;
    float v = 9 * xyz.data[im.h * im.w + i] / den;
    luv.data[im.h*im.w+i] = 13 * luv.data[i] * (u-un); // u
    luv.data[2*im.h*im.w+i] = 13 * luv.data[i] * (v-vn); // v
  }

  // cout << "LUV" << endl;
  // cout << luv(0, 0, 0) << endl;
  // cout << luv(0, 0, 1) << endl;
  // cout << luv(0, 0, 2) << endl;
  // Luv -> LCH
  // L -> U -> V
  for (int i = 0; i < height * weight; ++i){
    float u = luv.data[height * weight + i]; // second channel
    float v = luv.data[2 * height * weight + i]; // thrid channel
    im.data[i] = luv.data[i]; // L
    im.data[height * weight + i] = sqrtf(u*u + v*v); // C
    // atan2 (y,x) * 180 / PI;
    im.data[height * weight * 2+ i] = atan2(v, u);
  }
  // cout << "LCH" << endl;
  // cout << im(0, 0, 0) << endl;
  // cout << im(0, 0, 1) << endl;
  // cout << im(0, 0, 2) << endl;


  // NOT_IMPLEMENTED();
  }

// HW0 #9
// Image& im: input image to be modified in-place
void lch_to_rgb(Image& im)
  {
  assert(im.c==3 && "only works for 3-channels images");
  
  // TODO: Convert all pixels from LCH format to RGB format
  
  int weight = im.w;
  int height = im.h;
  int ch = im.c;

  // cout << "Reverse Im" << endl;
  // cout << im(0, 0, 0) << endl;
  // cout << im(0, 0, 1) << endl;
  // cout << im(0, 0, 2) << endl;

  Image rgb = Image(weight, height, ch);
  Image xyz = Image(weight, height, ch);
  Image luv = Image(weight, height, ch);
  int i, j, k;

  // LCH -> Luv
  for (int i = 0; i < weight * height; i++){
    // float l = im.data[i];
    float c = im.data[weight * height + i];
    float h = im.data[weight * height * 2 + i];
    luv.data[i] = im.data[i]; // L
    luv.data[im.h*im.w+i] = cosf(h) * c; // u
    luv.data[2*im.h*im.w+i] = sinf(h) * c; // v
  }

  // cout << "Reverse LUV" << endl;
  // cout << luv(0, 0, 0) << endl;
  // cout << luv(0, 0, 1) << endl;
  // cout << luv(0, 0, 2) << endl;
  
  // Luv -> XYZ
  // float xn = 0.312713;
  // float yn = 0.329016;
  float y_n = 1.0;

  // 某網站
  // float un = 0.2009;
  // float vn = 0.4610;

  float un = 0.19793943;
  float vn = 0.46831096;
 
  for (int i = 0; i < weight * height; i++){
    float u, v;
    
    if (!luv.data[i]){
       xyz.data[i] = 0.0;
       xyz.data[im.h*im.w+i] = 0.0;
       xyz.data[2*im.h*im.w+i] = 0.0;
       continue;
     }
    u = luv.data[im.h*im.w+i] / (13*luv.data[i]) + un;
    v = luv.data[2*im.h*im.w+i] / (13*luv.data[i]) + vn;
    if (luv.data[i] > 8){
      xyz.data[im.h*im.w+i] = y_n * powf((luv.data[i]+16)/116, 3);
    } else{
      xyz.data[im.h*im.w+i] = y_n * luv.data[i] * powf(3.0/29, 3);
    }
    xyz.data[i] = (9.0/4) * xyz.data[im.h*im.w+i] * u / v;
    xyz.data[2*im.h*im.w+i] = xyz.data[im.h*im.w+i] * (12-3*u-20*v) / (4*v);
  }

  // cout << "xyz" << endl;
  // cout << xyz.data[306961] << endl;

  // cout << "Reverse xyz" << endl;
  // cout << xyz(0, 0, 0) << endl;
  // cout << xyz(0, 0, 1) << endl;
  // cout << xyz(0, 0, 2) << endl;
  // XYZ -> RGB

  float mat[3][3] = {{0.412453, 0.357580, 0.180423},
                          {0.212671, 0.715160, 0.072169}, 
                          {0.019334, 0.119193, 0.950277}};
  float inverse_scale[3][3];
  float determinant = 0;
  //finding determinant
    for(i = 0; i < 3; i++)
        determinant = determinant + (mat[0][i] * (mat[1][(i+1)%3] * mat[2][(i+2)%3] - mat[1][(i+2)%3] * mat[2][(i+1)%3]));
    
    
    for(i = 0; i < 3; i++){
        for(j = 0; j < 3; j++)
            inverse_scale[i][j] = ((mat[(j+1)%3][(i+1)%3] * mat[(j+2)%3][(i+2)%3]) - (mat[(j+1)%3][(i+2)%3] * mat[(j+2)%3][(i+1)%3]))/ determinant;
    }


  for (i = 0; i < im.c; ++i){
    for (j = 0; j < im.c; ++j){
      for (k = 0; k < im.h*im.w; ++k){
        rgb.data[i*im.h*im.w+k] += inverse_scale[i][j] * xyz.data[j*im.h*im.w+k];
      }
    }
  }

  // cout << "Reverse rgb" << endl;
  // cout << rgb(0, 0, 0) << endl;
  // cout << rgb(0, 0, 1) << endl;
  // cout << rgb(0, 0, 2) << endl;

  // RGB -> sRGB
  for (i = 0; i < im.c*im.h*im.w; ++i){
    if (rgb.data[i] > 0.0031308){
      im.data[i] = (1.055) * powf(rgb.data[i], 1.0/2.4) - 0.055;
    } else{
      im.data[i] = 12.92 * rgb.data[i];
    }
  }


  // cout << "Reverse Im" << endl;
  // cout << im(0, 0, 0) << endl;
  // cout << im(0, 0, 1) << endl;
  // cout << im(0, 0, 2) << endl;

  // NOT_IMPLEMENTED();


  
  }



// Implementation of member functions
void Image::clamp(void) { clamp_image(*this); }
void Image::shift(int c, float v) { shift_image(*this,c,v); }
void Image::scale(int c, float v) { scale_image(*this,c,v); }

void Image::HSVtoRGB(void) { hsv_to_rgb(*this); }
void Image::RGBtoHSV(void) { rgb_to_hsv(*this); }
void Image::LCHtoRGB(void) { lch_to_rgb(*this); }
void Image::RGBtoLCH(void) { rgb_to_lch(*this); }
