#include <cmath>
#include "image.h"

using namespace std;

// HW1 #1
// float x,y: inexact coordinates
// int c: channel
// returns the nearest neibor to pixel (x,y,c)
float Image::pixel_nearest(float x, float y, int c) const
  {
  // Since you are inside class Image you can
  // use the member function pixel(a,b,c)
  
  // TODO: Your code here

  
  return clamped_pixel(roundf(x), roundf(y), c);

  
  // NOT_IMPLEMENTED();
  }

// HW1 #1
// float x,y: inexact coordinates
// int c: channel
// returns the bilinearly interpolated pixel (x,y,c)
float Image::pixel_bilinear(float x, float y, int c) const
  {
  // Since you are inside class Image you can
  // use the member function pixel(a,b,c)
  
  float v1, v2, v3, v4, a1, a2, a3, a4, q;

  v1 = clamped_pixel(floorf(x), floorf(y), c);
  v2 = clamped_pixel(ceilf(x), floorf(y), c);
  v3 = clamped_pixel(floorf(x), ceilf(y), c);
  v4 = clamped_pixel(ceilf(x), ceilf(y), c);

  a1 = (ceilf(x) - x) * (ceilf(y)- y);
  a2 = (x - floorf(x)) * (ceilf(y)- y);
  a3 = (ceilf(x) - x) * (y - floorf(y));
  a4 = (x - floorf(x)) * (y - floorf(y));
  
  q = v1 * a1 + v2 * a2 + v3 * a3 + v4 * a4;
  
  return q;
  }

// HW1 #1
// int w,h: size of new image
// const Image& im: input image
// return new Image of size (w,h,im.c)
Image nearest_resize(const Image& im, int w, int h)
  {
  Image ret(w,h,im.c);

  // w and h are int so we should do float transform
  float w_scale = (float) im.w / w;
  float h_scale = (float) im.h / h;

  // w = new image weight, h = new image height
  for (int i = 0; i < im.c; i++) {
    for (int j = 0; j < h; j++) {
      for (int k = 0; k < w; k++) {
        float w_map = w_scale * (k + 0.5) - 0.5;
        float h_map = h_scale * (j + 0.5) - 0.5;

        float value = im.pixel_nearest(w_map, h_map, i);
        ret.set_pixel(k, j, i, value);
      }
    }
  }
  
  
  return ret;
  }


// HW1 #1
// int w,h: size of new image
// const Image& im: input image
// return new Image of size (w,h,im.c)
Image bilinear_resize(const Image& im, int w, int h)
  {
  
  Image ret(w,h,im.c);

  // w and h are int so we should do float transform
  float w_scale = (float) im.w / w;
  float h_scale = (float) im.h / h;

  // w = new image weight, h = new image height
  for (int i = 0; i < im.c; i++) {
    for (int j = 0; j < h; j++) {
      for (int k = 0; k < w; k++) {
        float w_map = w_scale * (k + 0.5) - 0.5;
        float h_map = h_scale * (j + 0.5) - 0.5;

        float value = im.pixel_bilinear(w_map, h_map, i);
        ret.set_pixel(k, j, i, value);
      }
    }
  }
  
  
  return ret;
  }


