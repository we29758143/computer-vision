#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <iostream> 
#include "image.h"

#define M_PI 3.14159265358979323846

// HW1 #2.1
// Image& im: image to L1-normalize
void l1_normalize(Image& im)
  {
  
  float sum = 0;

  for (int i = 0; i < im.c; i++) {
    for (int j = 0; j < im.h; j++) {
      for (int k = 0; k < im.w; k++) {
        sum += im.pixel(k, j, i);
      }
    }
  }

  for (int i = 0; i < im.c; i++) {
    for (int j = 0; j < im.h; j++) {
      for (int k = 0; k < im.w; k++) {
        im.set_pixel(k, k, i, im.pixel(k, j, i) / sum);
      }
    }
  }
  
  }

// HW1 #2.1
// int w: size of filter
// returns the filter Image of size WxW
Image make_box_filter(int w)
  {
  assert(w%2); // w needs to be even
  
  // TODO: Implement the filter
  Image filter = Image(w, w);

  for (int i = 0; i < filter.w; i++) {
    for (int j = 0; j < filter.w; j++) {
      filter.set_pixel(j, i, 0, 1.0 / (w * w));
    }
  }
  
  return filter;
  }

// HW1 #2.2
// const Image&im: input image
// const Image& filter: filter to convolve with
// bool preserve: whether to preserve number of channels
// returns the convolved image
Image convolve_image(const Image& im, const Image& filter, bool preserve)
  {
  assert(filter.c==1);

  Image ret;
  int w_offset = filter.w / 2;
  int h_offset = filter.h / 2;

  if (preserve) {
    ret = Image(im.w, im.h, im.c);
    for (int i = 0; i < im.c; i++) {
      for (int j = 0 - h_offset; j < im.h - h_offset; j++) {
        for (int k = 0 - w_offset; k < im.w - w_offset; k++) {
          float origin = 0;
          for (int m = 0; m < filter.h; m++) {
            for (int n = 0; n < filter.w; n++) {
              origin += im.clamped_pixel(k + n, j + m, i) * filter(n, m, 0);
            }
          }
          // if (origin > 1) {
          //   origin = 1.0;
          // } else if (origin < 0) {
          //   origin = 0.0;
          // }
          ret.set_pixel(k + w_offset, j + h_offset, i, origin);

        }
      }
    }
  } else {
    ret = Image(im.w, im.h, 1);
  
    for (int j = 0 - h_offset; j < im.h - h_offset; j++) {
      for (int k = 0 - w_offset; k < im.w - w_offset; k++) {
        float origin = 0;
        for (int i = 0; i < im.c; i++) {
          for (int m = 0; m < filter.h; m++) {
            for (int n = 0; n < filter.w; n++) {
              origin += im.clamped_pixel(k + n, j + m, i) * filter(n, m, 0);
            }
          }
        }
        
        // if (origin > 1) {
        //   origin = 1; 
        // } else if (origin < 0) {
        //   origin = 0;
        // }
        ret.set_pixel(k + w_offset, j + h_offset, 0, origin);
      }
    }
    
  }


  
  // This is the case when we need to use the function clamped_pixel(x,y,c).
  // Otherwise you'll have to manually check whether the filter goes out of bounds
  
  // TODO: Make sure you set the sizes of ret properly. Use ret=Image(w,c,c) to reset ret
  // TODO: Do the convolution operator

  
  
  // Make sure to return ret and not im. This is just a placeholder
  return ret;
  }

// HW1 #2.3
// returns basic 3x3 high-pass filter
Image make_highpass_filter()
  {
  // TODO: Implement the filter
  Image filter = Image(3, 3);

  int arr[] = {0, -1, 0, -1, 4, -1, 0, -1, 0};
  int index = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      filter.set_pixel(j, i, 0, arr[index++]);
    }
  }
  
  return filter;
  
  }

// HW1 #2.3
// returns basic 3x3 sharpen filter
Image make_sharpen_filter()
  {
  // TODO: Implement the filter
  Image filter = Image(3, 3);

  int arr[] = {0, -1, 0, -1, 5, -1, 0, -1, 0};
  int index = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      filter.set_pixel(j, i, 0, arr[index++]);
    }
  }
  
  return filter;
  
  }

// HW1 #2.3
// returns basic 3x3 emboss filter
Image make_emboss_filter()
  {
  // TODO: Implement the filter
  Image filter = Image(3, 3);

  int arr[] = {-2, -1, 0, -1, 1, 1, 0, 1, 2};
  int index = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      filter.set_pixel(j, i, 0, arr[index++]);
    }
  }
  
  return filter;
  
  }

// HW1 #2.4
// float sigma: sigma for the gaussian filter
// returns basic gaussian filter
Image make_gaussian_filter(float sigma)
  {
  // TODO: Implement the filter
  int ceil_to_num = ((int)(ceil(sigma * 6)));
  int size = ceil_to_num % 2 == 0 ? ceil_to_num + 1 : ceil_to_num;
  float sum;
 
  Image filter = Image(size, size, 1);

  for (int i = -filter.h / 2; i < filter.h / 2; i++) {
    for (int j = -filter.w / 2; j < filter.w / 2; j++) {
      float val = (1 / (2.0 * M_PI * sigma * sigma)) * expf((- i * i - j * j) / (2.0 * sigma * sigma));
      sum += val;
      filter.set_pixel(j + filter.w / 2, i + filter.h / 2, 0, val);
    }
  }

  for (int i = -filter.h / 2; i < filter.h / 2; i++) {
    for (int j = -filter.w / 2; j < filter.w / 2; j++) {
      float val = filter.pixel(j + filter.w / 2, i + filter.h / 2, 0);
      filter.set_pixel(j + filter.w / 2, i + filter.h / 2, 0, 1 * val / sum);
    }
  }
  
  return filter;

  }


// HW1 #3
// const Image& a: input image
// const Image& b: input image
// returns their sum
Image add_image(const Image& a, const Image& b)
  {
  assert(a.w==b.w && a.h==b.h && a.c==b.c); // assure images are the same size
  
  Image ret = Image(a.w, a.h, a.c);

  for (int i = 0; i < a.c; i++) {
    for (int j = 0; j < a.h; j++) {
      for (int k = 0; k < a.w; k++) {
        ret.set_pixel(k, j, i, a.pixel(k, j, i) + b.pixel(k, j, i));
      }
    }
  }

  
  return ret;
  }

// HW1 #3
// const Image& a: input image
// const Image& b: input image
// returns their difference res=a-b
Image sub_image(const Image& a, const Image& b)
  {
  assert(a.w==b.w && a.h==b.h && a.c==b.c); // assure images are the same size
  
  Image ret = Image(a.w, a.h, a.c);

  for (int i = 0; i < a.c; i++) {
    for (int j = 0; j < a.h; j++) {
      for (int k = 0; k < a.w; k++) {
        ret.set_pixel(k, j, i, a.pixel(k, j, i) - b.pixel(k, j, i));
      }
    }
  }
  
  return ret;

  }

// HW1 #4.1
// returns basic GX filter
Image make_gx_filter()
  {
  Image filter = Image(3, 3);

  int arr[] = {-1, 0, 1, -2, 0, 2, -1, 0, 1};
  int index = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      filter.set_pixel(j, i, 0, arr[index++]);
    }
  }
  
  return filter;
  }

// HW1 #4.1
// returns basic GY filter
Image make_gy_filter()
  {
  Image filter = Image(3, 3);

  int arr[] = {-1, -2, -1, 0, 0, 0, 1, 2, 1};
  int index = 0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      filter.set_pixel(j, i, 0, arr[index++]);
    }
  }
  
  return filter;
  }

// HW1 #4.2
// Image& im: input image
void feature_normalize(Image& im)
  {
  assert(im.w*im.h); // assure we have non-empty image
  
  float min = 1.0;
  float max = 0.0;

  for (int i = 0; i < im.h; i++) {
    for (int j = 0 ; j < im.w; j++) {
      if (im.pixel(j, i, 0) > max) {
        max = im.pixel(j, i, 0);
      }
    }
  }

  for (int i = 0; i < im.h; i++) {
    for (int j = 0 ; j < im.w; j++) {
      if (im.pixel(j, i, 0) < min) {
        min = im.pixel(j, i, 0);
      }
    }
  }

  float range = max - min;

  if (range) {
    for (int i = 0; i < im.h; i++) {
      for (int j = 0; j < im.w; j++) {
        im.set_pixel(j, i, 0, (im.pixel(j, i, 0) - min) / range);
      }
    }
  } else {
    for (int i = 0; i < im.h; i++) {
      for (int j = 0; j < im.w; j++) {
        im.set_pixel(j, i, 0, 0);
      }
    }
  }
  
  }


// Normalizes features across all channels
void feature_normalize_total(Image& im)
  {
  assert(im.w*im.h*im.c); // assure we have non-empty image
  
  int nc=im.c;
  im.c=1;im.w*=nc;
  
  feature_normalize(im);
  
  im.w/=nc;im.c=nc;
  
  }


// HW1 #4.3
// Image& im: input image
// return a pair of images of the same size
pair<Image,Image> sobel_image(const Image& im)
  {
  // TODO: Your code here
  Image mag = Image(im.w, im.h, 1);
  Image dir = Image(im.w, im.h, 1);

  Image gx_filter = make_gx_filter();
  Image gy_filter = make_gy_filter();

  Image Gx = convolve_image(im, gx_filter, false);
  Image Gy = convolve_image(im, gy_filter, false);

  for (int i = 0; i < im.h; i++) {
    for (int j = 0; j < im.w; j++) {
      float Gx_val = Gx.pixel(j, i);
      float Gy_val = Gy.pixel(j, i);
      mag.set_pixel(j, i, 0, sqrtf(Gx_val * Gx_val + Gy_val * Gy_val));
      dir.set_pixel(j, i, 0, atan2f(Gy_val , Gx_val));
    }
  }

  
  return {mag,dir};
  }


// HW1 #4.4
// const Image& im: input image
// returns the colorized Sobel image of the same size
Image colorize_sobel(const Image& im)
  {
  
  // TODO: Your code here
  Image gau_filter = make_gaussian_filter(4);
  Image blur = convolve_image(im, gau_filter, true);
  blur.clamp();

  pair<Image,Image> res = sobel_image(blur);
  
  Image mag = res.first;
  Image theta = res.second;

  feature_normalize(mag);

  for (int i = 0; i < theta.h; i++) {
    for (int j = 0; j < theta.w; j++) {
      theta.set_pixel(j, i, 0, theta.pixel(j, i) / 2 * M_PI + 0.5);
    }
  }
  feature_normalize(theta);
  


  Image ret = Image(im.w, im.h, 3);
  for (int i = 0; i < im.h; i++) {
    for (int j = 0; j < im.w; j++) {
      // cout << theta.pixel(j,i,0) << endl;
      ret.set_pixel(j, i, 0, theta(j, i, 0));
    }
  }

  for (int ch = 1; ch < 3; ch++) {
    for (int i = 0; i < im.h; i++) {
      for (int j = 0; j < im.w; j++) {
        ret.set_pixel(j, i, ch, mag.pixel(j, i, 0));
      }
    }
  }

  hsv_to_rgb(ret);

  
  
  return ret;
  }


float distance(int x, int y, int i, int j) {
    return float(sqrt(pow(x - i, 2) + pow(y - j, 2)));
}

double gaussian(float x, double sigma) {
    return exp(-(pow(x, 2))/(2 * pow(sigma, 2))) / (2 * M_PI * pow(sigma, 2));

}

void applyBilateralFilter(Image im, Image& filter_image, int x, int y, int ch, int size, double sigma1, double sigma2) {
    double iFiltered = 0;
    double wP = 0;
    int neighbor_x = 0;
    int neighbor_y = 0;
    int half = size / 2;
    
    
    for(int i = 0; i < size; i++) {
        for(int j = 0; j < size; j++) {

            neighbor_x = x - (half - i);
            neighbor_y = y - (half - j);
            // cout <<neighbor_y << neighbor_x <<endl;
            double gi = gaussian(im.clamped_pixel(neighbor_x, neighbor_y, ch) - im.clamped_pixel(x, y, ch), sigma2);
           
            double gs = gaussian(distance(x, y, neighbor_x, neighbor_y), sigma1);
            double w = gi * gs;
            // iFiltered = iFiltered + source.at<uchar>(neighbor_x, neighbor_y) * w;
            iFiltered = iFiltered + im.clamped_pixel(neighbor_x, neighbor_y, ch) * w;
           
            wP = wP + w;
        }
    }
    iFiltered = iFiltered / wP;
    //cout << iFiltered <<endl;
    // filteredImage.at<double>(x, y) = iFiltered;
    filter_image.set_pixel(x, y, ch, iFiltered);
  
}

Image bilateralFilterOwn(Image im, double sigma1, double sigma2) {
    // Mat filteredImage = Mat::zeros(source.rows,source.cols,CV_64F);

    int ceil_to_num = ((int)(ceil(sigma1 * 6)));
    int size = ceil_to_num % 2 == 0 ? ceil_to_num + 1 : ceil_to_num;
    float sum;

    Image filter_image = Image(im.w, im.h, im.c);

    int width = im.w;
    int height = im.h;
    for(int ch = 0; ch < im.c; ch++){
      for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            applyBilateralFilter(im, filter_image, j, i, ch, size, sigma1, sigma2);
          }
      }
    }
   
    return filter_image;
}


// HW1 #4.5
// const Image& im: input image
// float sigma1,sigma2: the two sigmas for bilateral filter
// returns the result of applying bilateral filtering to im
Image bilateral_filter(const Image& im, float sigma1, float sigma2)
  {


  Image bi_filter_image = bilateralFilterOwn(im, sigma1, sigma2);
  
  
  return bi_filter_image;
  }



// HELPER MEMBER FXNS

void Image::feature_normalize(void) { ::feature_normalize(*this); }
void Image::feature_normalize_total(void) { ::feature_normalize_total(*this); }
void Image::l1_normalize(void) { ::l1_normalize(*this); }

Image operator-(const Image& a, const Image& b) { return sub_image(a,b); }
Image operator+(const Image& a, const Image& b) { return add_image(a,b); }
