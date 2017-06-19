//
// author: Wu Shensen
//

#include <unistd.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>

#define PNG_DEBUG 3
#include <png.h>

typedef struct {
  unsigned long height;
  unsigned long width;
  unsigned long dpi;
  png_byte bit_depth;
  png_byte color_type;
  png_byte interlace;
  png_bytepp* rgba_storage;    // 3-D: channels * height * width: RGBA 4-channels
  float*** lab_storage;    // 3-D: channels * height * width: LAB 3-channels
} png_image, *png_image_p;




void debug(const char *s, ...);

png_image_p readPngFile(char* file_name);
int writePngFile(char* file_name, png_image_p output_image_p);

void convertRawPNGData(png_bytepp src_storage, png_bytepp* dest_storage,
                       unsigned long height, unsigned long width);
void convertBackRawPNGData(png_bytepp* src_storage, png_bytepp dest_storage,
                           unsigned long height, unsigned long width);

static inline void rgb2lab(png_bytep rgb, float* lab);
static inline void lab2rgb(float* lab, png_bytep rgb);

void convertRGBA2LABA(png_bytepp* src_rgba, float*** dest_lab,
                      unsigned long height, unsigned long width);
void convertLABA2RGBA(float*** src_lab, png_bytepp* dest_rgba,
                      unsigned long height, unsigned long width);





int main(int argc, char **argv) {
//  if (argc != 3)
//    debug("Usage: program_name <file_in> <file_out>");

//  read_png_file(argv[1]);
//  process_file();
//  write_png_file(argv[2]);

  png_image_p input_image_ptr = readPngFile("/data2/intellid/data/test/transparent_test.png");

  convertRGBA2LABA(input_image_ptr->rgba_storage, input_image_ptr->lab_storage,
                   input_image_ptr->height, input_image_ptr->width);

  convertLABA2RGBA(input_image_ptr->lab_storage, input_image_ptr->rgba_storage,
                   input_image_ptr->height, input_image_ptr->width);


  writePngFile("/data2/intellid/data/test/output.png", input_image_ptr);

  return 0;
}

png_image_p readPngFile(char* file_name) {
  /* open file and test for it being a png */
  FILE* fp = fopen(file_name, "rb");
  if (!fp) {
    debug("[read_png_file] File %s could not be opened for reading", file_name);
    return NULL;
  }

  unsigned char header[8];    // 8 is the maximum size that can be checked
  fread(header, 1, 8, fp);

  if (png_sig_cmp(header, 0, 8)) {
    fclose(fp);
    debug("[read_png_file] File %s is not recognized as a PNG file", file_name);
    return NULL;
  }

  /* initialize stuff */
  png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr) {
    fclose(fp);
    debug("[read_png_file] png_create_read_struct failed");
    return NULL;
  }

  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    png_destroy_read_struct(&png_ptr, png_infopp_NULL, png_infopp_NULL);
    fclose(fp);
    debug("[read_png_file] png_create_info_struct failed");
    return NULL;
  }

  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
    fclose(fp);
    debug("[read_png_file] Set Error Handler during init_io");
    return NULL;
  }

  png_init_io(png_ptr, fp);
  png_set_sig_bytes(png_ptr, 8);
  png_read_info(png_ptr, info_ptr);

  png_image_p input_image_p = (png_image*) malloc(sizeof(png_image));
  input_image_p->rgba_storage = NULL;
  input_image_p->lab_storage = NULL;
  input_image_p->width = png_get_image_width(png_ptr, info_ptr);
  input_image_p->height = png_get_image_height(png_ptr, info_ptr);
  input_image_p->dpi = png_get_pixels_per_meter(png_ptr, info_ptr);
  input_image_p->bit_depth = png_get_bit_depth(png_ptr, info_ptr);
  input_image_p->color_type = png_get_color_type(png_ptr, info_ptr);
  input_image_p->interlace = png_get_interlace_type(png_ptr, info_ptr);

  /* read file */
  if (setjmp(png_jmpbuf(png_ptr))) {
    free(input_image_p);
    png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
    fclose(fp);
    debug("[read_png_file] Set Error Handler during read_image");
    return NULL;
  }

  // read image pixels data
  png_bytepp raw_image_data_pp = (png_bytepp) malloc(sizeof(png_bytep) * input_image_p->height);
  for (int i = 0; i < input_image_p->height; i++) {
    raw_image_data_pp[i] = (png_bytep) malloc(png_get_rowbytes(png_ptr,info_ptr));
  }
  png_read_image(png_ptr, raw_image_data_pp);

  //initialize rgba storage
  input_image_p->rgba_storage = (png_bytepp*) malloc(sizeof(png_bytepp) * 4);
  for (int i = 0; i < 4; i++) {
    input_image_p->rgba_storage[i] = (png_bytepp) malloc(sizeof(png_bytep) * input_image_p->height);
    for (int j = 0; j < input_image_p->height; j++) {
      input_image_p->rgba_storage[i][j] = (png_bytep) malloc(sizeof(png_byte) * input_image_p->width);
    }
  }
  //initialize lab storage
  input_image_p->lab_storage = (float***) malloc(sizeof(float**) * 3);
  for (int i = 0; i < 3; i++) {
    input_image_p->lab_storage[i] = (float**) malloc(sizeof(float*) * input_image_p->height);
    for (int j = 0; j < input_image_p->height; j++) {
      input_image_p->lab_storage[i][j] = (float*) malloc(sizeof(float) * input_image_p->width);
    }
  }
  //convert raw image from height * width(combined by 4-channels) view to
  //4-channels * height * width.
  convertRawPNGData(raw_image_data_pp, input_image_p->rgba_storage,
                    input_image_p->height, input_image_p->width);

  //free resources.
  for (int y = 0; y < input_image_p->height; y++) {
    free(raw_image_data_pp[y]);
  }
  free(raw_image_data_pp);
  png_destroy_read_struct(&png_ptr, &info_ptr, png_infopp_NULL);
  fclose(fp);

  return input_image_p;
}

int writePngFile(char* file_name, png_image_p output_image_p) {
  /* create file */
  FILE *fp = fopen(file_name, "wb");
  if (!fp) {
    debug("[write_png_file] File %s could not be opened for writing", file_name);
    return -1;
  }

  /* initialize stuff */
  png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

  if (!png_ptr) {
    fclose(fp);
    debug("[write_png_file] png_create_write_struct failed");
    return -1;
  }

  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    png_destroy_write_struct(&png_ptr, png_infopp_NULL);
    fclose(fp);
    debug("[write_png_file] png_create_info_struct failed");
    return -1;
  }

  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    debug("[write_png_file] Error during init_io");
    return -1;
  }

  png_init_io(png_ptr, fp);

  /* write header */
  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    debug("[write_png_file] Error during writing header");
    return -1;
  }

  png_set_IHDR(png_ptr, info_ptr, output_image_p->width, output_image_p->height,
               output_image_p->bit_depth, output_image_p->color_type, output_image_p->interlace,
               PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

  png_set_pHYs(png_ptr, info_ptr, output_image_p->dpi, output_image_p->dpi, PNG_RESOLUTION_METER);

  png_write_info(png_ptr, info_ptr);

  /* write bytes */
  if (setjmp(png_jmpbuf(png_ptr))) {
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    debug("[write_png_file] Error during writing bytes");
    return -1;
  }

  png_bytepp image_bytes_p = (png_bytepp) malloc(sizeof(png_bytep) * output_image_p->height);
  for (int i = 0; i < output_image_p->height; i++) {
    image_bytes_p[i] = (png_bytep) malloc(4 * output_image_p->width);
  };

  convertBackRawPNGData(output_image_p->rgba_storage, image_bytes_p,
                        output_image_p->height, output_image_p->width);

  png_write_image(png_ptr, image_bytes_p);


  /* end write */
  if (setjmp(png_jmpbuf(png_ptr))) {
    for (int y = 0; y < output_image_p->height; y++) {
      free(image_bytes_p[y]);
    }
    free(image_bytes_p);
    png_destroy_write_struct(&png_ptr, &info_ptr);
    fclose(fp);
    debug("[write_png_file] Error during end of write");
    return -1;
  }

  png_write_end(png_ptr, NULL);

  /* cleanup heap allocation */
  for (int y = 0; y < output_image_p->height; y++) {
    free(image_bytes_p[y]);
  }
  free(image_bytes_p);
  png_destroy_write_struct(&png_ptr, &info_ptr);
  fclose(fp);
  return 0;
}

void convertRawPNGData(png_bytepp src_storage, png_bytepp* dest_storage,
                       unsigned long height, unsigned long width) {
  for (int i = 0; i < height; i++) {
    png_bytep row = src_storage[i];
    for (int j = 0; j < width; j++) {
      png_bytep ptr = &(row[j*4]);
      dest_storage[0][i][j] = ptr[0];
      dest_storage[1][i][j] = ptr[1];
      dest_storage[2][i][j] = ptr[2];
      dest_storage[3][i][j] = ptr[3];
    }
  }
}

void convertBackRawPNGData(png_bytepp* src_storage, png_bytepp dest_storage,
                           unsigned long height, unsigned long width) {
  for (int i = 0; i < height; i++) {
    png_bytep row = dest_storage[i];
    for (int j = 0; j < width; j++) {
      png_bytep ptr = &(row[j*4]);
      ptr[0] = src_storage[0][i][j];
      ptr[1] = src_storage[1][i][j];
      ptr[2] = src_storage[2][i][j];
      ptr[3] = src_storage[3][i][j];
    }
  }
}

static inline void rgb2lab(png_bytep rgb, float* lab) {
  float r, g, b, x, y, z, fx, fy, fz;
  r = rgb[0] / 255.0f;
  g = rgb[1] / 255.0f;
  b = rgb[2] / 255.0f;

  r = r > 0.04045 ? (float) pow((r + 0.055f) / 1.055f, 2.4f) : r / 12.92f;
  g = g > 0.04045 ? (float) pow((g + 0.055f) / 1.055f, 2.4f) : g / 12.92f;
  b = b > 0.04045 ? (float) pow((b + 0.055f) / 1.055f, 2.4f) : b / 12.92f;

  //Observer. = 2°, Illuminant = D65
  x = (r * 0.4124f + g * 0.3576f + b * 0.1805f) * 100.0f / 95.047f;
  y = r * 0.2126f + g * 0.7152f + b * 0.0722f;
  z = (r * 0.0193f + g * 0.1192f + b * 0.9505f) * 100.0f / 108.883f;

  //xyz to lab
  fx = x > 0.008856f ? (float) pow(x, 0.333333f) : (7.787f * x + 0.137931f);
  fy = y > 0.008856f ? (float) pow(y, 0.333333f) : (7.787f * y + 0.137931f);
  fz = z > 0.008856f ? (float) pow(z, 0.333333f) : (7.787f * z + 0.137931f);

  lab[0] = 116.0f * fy - 16.0f;
  lab[1] = 500.0f * (fx - fy);
  lab[2] = 200.0f * (fy - fz);
}

static inline void lab2rgb(float* lab, png_bytep rgb) {
  float r, g, b, x, y, z, fx, fy, fz;
  fy = (lab[0] + 16.0f) / 116.0f;
  fx = lab[1] / 500.0f + fy;
  fz = fy - lab[2] / 200.0f;

  x = pow(fx, 3.0f) > 0.008856f ? (float) pow(fx, 3.0f) : fx / 7.787f - 0.017713f;
  y = pow(fy, 3.0f) > 0.008856f ? (float) pow(fy, 3.0f) : fy / 7.787f - 0.017713f;
  z = pow(fz, 3.0f) > 0.008856f ? (float) pow(fz, 3.0f) : fz / 7.787f - 0.017713f;

  x = 0.95047f * x;
  z = 1.08883f * z;

  r = x *  3.2406f + y * -1.5372f + z * -0.4986f;
  g = x * -0.9689f + y *  1.8758f + z *  0.0415f;
  b = x *  0.0557f + y * -0.2040f + z *  1.0570f;

  r =  r > 0.0031308 ? (float) (1.055 * pow(r, 0.416667f) - 0.055f) : 12.92f * r;
  g =  g > 0.0031308 ? (float) (1.055 * pow(g, 0.416667f) - 0.055f) : 12.92f * g;
  b =  b > 0.0031308 ? (float) (1.055 * pow(b, 0.416667f) - 0.055f) : 12.92f * b;

  rgb[0] = (png_byte) (r * 255);
  rgb[1] = (png_byte) (g * 255);
  rgb[2] = (png_byte) (b * 255);
}

void convertRGBA2LABA(png_bytepp* src_rgba, float*** dest_lab,
                      unsigned long height, unsigned long width) {
  //temp variable.
  png_byte rgb_pixel[3];
  float lab_pixel[3];
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      rgb_pixel[0] = src_rgba[0][i][j];
      rgb_pixel[1] = src_rgba[1][i][j];
      rgb_pixel[2] = src_rgba[2][i][j];

      rgb2lab(rgb_pixel, lab_pixel);

      dest_lab[0][i][j] = lab_pixel[0];
      dest_lab[1][i][j] = lab_pixel[1];
      dest_lab[2][i][j] = lab_pixel[2];
    }
  }
}

void convertLABA2RGBA(float*** src_lab, png_bytepp* dest_rgba,
                      unsigned long height, unsigned long width) {
  //temp variable.
  png_byte rgb_pixel[3];
  float lab_pixel[3];
  for (int i = 0; i < height; i++) {
    for (int j = 0; j < width; j++) {
      //lab to xyz
      lab_pixel[0] = src_lab[0][i][j];
      lab_pixel[1] = src_lab[1][i][j];
      lab_pixel[2] = src_lab[2][i][j];

      lab2rgb(lab_pixel, rgb_pixel);

      dest_rgba[0][i][j] = rgb_pixel[0];
      dest_rgba[1][i][j] = rgb_pixel[1];
      dest_rgba[2][i][j] = rgb_pixel[2];
    }
  }
}

void debug(const char *s, ...) {
  va_list args;
  va_start(args, s);
  vfprintf(stderr, s, args);
  fprintf(stderr, "\n");
  va_end(args);
  return;
}
