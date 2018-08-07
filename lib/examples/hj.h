#include "srslte/srslte.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define PI 3.141592653589793

void DumpHex(const void* data, size_t size) {
  char ascii[17];
  size_t i, j;
  ascii[16] = '\0';
  for (i = 0; i < size; ++i) {
    printf("%02X ", ((unsigned char*)data)[i]);
    if (((unsigned char*)data)[i] >= ' ' && ((unsigned char*)data)[i] <= '~') {
      ascii[i % 16] = ((unsigned char*)data)[i];
    } else {
      ascii[i % 16] = '.';
    }
    if ((i+1) % 8 == 0 || i+1 == size) {
      printf(" ");
      if ((i+1) % 16 == 0) {
        printf("|  %s \n", ascii);
      } else if (i+1 == size) {
        ascii[(i+1) % 16] = '\0';
        if ((i+1) % 16 <= 8) {
          printf(" ");
        }
        for (j = (i+1) % 16; j < 16; ++j) {
          printf("   ");
        }
        printf("|  %s \n", ascii);
      }
    }
  }
}

void read_file(cf_t *buff, char *file_name){
  FILE *fp = NULL;
  int cnt = 0;
  int i = 0;
  int total = 0;
  if (!buff) {
    printf("buff null\n");
  }
  fp = fopen(file_name, "rb");
  if (!fp) {
    perror("fopen");
    exit(-1);
  }

  while(!feof(fp)) {
  //for (int i=0; i<2;i++) {
    cnt = fread(buff+i, sizeof(cf_t), 1, fp);
    total += cnt;
    i++;
  }
  printf("%d\n", total);
  //DumpHex(buff,16);
  fclose(fp);
  /*
  int cnt2 =0;
  fp = fopen("temp_zz", "wb");
  cnt2 = fwrite(buff, sizeof(cf_t), total, fp);
  printf("cnt2: %d\n",cnt2);
  fclose(fp);
  */
}

void freq_offset_apply(cf_t *input_buffer, cf_t *output_buffer, int num_samples, float samp_rate, float freq) {
  float a1,b1;
  float a2,b2;
  float sin_, cos_;
  int i = 0;

  for (i = 0; i < num_samples; i++) {
    a1 = creal(input_buffer[i]);
    b1 = cimag(input_buffer[i]);

    sin_ = sin(2*PI*freq*((float)i)/samp_rate);
    cos_ = cos(2*PI*freq*((float)i)/samp_rate);

    a2 = a1*cos_ + (-1*b1*sin_);
    b2 = a1*sin_ + b1*cos_;

    output_buffer[i] = a2 + b2*_Complex_I;
  }
}


