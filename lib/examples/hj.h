#include "srslte/srslte.h"
#include <stdio.h>

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

