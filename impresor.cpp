#include "impresor.h"

void out(double value, std::string ARCHIVO) {
  std::ofstream file;
  file.open (ARCHIVO.c_str(), std::fstream::app);
  file<<value<<"\n";
  file.close();
}
// Imprime datos de los arreglos vectoriales
void array_print(const std::vector< double >& V, std::string ARCHIVO, unsigned int colsize, double scale){
  std::fstream file;
  file.open (ARCHIVO.c_str(), std::fstream::out);
  unsigned int rowsize = V.size()/colsize;

  for(unsigned int i = 0; i<rowsize; i++){
    for(unsigned int j = 0; j<colsize;j++)
      file<<V[i*colsize+j]/scale<<" ";
    file<<"\n";
  }
  file.close();
}

//Numpy Arrays
void create_metadata(std::string ARCHIVO, char * descr, int fortran_order,
		     const std::vector<unsigned int> shape)
{
  char preamble[PREAMBLE_LEN], header[MAX_HDR_LEN];
  unsigned char byte;
  uint16_t      hdrlen;
  unsigned int n, l, topad, n1, n2, mtd_len;
  /*
   *    See numpy/lib/format.py for details of the .npy file format.
   */
  strcpy(header, "{'descr': '");
  strcat(header, descr);
  strcat(header, "', 'fortran_order': ");
  if ( fortran_order )
    strcat(header, "True, ");
  else
    strcat(header, "False, ");
  strcat(header, "'shape': (");
  for (unsigned int m=0; m<shape.size(); m++ ) {
    l = strlen(header);
    if ( shape[m] < 0 ) {
      printf("shape[%d] = %d is negative!\n", m, shape[m]);
      abort();
    }
    if ( l + MAX_INT_LEN + 4 >= MAX_HDR_LEN ) {
      printf("header too long\n");
      abort();
    }
    sprintf(header+l, "%d,", shape[m]);
  }
  l = strlen(header);
  if ( shape.size() > 1 ) header[l-1] = '\0'; // remove comma
  strcat(header, "), }");

  l = strlen(header);
  topad = 16 - ( PREAMBLE_LEN + l + 1 ) % 16;
  if ( l + topad + 1 > MAX_HDR_LEN ) {
    printf("header too long\n");
    abort();
  }
  for ( unsigned int m=0; m<topad; m++ ) header[l+m] = ' ';
  l += topad;
  header[l] = '\n';
  header[++l] = '\0';

  strcpy(preamble, MAGIC);
  n = strlen(preamble);
  byte = MAJOR;
  preamble[n++] = byte;
  byte = MINOR;
  preamble[n++] = byte;
  hdrlen = htole16(l);
  memcpy((void*) preamble+n, (void*) &hdrlen, 2);

  mtd_len = PREAMBLE_LEN + l;
  if ( mtd_len % 16 != 0 ) {
    printf("formatting error: metadata length %d not divisible by 16\n",
	   mtd_len);
    abort();
  }

  FILE *fp;
  fp = fopen(ARCHIVO.c_str(), "w");
  n1 = fwrite(preamble, sizeof(char), PREAMBLE_LEN, fp);
  n2 = fwrite(header,   sizeof(char), l, fp);
  fclose(fp);
}