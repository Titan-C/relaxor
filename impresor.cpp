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

// Imprime datos de los arreglos matricales
void array_print_bin(const std::vector< double * >& V, std::string ARCHIVO, unsigned int cols){
  std::ofstream file;
  file.open (ARCHIVO.c_str());
  
  for (unsigned int i=0;i<V.size();i++)
    file.write((char * )&V[i][0],cols*sizeof(double));
  
  file.close();
}
void import_data(std::vector< std::vector< double > >& M, std::string ARCHIVO, unsigned int filas, unsigned int columnas) {
  std::fstream file;
  file.open (ARCHIVO.c_str(), std::fstream::in);
  M.resize(filas);
  for(unsigned int i=0;i<filas;i++){
    M[i].resize(columnas);
    for(unsigned int j=0;j<columnas;j++)
      file>>M[i][j];
  }
  file.close();
}

//Numpy Arrays
int create_metadata(char preamble[PREAMBLE_LEN], char header[MAX_HDR_LEN],
                    char * descr, int fortran_order,
		    const std::vector<unsigned int> shape)
{
    unsigned char byte;
    uint16_t      hdrlen;
    unsigned int n, l, topad;
    /*
    See numpy/lib/format.py for details of the .npy file format.
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
    return l;
}

void npy_save(char* fname, char * descr, int fortran_order,
              const std::vector<unsigned int> shape, size_t sz, void* data)
{
    char preamble[PREAMBLE_LEN], header[MAX_HDR_LEN];
    FILE *fp;
    int l, n1, n2, n3, mtd_len;

    l = create_metadata(preamble, header, descr, fortran_order, shape);
    mtd_len = PREAMBLE_LEN + l;
    if ( mtd_len % 16 != 0 ) {
        printf("formatting error: metadata length %d not divisible by 16\n",
	       mtd_len);
        abort();
    }
    fp = fopen(fname, "w");
    n1 = fwrite(preamble, sizeof(char), PREAMBLE_LEN, fp);
    n2 = fwrite(header,   sizeof(char), l, fp);
    l = 1;
    for (unsigned int m=0; m<shape.size(); m++ ) l *= shape[m];
    n3 = fwrite(data, sz, l, fp);
    fclose(fp);
}

void npy_save_double(char* fname, int fortran_order,
                     const std::vector<unsigned int> shape, double* data)
{
    char descr[5];
    descr[0] = ENDIAN_CHAR;
    descr[1] = 'f';
    sprintf(descr+2, "%d", (int) sizeof(double));
    npy_save(fname, descr, fortran_order, shape, sizeof(double), data);
}