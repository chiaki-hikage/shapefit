#ifndef PSF_H
#define PSF_H
#define KOLB2FWHM 2.9207
#define STD2FWHM 2.354

#include<iostream>
#include<math.h>
#include<rfftw.h>

using namespace std;

template<class Type>
class Psf {

 private:
  int ni,nj,ntot;
  Type bpsf;

 public:
  Psf(int n, Type s): ni(n), nj(n), bpsf(s), ntot(n*n) {}
  ~Psf(){}

  int nps(){
    return 2;
  }

  void smooth_kolmogorov_dgauss(Type **f, Type len) {
    int i,j,np;
    Type k,filter;
    westward_ho(f);
    for (i=0;i<=ni/2;i++) for (j=0;j<=nj/2;j++) {
      k=sqrt((Type)(i*i+j*j))*2*pi/len;
      filter=0.; for (np=0;np<2;np++) filter+=kolmogorov_dgauss(np,k);
      f[i][2*j]*=filter, f[i][2*j+1]*=filter;
      if ((i>0) && (i<ni/2)) f[ni-i][2*j]*=filter, f[ni-i][2*j+1]*=filter;
    }
    eastward_ho(f);
  }

  void smooth_kolmogorov(Type **f, Type len) {
    int i,j,np;
    Type k,filter;
    westward_ho(f);
    for (i=0;i<=ni/2;i++) for (j=0;j<=nj/2;j++) {
      k=sqrt((Type)(i*i+j*j))*2*pi/len;
      filter=kolmogorov(k);
      f[i][2*j]*=filter, f[i][2*j+1]*=filter;
      if ((i>0) && (i<ni/2)) f[ni-i][2*j]*=filter, f[ni-i][2*j+1]*=filter;
    }
    eastward_ho(f);
  }

  void smooth_gauss(Type **f, Type len) {
    int i,j,np;
    Type k,filter;
    westward_ho(f);
    for (i=0;i<=ni/2;i++) for (j=0;j<=nj/2;j++) {
      k=sqrt((Type)(i*i+j*j))*2*pi/len;
      filter=gauss(1.41421356*k);
      f[i][2*j]*=filter, f[i][2*j+1]*=filter;
      if ((i>0) && (i<ni/2)) f[ni-i][2*j]*=filter, f[ni-i][2*j+1]*=filter;
    }
    eastward_ho(f);
  }
  
  void westward_ho(Type **f) {
    int i,j;
    static rfftwnd_plan plan=0;
    
    //add a Nyquist plane to the array
    add_nyquist(f);
 
    // plan the Fourier transform
    if(!plan) plan=rfftw2d_create_plan(ni,nj,FFTW_REAL_TO_COMPLEX,
				       FFTW_ESTIMATE|FFTW_IN_PLACE);

    // compute the fastest Fourier transform of the array
    rfftwnd_one_real_to_complex(plan,(fftw_real*)f[0],(fftw_complex*)NULL);


    // get the correct normalization
    for(i=0;i<ni;i++) for(j=0;j<nj+2;j++) f[i][j]/=ntot;

    return;
  }

  void eastward_ho(Type **f) {
    int i,j;
    static rfftwnd_plan plan=0;

    /* plan the Fourier transform */
    if(!plan) plan=rfftw2d_create_plan(ni,nj,FFTW_COMPLEX_TO_REAL,
				       FFTW_ESTIMATE|FFTW_IN_PLACE);

    /* compute the reverse fastest Fourier transform of the array */
    rfftwnd_one_complex_to_real(plan,(fftw_complex*)f[0],(fftw_real*)NULL);
    
    /* remove the Nyquist plane */
    remove_nyquist(f);

    return;
  }

  void add_nyquist(Type **f) {

    int i,j;
    
    if(!f[0]) throw ios::failure("not enough memory for Nyquist plane");

    /* rearrange the values */
    for(i=ni-1;i>=0;i--) for(j=nj-1;j>=0;j--)
      f[0][j+i*(nj+2)]=f[0][j+i*nj];
    for(i=ni-1;i>=0;i--) for(j=nj+1;j>=nj;j--)
      f[0][j+i*(nj+2)]=0;
    
    /* calculate new pointers */
    for(i=1;i<ni;i++) f[i]=f[i-1]+nj+2;

    return;
  }

  void remove_nyquist(Type **f) {
    int i,j;

    /* rearrange the array values */
    for(i=0;i<ni;i++) for(j=0;j<nj;j++) f[0][j+i*nj]=f[0][j+i*(nj+2)];

    if(!f[0]) throw ios::failure("not enough memory for Nyquist plane");

    /* calculate new pointers */
    for(i=1;i<ni;i++) f[i]=f[i-1]+nj;

    return;
  }

  Type kolmogorov_dgauss(const int& n,const Type& k){
    if (n==0) {
      return 0.9*gauss(1.41421356*k);
    } else if (n==1) {
      return 0.09*gauss(2.82842712*k);
    } else { return 1; }
  }

  Type kolmogorov(const Type& k) {
    return exp(-pow(k*bpsf/KOLB2FWHM,5./3.));
  }

  Type gauss(const Type& k) {
    return exp(-pow(k*bpsf/STD2FWHM,2)/2.);
  }

};

#undef KOLB2FWHM
#undef STD2FWHM
#endif // PSF_H
