//main.h
#ifndef MAIN_H
#define MAIN_H
#define pi 3.141592653589793238462643
#define degrad pi/180

#include<iostream>
using namespace std;

template<class Type>
class Parlist {

 public:
  
  Parlist(){}
  ~Parlist(){}

  Type len;          // side lengh
  int npix;          // grid number of each side
  int ntry;          // number of realizations
  Type I0;           // brightness at r=r0
  Type r0;           // half light radius
  Type iser;         // sersic index
  Type bpsf;         // FWHM of PSF
  Type sn;           // S/N of image
  Type eps;          // input value of ellipticity
  Type phi;          // input value of orientation angle [radian]
  int ntmax;         // order of Taylor expansion (ntmax<=2)
  char*outf;         // file of gamma's (fit & input)
  char*fname;        // file of input image (binary)
  char*fname2;       // (option) file of images (fit & input) (ascii)
  char*ranellip;     // if you want to give random ellipticity, say "y" 

/* set default values for global variables */				   
  void defaultval(){
    len=1.;               
    npix=64;
    ntry=1;
    I0=1.;
    r0=0.05; r0*=len;
    iser=1.;
    bpsf=0.05; bpsf*=len;
    sn=1.e3;
    eps=0.3;
    phi=30.; phi*=degrad;
    ntmax=1;
    outf="fitval.dat";
    fname="img_bin.dat";
    fname2="img.dat";
    ranellip="n";
  }
  
  /* random command line into the global parameter list */
  void argument(int argc,char*argv[]) {
    char c;
   
    while((char)EOF!=(c=getopt(argc,argv,"l:p:n:I:r:i:s:b:e:f:t:o:d:c:a:"))) 
      switch(c) {
      case 'l': len      = atof(optarg); break;
      case 'p': npix     = atof(optarg); break;
      case 'n': ntry     = atoi(optarg); break;
      case 'I': I0       = atof(optarg); break;
      case 'r': r0       = atof(optarg); r0*=len; break;
      case 'i': iser     = atof(optarg); break;
      case 's': sn       = atof(optarg); break;
      case 'b': bpsf     = atof(optarg); bpsf*=len; break;
      case 'e': eps      = atof(optarg); break;
      case 'f': phi      = atof(optarg); phi*=degrad; break;
      case 't': ntmax    = atoi(optarg); break;
      case 'o': outf     =      optarg ; break;
      case 'd': fname    =      optarg ; break;
      case 'c': fname2   =      optarg ; break;
      case 'a': ranellip =      optarg ; break;
      }
  }


  
};

#undef degrad
#undef  pi
#endif // MAIN_H
