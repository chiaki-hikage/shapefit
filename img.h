//img.h
#ifndef IMG_H
#define IMG_H

#include<iostream>
#include<fstream>
#include<math.h>
#include"psf.h"
#include"normsinv.h"
#include"ran2.h"

using namespace std;

template<class Type>
class Img {

 private:
  const int ni,nj,npix;
  const Type len,gsize;
  Type **f;

 public:

  Img(const int& n,const Type& l): ni(n), nj(n), npix(n), len(l), gsize(l/n) {
    f=new Type*[ni+2];
    f[0]=new Type[(ni+2)*nj];
    for(int i=1;i<ni+2;i++) f[i]=f[0]+i*nj;
  }

  ~Img(){
    delete [] f;
  }

  void sersic(const Type& phi,const Type& eps,const Type& r0,
	      const Type& I0,const Type& n) const {

    Type b1=1.6783, b4=7.6693;
    Type p[2]={len/2,len/2},g[2],bn;

    bn=b1+(b4-b1)*(n-1.)/3.;
    for(int i=0;i<ni;i++) for(int j=0;j<nj;j++) {
      g[0]=(i+0.5)*gsize;
      g[1]=(j+0.5)*gsize;
      for (int k=0;k<2;k++) g[k]-=p[k];
      rotate(g,phi);
      g[0]*=sqrt(1.-eps);
      g[1]*=sqrt(1.+eps);
      //!!!
      f[i][j]=pow(b1/r0,2)*exp(-bn)*
	I0*exp(-bn*(pow(norm(g,2)/r0,1./n)-1));
    }
  }

  void addpsf(const Type& bpsf) const {
    Psf<Type> psf(npix,bpsf);
    psf.smooth_kolmogorov_dgauss(f,len);
  }

  void addnoise(long *idum,const Type& sig) const {
    for (int i=0;i<ni;i++) for (int j=0;j<nj;j++) 
      f[i][j]+=sig*normsinv(ran2(idum));
  }

  double getsn(const Type& sig) const {
    double sn=0;
    for (int i=0;i<ni;i++) for (int j=0;j<nj;j++) sn+=pow(f[i][j]/sig,2);
    return sqrt(sn);
  }

  void constval(const Type& x) const {
    for (int j=0;j<nj;j++) for (int i=0;i<ni;i++) f[i][j]=x;
  }

  void ofwritebin(const char *string) const {
    ofstream ofs;
    ofs.open(string,ios::out|ios::binary); 
    if (!ofs) throw ios::failure("Failed to open output file");
    ofs.write((char*)f,sizeof(char)*ni*nj);
    ofs.close();
  }

  void ofwriteasc(const char *string) const {
    ofstream ofs;
    ofs.open(string,ios::out); 
    if (!ofs) throw ios::failure("Failed to open output file");
    for (int j=0;j<nj;j++) {
      for (int i=0;i<ni;i++) ofs << i << " " << j << " " << f[i][j] << endl;
      ofs << endl;
    }
    ofs.close();
  }

  void ifreadbin(const char *string) const {
    ifstream ifs;
    ifs.open(string,ios::in|ios::binary); 
    if (!ifs) throw ios::failure("Failed to open input file");
    ifs.read((char*)f,sizeof(char)*ni*nj);
    ifs.close();
  }

  void initcenter(Type *p) const {

    // initial center position is given by the center of brightest 3x3 grids 

    double max=0.,sm;
    int i,j,is,js;
    for (i=1;i<ni-1;i++) for (j=1;j<nj-1;j++) {
      sm=0.;
      for (is=-1;is<=1;is++) for (js=-1;js<=1;js++) sm+=f[i+is][j+js];
      if (sm>max) {
	p[0]=(i-0.5)*gsize;
	p[1]=(j-0.5)*gsize;
	max=sm;
      }
    }
    // cout << "init center x=" << p[0] << ",y=" << p[1] << endl;
  }

  void findcenter(Type *cpos) const {
    int i,j,k,ih,jh,iter=0;
    float dif=gsize;
    const int fac=10,itermax=10;
    const float hgsize=gsize/(float)fac;
    const float ftol=1.e-4;
    double sm;
    Type ipos[2],gpos[2],hgpos[2],r;

    initcenter(ipos);
    halfradius(ipos,&r);

    while(dif>gsize*ftol){
      iter++;
      sm=0.;
      // initialize new center position
      for (k=0;k<2;k++) *(cpos+k)=0.;
      for (j=0;j<nj;j++) for (i=0;i<ni;i++) {
	gpos[0]=(i+0.5)*gsize;
	gpos[1]=(j+0.5)*gsize;
	if (dist(gpos,ipos,2)<r+gsize) {	    
	  for (jh=0;jh<fac;jh++) for (ih=0;ih<fac;ih++) {
	    hgpos[0]=gpos[0]+(ih-fac/2+0.5)*hgsize;
	    hgpos[1]=gpos[1]+(jh-fac/2+0.5)*hgsize;
	    if (dist(hgpos,ipos,2)<r) {
	      for (k=0;k<2;k++) cpos[k]+=f[i][j]*hgpos[k];
	      sm+=f[i][j];
	    }
	  }
	}
      }
      for (k=0;k<2;k++) cpos[k]/=sm;
      dif=dist(ipos,cpos,2);
      // define a new center position
      for (k=0;k<2;k++) ipos[k]=cpos[k];
      if (iter==itermax) break;
    }
    cout << "center pos: x=" << ipos[0] << ",y=" << ipos[1] << endl;
  }

  void rotate(Type *x,const Type& phi) const {
    Type *temp=new Type[2];
    temp[0]=x[0]*cos(phi)+x[1]*sin(phi);
    temp[1]=-x[0]*sin(phi)+x[1]*cos(phi);
    for (int i=0;i<2;i++) x[i]=temp[i];
    delete [] temp;
  }

  void halfradius(Type*p,Type*r) const {
    int i,j,k,ir=0;
    Type r2;
    const int ikmax=20;
    const double dr=0.1,dk=2*pi/(double)ikmax;
    double sm,sh=0.;

    for (sm=0.,i=0;i<ni;i++) for (j=0;j<nj;j++) sm+=f[i][j];

    while(sh<sm/2){
      r2=dr*(ir+0.5);
      for(k=0;k<ikmax;k++){
	i=(int)(p[0]/gsize+r2*cos(k*dk))+1;
	j=(int)(p[1]/gsize+r2*sin(k*dk))+1;
	sh+=f[i][j]*r2*dr*dk;
      }
      ir++;
    }
    *r=r2*gsize;
    //cout << "half radii (w/ psf): " << *r << endl;
  }

  void halfradius_psfcorrect(Type*r,const Type& bpsf) const {
    if (*r>bpsf) *r=sqrt(*r**r-pow(bpsf,2));
    cout << "half radii (w/o psf): " << *r << endl;
  }

  Type**getaddress() const {
    return f;
  }

  Type dist(Type*x,Type*y,const int& n) const {
    Type sm=0.;
    for (int i=0;i<n;i++) sm+=(x[i]-y[i])*(x[i]-y[i]);
    return sqrt(sm);
  }
  
  Type norm(double*x,const int& n) const {
    int i;
    Type sm;
    for (sm=0.0,i=0;i<n;i++) sm+=x[i]*x[i];
    return sqrt(sm);
  }

  Type**initmat2(const int& n1,const int& n2) const {
    int i,j;
    Type**f=new Type*[n1];
    f[0]=new Type[n1*n2];
    for (i=1;i<n1;i++) f[i]=f[0]+i*n2;
    for (i=0;i<n1;i++) for (j=0;j<n2;j++) f[i][j]=0;
    return f;
  }

};

#endif // IMG_H
