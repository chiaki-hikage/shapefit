//ellipticity.h
#ifndef ELLIPTICITY_H
#define ELLIPTICITY_H
#define pi 3.141592653589793238462643

#include<iostream>
#include"img.h"
#include"besselfit.h"
using namespace std;

template<class Type>
class Ellipticity {

 private:
  const int npix;
  const Type len,r0,iser,I0,bpsf;
  Type eps,phi;

 public: 
  Ellipticity(const int& n, const Type& l, const Type& i, const Type& r, 
	      const Type& is, const Type& s):
    npix(n),len(l), I0(i), r0(r),iser(is),bpsf(s) {}    

  ~Ellipticity(){}
  
  void giveellip(const Type& p, const Type& e){phi=p,eps=e;}

  void randomellip(long*idum) {
    Type epsmin=0.,epsmax=0.5;
    eps=epsmin+ran2(idum)*(epsmax-epsmin);
    phi=pi*(ran2(idum)-0.5);
  }

  void ellipimg(char*outf,const Type& sn,long*idum,Type& pnoi) const {
    Img<Type> img(npix,len);
    img.sersic(phi,eps,r0,I0,iser);
    img.addpsf(bpsf);
    pnoi=img.getsn(1.)/sn;
    img.addnoise(idum,pnoi);
    img.ofwritebin(outf);
  }

  void besselfit(char*inf,Type*gamfit,const int& ntmax,char*outf,const Type& pnoi) const {
    const int itermax=20; // maximum iteration number
    const int nb=2;       // number of basis functions for images (default 2)
    int i,iter=-1;
    Img<Type> gal(npix,len),noi(npix,len);
    Besselfit<Type> ellip(npix,len,ntmax,nb);
    Type cp[2],**f,**n,*r=new Type[itermax+1];

    gal.ifreadbin(inf);
    gal.findcenter(cp);
    gal.halfradius(cp,r);
    gal.halfradius_psfcorrect(r,bpsf);
    noi.constval(pnoi);

    f=gal.getaddress();
    n=noi.getaddress();

    while (iter<=itermax) {
      iter++;
      ellip.gridtouse(cp,r+iter);
      ellip.makebasis(cp,r+iter,bpsf);
      ellip.linearfit(f,n);
      ellip.nonlinearfit();
      if (ellip.iterate(r+iter,iter)) break;
    }
    if (*outf) ellip.outimg(outf,cp,r+iter,bpsf,f);
    ellip.getellip(gamfit);
    delete [] r;
  }

  void outfitval(Type*gamfit,const string& flag,char*outf,const double& period) const {
    int i;
    Type gamin[2],dif[2];
    static Type mdif[2],edif[2];
    static ofstream ofs;

    if (flag=="b") {
      ofs.open(outf,ios::out);
      ofs << "# gam1(fit), gam2(fit), gam1(input), gam2(input), comp time" << endl;
    }

    gamin[0]=eps*cos(2*phi); 
    gamin[1]=eps*sin(2*phi);
    ofs << gamfit[0] << " " << gamfit[1] 
      	<< " " << gamin[0] << " " << gamin[1] << " " << period << endl;

    cout << "input gam1: " << gamin[0] << ", gam2: " << gamin[1] << endl;
    cout << "fitted gam1: " << gamfit[0] << ", gam2: " << gamfit[1] << endl; 

    for (int n=0;n<2;n++) dif[n]=fabs(gamfit[n])-fabs(gamin[n]); 
    averg(dif,2,mdif,edif,flag);
    if (flag=="e") {
      cout << "Difference between abs of fit gamma and abs of input gamma" << endl;
      ofs << "#Difference between abs of fit gamma and abs of input gamma" << endl;
      for (i=0;i<2;i++) {
	cout << "gam" << i+1 << " mean: " << mdif[i] << ", std: " << edif[i] << endl;
	ofs << "#gam" << i+1 << " mean: " << mdif[i] << ", std: " << edif[i] << endl;
      }
      ofs.close();
    }

  }
  
  void averg(Type*a,const Type& n,Type*m,Type*e,const string& flag) const {
    int i;
    static int num;
    if (flag=="b") {
      num=1;
      for (i=0;i<n;i++) m[i]=a[i], e[i]=a[i]*a[i];
    } else if (flag=="e") {
      for (i=0;i<n;i++) {
	m[i]+=a[i], m[i]/=num;
	e[i]+=a[i]*a[i], e[i]/=num, e[i]-=m[i]*m[i];
	if (e[i]<=0.) e[i]=0.; else e[i]=sqrt(e[i]);
      }
    } else {
      num++;
      for (i=0;i<n;i++) m[i]+=a[i], e[i]+=a[i]*a[i];
    }
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

#undef pi
#endif // ELLIPTICITY_H
