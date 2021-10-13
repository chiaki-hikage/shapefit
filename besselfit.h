#ifndef BESSELFIT_H
#define BESSELFIT_H

#include<iostream>
#include<math.h>
#include"psf.h"
using namespace std;

extern"C" {
  void fhti_(int*i,double*d,double*d,double*d,double*d,int*i,double*d,bool*b);
  void fht_(int*i,double*d,int*i,double*d);
}

template<class Type>
class Besselfit {

 private:
  const int ni,nj,npix,npara,ntmax,rmax,nb;
  int nfit,**ug;
  const Type len,gsize;
  Type cnu[2],nu[2],gam[2];
  Type ***ba, **m, *aobs, *afit, *para, *paraf, chi, chimin;
				
 public:
  Besselfit(const int& n, const Type& l, const int& nt, const int& b): 
    ni(n),nj(n),npix(n),len(l),gsize(l/n),ntmax(nt),
    nb(b),npara(b+3),rmax(7) {
    cnu[0]=1.67835, nu[0]=0.5, gam[0]=0.886227;
    cnu[1]=0.35, nu[1]=-0.85, gam[1]=6.220274;
    nfit=nb*(ntmax+1)*(ntmax+1);
    aobs=new Type[nfit];
    afit=new Type[nfit];
    m=initmat2(nfit,nfit);
    para=new Type[npara];
    paraf=new Type[npara];
    ug=initmat2_int(ni,nj);
    chimin=1.e30;
  }
  ~Besselfit(){
    delete [] aobs;
    delete [] afit;
    delete [] *m, delete [] m;
    delete [] para;
    delete [] paraf;
    delete [] *ug,delete [] ug;
  }

  void gridtouse(Type*cp,Type*r) const {

    // use grids within "rmax" times half-light radius from the centroid
    Type *g=new Type[2];
    for(int i=0;i<ni;i++) for(int j=0;j<nj;j++){
      g[0]=(i+0.5)*gsize;
      g[1]=(j+0.5)*gsize;
      if (dist(g,cp,2)<(*r)*rmax){
	ug[i][j]=1;
      } else {
	ug[i][j]=0;
      };
    }
    delete [] g;
  }

  void makebasis(Type*cp,Type*r,Type bpsf) {
    int i,j,ifit=0,l,nk=256;
    const double lkmin=log10(2*pi/len)-4.,lkmax=log10(2*pi/len*ni)+4,
      lkc=(lkmin+lkmax)/2, dlogk=(lkmax-lkmin)/nk, nc=(nk+1.)/2.;
    double kabs,kdash;
    Type *f=new Type[nk], *f2=new Type[nk], *lr=new Type[nk];
    Type gp[2],dis,fval,cos1,sin1,cos2,sin2,cos4,sin4;
    Psf<Type> psf(npix,bpsf);
    
    ba=initmat3(nfit,ni,nj);
    for (int p=0;p<=ntmax;p++) for (int nt=p;nt<=ntmax;nt++) 
      for (int ib=0;ib<nb;ib++) {
	for (int np=0;np<psf.nps();np++) {

	  // PSF-convolved basis functions
	  for (i=0;i<nk;i++) {
	    kabs=pow(10.,lkc+(i+1-nc)*dlogk);
	    kdash=kabs*(*r)/ cnu[ib];
	    f[i]=psf.kolmogorov_dgauss(np,kabs)
	      *kabs*pow(kdash,2*nt)
	      /pow(1.+kdash*kdash,nu[ib]+nt+1);
	  }

	  // Fourier transform of the basis functions
	  bessk(lkc,dlogk,nk,p,f,lr,f2);
	  //char*outfunc="basisfunc.dat"; ofwritefunc(outfunc,lr,f,nk);

	  for(i=0;i<ni;i++) for (j=0;j<nj;j++) {
	    gp[0]=(i+0.5)*gsize;
	    gp[1]=(j+0.5)*gsize;
	    dis=dist(cp,gp,2);
	    splint(lr,f,f2,nk,log10(dis),&fval);

	    //real and imaginary parts are separated into different bases
	    if (p==0) {
	      ba[ifit][i][j]+=fval;
	    } else {
	      cos1=(gp[0]-cp[0])/dis;
	      sin1=(gp[1]-cp[1])/dis;
	      cos2=cos1*cos1-sin1*sin1;
	      sin2=2*sin1*cos1;
	      if (p==1) {
		ba[ifit][i][j]+=fval*cos2;
		ba[ifit+1][i][j]+=fval*sin2;
	      } else if (p==2) {
		cos4=cos2*cos2-sin2*sin2;
		sin4=2*sin2*cos2;
		ba[ifit][i][j]+=fval*cos4;
		ba[ifit+1][i][j]+=fval*sin4;
	      } else {
		throw ios::failure("Expansion is upto 2nd order");
	      }
	    }
	  }
	}
	if (p==0) {ifit++;} else {ifit+=2;};
      }
    delete [] f;
    delete [] f2;
    delete [] lr;
  }

  void linearfit(Type**f,Type**n) const {

    int i,j,p,q;
    Type *b=new Type[nfit];
    Type **im=initmat2(nfit,nfit);

    // computing noise matrix
    for (p=0;p<nfit;p++) for (q=p;q<nfit;q++) {
      m[p][q]=0.;
      for (i=0;i<ni;i++) for (j=0;j<nj;j++) 
	if (ug[i][j]) m[p][q]+=ba[p][i][j]*ba[q][i][j]/(n[i][j]*n[i][j]);
    }

    for (p=0;p<nfit;p++) for (q=0;q<p;q++) m[p][q]=m[q][p];

    invmat(m,im,nfit);

    // obtain best fit values of coefficients of basis functions
    for (p=0;p<nfit;p++) {
      b[p]=0.;
      for (i=0;i<ni;i++) for (j=0;j<nj;j++) 
	if (ug[i][j]) b[p]+=ba[p][i][j]*f[i][j]/(n[i][j]*n[i][j]);	
    }

    for (p=0;p<nfit;p++) {
      aobs[p]=0.;
      for (q=0;q<nfit;q++) aobs[p]+=im[p][q]*b[q];
    }
    //for (p=0;p<nfit;p++) cout << aobs[p] << endl;

    delete [] b;
    delete [] *im, delete [] im;

    // output basis functions if needed
    char*outf="basis.dat";
    Type **g=initmat2(ni,nj);
    for (i=0;i<ni;i++) for (j=0;j<nj;j++) for (p=0;p<nfit;p++) 
      g[i][j]+=aobs[p]*ba[p][i][j];
    //ofwriteasc(outf,g);
    delete [] *g,delete [] g;
    delete [] **ba, delete [] *ba, delete [] ba;
  }

  void getellip(Type*gfit) const {
    for (int n=0;n<2;n++) gfit[n]=para[nb+1+n];
  }

  void outimg(char*outf,Type*cp,Type*r,Type bpsf,Type**f) {
    // compare a fitted image with the input one
    int i,j,p;
    makebasis(cp,r,bpsf);
    Type **g=initmat2(ni,nj);
    for (i=0;i<ni;i++) for (j=0;j<nj;j++)
      for (p=0;p<1;p++) g[i][j]+=afit[p]*ba[p][i][j];
    ofwriteasc2(outf,g,f);
    delete [] *g, delete [] g;
    delete [] **ba, delete [] *ba, delete [] ba;
  }

  void nonlinearfit() {
    int i,j;
    int *nfunk=new int;
    const float ftol=1.e-8;
    Type chi=1.e30,**parai,*chii;

    chii=new Type[npara+1];
    parai=initmat2(npara+1,npara);

    // give initial values of parameters
    initpara(parai);

    // find maximum likelihood
    for (i=0;i<=npara;i++) chii[i]=chisq(parai[i]);
    amoeba(parai,chii,npara,ftol,nfunk);
    for (i=0;i<=npara;i++){
      if (chii[i]<chi) {
	chi=chii[i]; 
	for (j=0;j<npara;j++) para[j]=parai[i][j];
      }
    }

    delete [] *parai, delete [] parai;
    delete [] chii;

    if (chi<chimin) chimin=chi;
    cout << "chisq, param (L[" << nb << "],del,gam1,gam2):" << endl;
    cout << chi << " ";
    for (i=0;i<npara;i++) cout << para[i] << " " ;
    cout << endl;
    //for (int i=0;i<nfit;i++) cout << afit[i] << endl;
  }

  void initpara(Type **parai) const {
    int i,j,ib;
    Type w1,del[nb];

    for (ib=0;ib<nb;ib++) {
      w1=1.+nu[ib];
      del[ib]=aobs[ib+nb] / aobs[ib] / w1;
    }
    for (ib=0;ib<nb;ib++) para[ib]=fabs(aobs[ib]);
    para[nb]=del[0];
    findepara("l");
    for (i=0;i<=npara;i++) for (j=0;j<npara;j++) parai[i][j]=para[j];
    for (i=1;i<=npara;i++) parai[i][i-1]=0.;
    
  }

  void findepara(char*c) const {
    int i,j,n,inmax=100,jnmax=100;
    Type ymin=1.e30,imin,jmin,di,dj;
    Type *g=new Type[2],y,*ptemp=new Type[npara];

    if (c=="l") {
      imin=-.8;
      di=1.6/(Type)inmax;
      jmin=-.8;
      dj=1.6/(Type)jnmax;
    } else {
      di=.16/(Type)inmax;
      imin=*(para+nb+1)-di*inmax/2;
      dj=.16/(Type)jnmax; 
      jmin=*(para+nb+2)-dj*jnmax/2;
    }

    for (n=0;n<=nb;n++) ptemp[n]=para[n];
    for (i=0;i<=inmax;i++) for (j=0;j<=jnmax;j++) {
      ptemp[nb+1]=imin+i*di;
      ptemp[nb+2]=jmin+j*dj;
      y=chisq(ptemp);
      if (y<ymin) {
	g[0]=ptemp[nb+1];
	g[1]=ptemp[nb+2];
	ymin=y;
      }
    }
    for (n=0;n<2;n++) para[nb+1+n]=g[n];
    delete [] ptemp;
    delete [] g;
  }

  Type chisq(Type *p) const {
    Type w1,w2,del,mdel,g1,g2,g3,g4;

    for (int ib=0;ib<nb;ib++) {
      afit[ib]=p[ib];
      w1=1.+nu[ib];
      del=p[nb];
      mdel=1.-del;
      g1=p[nb+1];
      g2=p[nb+2];
      afit[nb+ib]=p[ib]*w1*del;
      afit[nb*ntmax+ib*2+nb]=p[ib]*w1*g1*mdel;
      afit[nb*ntmax+ib*2+nb+1]=p[ib]*w1*g2*mdel;
      if (ntmax==2) {
	w2=(1.+nu[ib])*(2+nu[ib])/2.;
	g3=g1*g1-g2*g2;
	g4=2*g1*g2;
	afit[nb*ntmax+ib]=p[ib]*w2*(del*del+(g1*g1+g2*g2)*mdel*mdel/2);
	afit[nb*(ntmax+2)+ib*2+nb]=p[ib]*w2*g1*del*mdel*2;
	afit[nb*(ntmax+2)+ib*2+nb+1]=p[ib]*w2*g2*del*mdel*2;
	afit[nb*(ntmax+4)+ib*2+nb]=p[ib]*w2*g3*del*del/2;
	afit[nb*(ntmax+4)+ib*2+nb+1]=p[ib]*w2*g4*del*del/2;
      }
    }

    Type y=0;
    for (int i=0;i<nfit;i++) for (int j=0;j<nfit;j++) 
      y+=(aobs[i]-afit[i])*(aobs[j]-afit[j])*m[i][j];
    //!!! if nfit=npara, minimum chisq becomes 0, so add 1 
    return y+1;

  }

  Type**initmat2(const int& n1,const int& n2) const {
    int i,j;
    Type**f=new Type*[n1];
    f[0]=new Type[n1*n2];
    for (i=1;i<n1;i++) f[i]=f[0]+i*n2;
    for (i=0;i<n1;i++) for (j=0;j<n2;j++) f[i][j]=0;
    return f;
  }

  int**initmat2_int(const int& n1,const int& n2) const {
    int i,j;
    int**f=new int*[n1];
    f[0]=new int[n1*n2];
    for (i=1;i<n1;i++) f[i]=f[0]+i*n2;
    for (i=0;i<n1;i++) for (j=0;j<n2;j++) f[i][j]=0;
    return f;
  }

  Type***initmat3(const int& n1,const int& n2,const int& n3) const {
    int i,j,k;

    // pointers to first level
    Type***f=new Type**[n1];

    // pointers to second level
    f[0]=new Type*[n1*n2];
    for (i=1;i<n1;i++) f[i]=f[i-1]+n2;

    // pointers to third level
    f[0][0]=new Type[n1*n2*n3];
    for (i=1;i<n1;i++) f[i][0]=f[i-1][0]+n2*n3;
    for (i=0;i<n1;i++) for(j=1;j<n2;j++) f[i][j]=f[i][j-1]+n3;
    for (i=0;i<n1;i++) for(j=0;j<n2;j++) for(k=0;k<n3;k++) f[i][j][k]=(Type)0;
    return f;
  }

  Type****initmat4(const int& n1,const int& n2,const int& n3,const int& n4) const {
    int i,j,k,l;

    // pointers to first level
    Type****f=new Type***[n1];

    // pointers to second level
    f[0]=new Type**[n1*n2];
    for (i=1;i<n1;i++) f[i]=f[i-1]+n2;

    // pointers to third level
    f[0][0]=new Type*[n1*n2*n3];
    for (i=1;i<n1;i++) f[i][0]=f[i-1][0]+n2*n3;
    for (i=0;i<n1;i++) for(j=1;j<n2;j++) f[i][j]=f[i][j-1]+n3;

    // pointers to fourth level
    f[0][0][0]=new Type[n1*n2*n3*n4];
    for (i=1;i<n1;i++) f[i][0][0]=f[i-1][0][0]+n2*n3*n4;
    for (i=0;i<n1;i++) for(j=1;j<n2;j++) f[i][j][0]=f[i][j-1][0]+n3*n4;
    for (i=0;i<n1;i++) for(j=0;j<n2;j++) for (k=1;k<n3;k++)
      f[i][j][k]=f[i][j][k-1]+n4;

    for (i=0;i<n1;i++) for(j=0;j<n2;j++) 
      for(k=0;k<n3;k++) for(l=0;l<n4;l++) f[i][j][k][l]=0;

    return f;
  }

  void bessk(const Type& lkc,const Type& dlogk,int nk,
	     const int& p,Type *f,Type *lr,Type *f2) {


    bool flag=false;
    double kr=1.,rk,lrc,dlnk=dlogk*log(10.),nc=(nk+1.)/2.;
    int kropt=3,dir=1;
    double wsave[2*nk+2*(nk/2)+18];

    double p2=2*p;
    double q=0.;
    fhti_(&nk,&p2,&q,&dlnk,&kr,&kropt,wsave,&flag);

    if (!flag) throw ios::failure("Failed");
    lrc=log10(kr)-lkc;
    rk=pow(10,lkc-lrc);
    fht_(&nk,f,&dir,wsave);
    for (int i=0;i<nk;i++) {
      *(lr+i)=lrc+(i+1-nc)*dlogk;
      *(f+i)/=pow(10.,*(lr+i));
    }
    spline(lr,f,nk,f2);
  }

  bool iterate(Type *r, int& iter) {
    int n;
    Type del,delf;
    static int ntar;
    const float ftol=1.e-1; // iteration finishes if |delta| < ftol

    del=para[nb];
    if (fabs(del)<ftol) return true;

    if (iter==0) {
      ntar=0;
      for (int n=0;n<npara;n++) paraf[n]=para[n];
    } else {
      delf=paraf[nb];
      
      // if |delta| > former |delta| or chi > 2*chimin,  reset parameters

      if ((fabs(del)>fabs(delf)&&del*delf<0)||chi>chimin*2) {
	ntar++;
	*r=*(r-1);
	for (n=0;n<npara;n++) para[n]=paraf[n];
	del=delf;
      } else {
	ntar=0;
	for (n=0;n<npara;n++) paraf[n]=para[n];
      }
    } 

    // iteration continues with a new radius

    *(r+1)=(*r)*(ntar+sqrt(1-del))/(ntar+1);
    cout << "new radius: " << *(r+1) << endl;
    return false;
  }

  Type dist(Type*x,Type*y,const int& n) const {
    Type sm=0.;
    for (int i=0;i<n;i++) sm+=(*(x+i)-*(y+i))*(*(x+i)-*(y+i));
    return sqrt(sm);
  }


  void spline(Type*x, Type*y, const int n, Type*y2) const  {
    int i,k;
    double p,qn,sig,un,*u=new double[n],yp1,ypn;
    
    yp1=ypn=1.e30;

    *y2 = -0.5;
    *u=(3.0/(*(x+1)-*x))*((*(y+1)-*y)/(*(x+1)-*x)-yp1);
    for (i=1;i<n-1;i++) {
      sig=(*(x+i)-*(x+i-1))/(*(x+i+1)-*(x+i-1));
      p=sig*(*(y2+i-1))+2.0;
      *(y2+i)=(sig-1.0)/p;
      *(u+i)=(*(y+i+1)-*(y+i))/(*(x+i+1)-*(x+i)) - (*(y+i)-*(y+i-1))/(*(x+i)-*(x+i-1));
      *(u+i)=(6.0*(*(u+i))/(*(x+i+1)-*(x+i-1))-sig*(*(u+i-1)))/p;
    }
    qn=0.5;
    un=(3.0/(*(x+n-1)-*(x+n-2)))*(ypn-(*(y+n-1)-*(y+n-2))/(*(x+n-1)-*(x+n-2)));
    *(y2+n-1)=(un-qn*(*(u+n-2)))/(qn*(*(y2+n-2))+1.0);
    for (k=n-2;k>=0;k--) *(y2+k)=*(y2+k)*(*(y2+k+1))+*(u+k);
    delete [] u;
  }

  void splint(Type*xa, Type*ya, Type*y2a, const int n, const Type& x, Type *y) const {
    int klo,khi,k;
    Type h,b,a;
    
    klo=0;
    khi=n-1;
    while (khi-klo > 1) {
      k=(khi+klo) >> 1;
      if (*(xa+k) > x) khi=k;
      else klo=k;
    }
    h=*(xa+khi)-*(xa+klo);
    if (h == 0.0) throw ios::failure("Bad xa input to routine splint");
    a=(*(xa+khi)-x)/h;
    b=(x-*(xa+klo))/h;
    *y=a*(*(ya+klo))+b*(*(ya+khi))+((a*a*a-a)*(*(y2a+klo))+(b*b*b-b)*(*(y2a+khi)))*(h*h)/6.0;
  }

  void invmat(Type **m1, Type **m2, const int& n) const {
    
    int i,j,indx[n];
    Type d,evec[n];
    Type **m3=initmat2(n,n);
    for (i=0;i<n;i++) for (j=0;j<n;j++) m3[i][j]=m1[i][j];

    ludcmp(m3,n,indx,&d);
    for (j=0;j<n;j++) {
      for (i=0;i<n;i++) if (i==j) *(evec+i)=1.; else *(evec+i)=0.;
      lubksb(m3,n,indx,evec);
      for (int i=0;i<n;i++) m2[i][j]=*(evec+i);
     }
    delete [] *m3, delete [] m3;
  }
    
  void ludcmp(Type **a, const int& n, int *indx, Type *d) const {
    int i,imax,j,k;
    float const TINY=1.e-20;
    Type big,dum,sum,temp;
    Type *vv;

    vv=new Type[n];
    *d=1.0;
    for (i=0;i<n;i++) {
      big=0.0;
      for (j=0;j<n;j++)	if ((temp=fabs(a[i][j])) > big) big=temp;
      if (big == 0.0) throw ios::failure("Singular matrix in routine ludcmp");
      *(vv+i)=1.0/big;
    }
    for (j=0;j<n;j++) {
      for (i=0;i<j;i++) {
	sum=a[i][j];
	for (k=0;k<i;k++) sum -= a[i][k]*a[k][j];
	a[i][j]=sum;
      }
      big=0.0;
      for (i=j;i<n;i++) {
	sum=a[i][j];
	for (k=0;k<j;k++) sum -= a[i][k]*a[k][j];
	a[i][j]=sum;
	if ( (dum=*(vv+i)*fabs(sum)) >= big) {
	  big=dum;
	  imax=i;
	}
      }
      if (j != imax) {
	for (k=0;k<n;k++) {
	  dum=a[imax][k];
	  a[imax][k]=a[j][k];
	  a[j][k]=dum;
	}
	*d = -(*d);
	*(vv+imax)=*(vv+j);
      }
      indx[j]=imax;
      if (a[j][j] == 0.0) a[j][j]=TINY;
      if (j != n) {
	dum=1.0/(a[j][j]);
	for (i=j+1;i<n;i++) a[i][j] *= dum;
      }
    }
    delete [] vv;
  }

  void lubksb(Type **a, const int& n, int *indx, Type *b) const {
    int i,ii=-1,ip,j;
    Type sum;

    for (i=0;i<n;i++) {
      ip=*(indx+i);
      sum=*(b+ip);
      *(b+ip)=*(b+i);
      if (ii!=-1) for (j=ii;j<=i-1;j++) sum -= a[i][j]*(*(b+j));
      else if (sum) ii=i;
      *(b+i)=sum;
    }
    for (i=n-1;i>=0;i--) {
      sum=*(b+i);
      for (j=i+1;j<n;j++) sum -= a[i][j]*(*(b+j));
      *(b+i)=sum/a[i][i];
    }
  }

  void ofwritefunc(char*string,Type*x, Type*y, int n) const {
    ofstream ofs;
    ofs.open(string,ios::out); 
    if (!ofs) throw ios::failure("Failed to open output file");
    for (int i=0;i<n;i++) ofs << x[i] << " " << y[i] << endl;
    ofs.close();
  }

  void ofwriteasc(char*string,Type**f) const {
    ofstream ofs;
    ofs.open(string,ios::out); 
    if (!ofs) throw ios::failure("Failed to open output file");
    for (int j=0;j<nj;j++) {
      for (int i=0;i<ni;i++) ofs << i << " " << j << " " << f[i][j] << endl;
      ofs << endl;
    }
    ofs.close();
  }

  void ofwriteasc2(char*string,Type**f,Type**g) const {
    ofstream ofs;
    ofs.open(string,ios::out); 
    if (!ofs) throw ios::failure("Failed to open output file");
    for (int j=0;j<nj;j++) {
      for (int i=0;i<ni;i++) ofs << i << " " << j << " " << f[i][j] 
				 << " " << g[i][j] << endl;
      ofs << endl;
    }
    ofs.close();
  }

  void amoeba(Type **p, Type *y, int ndim, Type ftol, int *nfunk) const {
    int i,ihi,ilo,inhi,j,mpts=ndim+1,nmax=5000;
    Type rtol,sum,swap,ysave,ytry,*psum;
    
    psum=new Type[ndim];
    *nfunk=0;
    for (j=0;j<ndim;j++) {						\
      for (sum=0.0,i=0;i<mpts;i++) sum += p[i][j];			\
      psum[j]=sum;
    }
    for (;;) {
      ilo=0;
      ihi = y[0]>y[1] ? (inhi=1,0) : (inhi=0,1);
      for (i=0;i<mpts;i++) {
	if (y[i] <= y[ilo]) ilo=i;
	if (y[i] > y[ihi]) {
	  inhi=ihi;
	  ihi=i;
	} else if (y[i] > y[inhi] && i != ihi) inhi=i;
      }
      rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo]));
      if (rtol < ftol) {
	SWAP(y[0],y[ilo]);
	for (i=0;i<ndim;i++) SWAP(p[0][i],p[ilo][i]);
	break;
      }
      if (*nfunk >= nmax) throw ios::failure("NMAX exceeded");
      *nfunk += 2;
      ytry=amotry(p,y,psum,ndim,ihi,-1.0);
      if (ytry <= y[ilo])
	ytry=amotry(p,y,psum,ndim,ihi,2.0);
      else if (ytry >= y[inhi]) {
	ysave=y[ihi];
	ytry=amotry(p,y,psum,ndim,ihi,0.5);
	if (ytry >= ysave) {
	  for (i=0;i<mpts;i++) {
	    if (i != ilo) {
	      for (j=0;j<ndim;j++) p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
	      y[i]=chisq(psum);
	    }
	  }
	  *nfunk += ndim;
	  for (j=0;j<ndim;j++) {					\
	    for (sum=0.0,i=0;i<mpts;i++) sum += p[i][j];		\
	    psum[j]=sum;
	  }
	}
      } else --(*nfunk);
    }
    delete [] psum;
  }

/* (C) Copr. 1986-92 Numerical Recipes Software z!0(0. */

  Type amotry(Type **p, Type *y, Type *psum, int ndim, int ihi, Type fac) const
  {
    int j;
    Type fac1,fac2,ytry,*ptry;
    
    ptry=new Type[ndim];
    fac1=(1.0-fac)/ndim;
    fac2=fac1-fac;
    for (j=0;j<ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
    ytry=chisq(ptry);
    if (ytry < y[ihi]) {
      y[ihi]=ytry;
      for (j=0;j<ndim;j++) {
	psum[j] += ptry[j]-p[ihi][j];
	p[ihi][j]=ptry[j];
      }
    }
    delete [] ptry;
    return ytry;
  }
  /* (C) Copr. 1986-92 Numerical Recipes Software z!0(0. */

  void SWAP(Type a, Type b) const {
    Type c;
    c=a;
    a=b;
    b=c;
  }

};

#endif // BESSELFIT_H
