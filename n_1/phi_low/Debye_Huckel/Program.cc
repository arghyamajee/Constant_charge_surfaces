#include<iostream>
#include<cstdlib> //for using atof
#include<cmath> //for using fabs and M_PI
#include<limits> //for using numeric_limits
#include<vector>
using namespace std;

int main(int args, char*arg[])

{
  if(args!=8)
    {
      cerr<<"problem1 <ep> <s> <I> <dz> <N> <tol> <lB>" <<endl;
      return 1;
    }
  int N;
  double ep,s,I,dz,tol,phi_old,lB,k;
  ep=atof(arg[1]);
  s=atof(arg[2]);
  I=atof(arg[3]);
  dz=atof(arg[4]);
  N=atoi(arg[5]);
  tol=atof(arg[6]);
  lB=atof(arg[7]);

  double e=ep/(4*M_PI*lB);
  k=sqrt(8*M_PI*lB*I*6.022e-4/ep);

  double a1=1.+(k*k*dz*dz)/3.;
  double a2=2.*a1;
  double b=-1.+(k*k*dz*dz)/6.;
  double c=s*dz/e;

  vector <double> phi(N+1,0.);
  for (int i=1; true; ++i)
    {
      double Delta=0.;
      phi_old=phi[0];
      phi[0]=(1/a1)*(c-b*phi[1]);
      Delta=max(Delta,fabs(phi[0]-phi_old));
      for (int n=1; n<N; n++)
	{
	phi_old=phi[n];
	phi[n]=(-b/a2)*(phi[n-1]+phi[n+1]);
        Delta=max(Delta,fabs(phi[n]-phi_old));
        }
      phi_old=phi[N];
      phi[N]=(1/a1)*(c-b*phi[N-1]);
      Delta=max(Delta,fabs(phi[N]-phi_old));
      if (Delta<=tol)
	break;
    }

  cout.precision(numeric_limits<double>::digits10);
  cout << "# ep  = " << ep << endl;
  cout << "# e  = " << e << endl;
  cout << "# s  = " << s << endl;
  cout << "# I   = " << I << endl; 
  cout << "# dz   = " << dz << endl;
  cout << "# N   = " << N << endl; 
  cout << "# tol = " << tol << endl;
  cout << "# lB = " << lB << endl;
  cout << "#" << endl;
  cout << "# n   phi(z)" << endl;
  for (int n=0; n<=N; ++n)
    {
      cout << n  << "  "<< phi[n] << endl;
    }

  return 0;
}
