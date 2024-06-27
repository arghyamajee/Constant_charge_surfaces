#include<iostream>
#include<algorithm> //for using max function
#include<cstdlib> //for using atof function
#include<cmath> //for using math functions
#include<limits> //for using numeric limits
#include<vector> //for using vectors
using namespace std;

int N;
double e,si,I,dz,l_B,delta;

//define the Gradient: derivative of omega with respect to phi_s

double Gradient(int s, const vector<double> &phi)
{
  double sum=0.0;
  if(s==0)
  {
    double dphi = phi[s+1]-phi[s];
    double T1 = -6.022e-4*2*dz*I*l_B*l_B*((fabs(dphi)>numeric_limits<double>::epsilon())? (dphi*cosh(phi[s])-sinh(phi[s+1])+sinh(phi[s]))/(dphi*dphi) : 0.5*sinh(phi[s]));
    double T2 = e*l_B/(4*M_PI*dz)*(phi[s]-phi[s+1]);
    double T3 = -si*l_B*l_B;
    sum+=T1+T2+T3;
  }
  else if (s==N)
  {
    double dphi = phi[s]-phi[s-1];
    double T1 = 6.022e-4*2*dz*I*l_B*l_B*((fabs(dphi)>numeric_limits<double>::epsilon())? (dphi*cosh(phi[s])-sinh(phi[s])+sinh(phi[s-1]))/(dphi*dphi) : 0.5*sinh(phi[s]));
    double T2 = e*l_B/(4*M_PI*dz)*(phi[s]-phi[s-1]);
    double T3 = -si*l_B*l_B;
    sum+=T1+T2+T3;
  }
  else
  {
    double dphi = phi[s]-phi[s-1];
    double T1 = 6.022e-4*2*dz*I*l_B*l_B*((fabs(dphi)>numeric_limits<double>::epsilon())? (dphi*cosh(phi[s])-sinh(phi[s])+sinh(phi[s-1]))/(dphi*dphi) : 0.5*sinh(phi[s]));
    dphi = phi[s+1]-phi[s];
    double T2 = -6.022e-4*2*dz*I*l_B*l_B*((fabs(dphi)>numeric_limits<double>::epsilon())? (dphi*cosh(phi[s])-sinh(phi[s+1])+sinh(phi[s]))/(dphi*dphi) : 0.5*sinh(phi[s]));
    double T3 = e*l_B/(4*M_PI*dz)*(2*phi[s]-phi[s-1]-phi[s+1]);
    sum+=T1+T2+T3;
  }
  
  return sum;
}

int main (int args, char *arg[])
{
  if(args!=10)
  {
    cerr<<"Plu <e> <s> <I> <dz> <N> <l_B> <tol> <iter_max> <alpha>"<<endl;
    return 1;
  }
  e = atof(arg[1]);
  si = atof(arg[2]);
  I = atof(arg[3]);
  dz = atof(arg[4]);
  N = atoi(arg[5]);
  l_B = atof(arg[6]);
  double tol = atof(arg[7]);
  int iter_max = atoi(arg[8]);
  double alpha = atof(arg[9]);
  
  vector<double> phi(N+1,0.),phi_new(N+1,0.);
  
  for(int iter=0; iter<=iter_max; ++iter)
  {
    delta=0.0;
    for(int s=0; s<=N; ++s)
    {
      phi_new[s]=phi[s]-alpha*Gradient(s,phi);
      delta=max(delta,fabs(phi_new[s]-phi[s]));
    }
    if(delta<=tol)
      break;
    phi=phi_new;
  }
  
  cout.precision(numeric_limits<double>::digits10);
  cout << "#e=" << e <<endl;
  cout << "#s=" << si <<endl;
  cout << "#I=" << I <<endl;
  cout << "#dz=" << dz <<endl;
  cout << "#N=" << N <<endl;
  cout << "#l_B=" << l_B <<endl;
  cout << "#tol=" << tol <<endl;
  cout << "#iter_max=" << iter_max <<endl;
  cout << "#alpha=" << alpha <<endl;
  cout << "#delta=" << delta <<endl;
  cout << "#" <<endl;
  cout << "# s phi(s)" <<endl;
  
  for(int s=0; s<=N; ++s)
  {
    cout << s << " " << phi[s] << endl;
  }
  return 0;
}

