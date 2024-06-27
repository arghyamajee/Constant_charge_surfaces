#include<iostream>
#include<algorithm> //for using max function
#include<cstdlib> //for using atof
#include<cmath> //for using fabs
#include<fstream> //for using files
#include<limits> //for using numeric_limits
#include<vector> //for using vectors
using namespace std;

int M,N;
double A1, A2, B1, B2, C1, C2, D1, D2, E, F1, F2;

bool inside(int r, int s)
{
  if (abs(r)>M)
    return false;
  if ((s<0)||(s>N))
    return false;
  else
    return true;
}

bool corner(int r, int s)
{
  return((r==-M)&&(s==0))||
        ((r==-M)&&(s==N))||
        ((r== M)&&(s==0))||
        ((r== M)&&(s==N));
}

bool edge(int r, int s)
{
  return((abs(r)<M)&&(s==0))||
        ((abs(r)<M)&&(s==N))||
        ((abs(r)==M)&&(s>0)&&(s<N));
}

bool bulk(int r, int s)
{
  return(abs(r)<M)&&(s>0)&&(s<N);
}

bool wall(int r, int s)
{
  return(s==0)||(s==N);
}

bool le(int r, int s)
{
  return((r==-M)&&(s>0)&&(s<N));
}

double theta1(int r, int s)
{
  if(r>0) 
    return 0.;
  double weight=((r==0)? 0.5 : 1.);
  if(bulk(r,s)||le(r,s))
    return 4.*E*weight;
  else
   return 2.*E*weight;
}

double theta2(int r, int s)
{
  if((s>0)&&(s<N))
    return 0.;
  if(r>0)
    return 2.*F1;
  else if(r<0)
   return 2.*F2;
  else
    return F1+F2;
}

struct Entry
{
  double coeff;
  int r;
  int s;
};

int main (int args, char *arg[])
{
   if(args!=17)
   {
      cerr<<"Program <e1p> <e2p> <s1> <s2> <I1> <I2> <dx> <dz> <phi_b> <M> <N> <file_phi_low> <file_phi_up> <tol> <iter_max> <l_B>" <<endl;
      return 1;
   }

   double e1p = atof(arg[1]);
   double e2p = atof(arg[2]);
   double s1 = atof(arg[3]);
   double s2 = atof(arg[4]);
   double I1 = atof(arg[5]);
   double I2 = atof(arg[6]);
   double dx = atof(arg[7]);
   double dz = atof(arg[8]);
   double phi_b = atof(arg[9]);
   M = atoi(arg[10]);
   N = atoi(arg[11]);
   double tol = atof(arg[14]);
   int iter_max = atoi(arg[15]);
   double l_B = atof(arg[16]);
   
   double k1=sqrt(8*M_PI*l_B*I1*6.022e-4/e1p);
   double k2=sqrt(8*M_PI*l_B*I2*6.022e-4/e2p);

   double e1=e1p/(4*M_PI*l_B);
   double e2=e2p/(4*M_PI*l_B);

   A1=(e1*k1*k1*dx*dz)/9.+(e1*dz)/(3.*dx)+(e1*dx)/(3.*dz);
   B1=(e1*k1*k1*dx*dz)/18.-(e1*dz)/(3.*dx)+(e1*dx)/(6.*dz);
   C1=(e1*k1*k1*dx*dz)/18.+(e1*dz)/(6.*dx)-(e1*dx)/(3.*dz);
   D1=(e1*k1*k1*dx*dz)/36.-(e1*dz)/(6.*dx)-(e1*dx)/(6.*dz);

   A2=(e2*k2*k2*dx*dz)/9.+(e2*dz)/(3.*dx)+(e2*dx)/(3.*dz);
   B2=(e2*k2*k2*dx*dz)/18.-(e2*dz)/(3.*dx)+(e2*dx)/(6.*dz);
   C2=(e2*k2*k2*dx*dz)/18.+(e2*dz)/(6.*dx)-(e2*dx)/(3.*dz);
   D2=(e2*k2*k2*dx*dz)/36.-(e2*dz)/(6.*dx)-(e2*dx)/(6.*dz);

   E=(e2*k2*k2*phi_b*dx*dz)/4.;

   F1=(s1*dx)/2.;
   F2=(s2*dx)/2.;

   vector<double> phi_low(N+1, 0.);
   ifstream file_phi_low(arg[12]);   // Open file "arg[12]" for reading

   // Omit the header
   while (file_phi_low.peek() == '#')
     file_phi_low.ignore(numeric_limits<int>::max(), '\n');

   // Read in two columns
   for (int s = 0; s <= N; ++s)
     {
       double tmp1, tmp2;
       file_phi_low >> tmp1 >> tmp2;
       if (!file_phi_low.good())
       {
         cerr << "Error in file \"" << arg[12] << "\"" << endl;
         return 2;
       }
       phi_low[s] = tmp2 + phi_b;
     }
   
   file_phi_low.clear();
   file_phi_low.close();


   vector<double> phi_up(N+1, 0.);
   ifstream file_phi_up(arg[13]);   // Open file "arg[13]" for reading

   // Omit the header
   while (file_phi_up.peek() == '#')
     file_phi_up.ignore(numeric_limits<int>::max(), '\n');

   // Read in two columns
   for (int s = 0; s <= N; ++s)
     {
       double tmp1, tmp2;
       file_phi_up >> tmp1 >> tmp2;
       if (!file_phi_up.good())
	   {
	     cerr << "Error in file \"" << arg[13] << "\"" << endl;
         return 2;
       }
       phi_up[s] = tmp2;
     }

   file_phi_up.clear();
   file_phi_up.close();

   
   vector< vector<Entry> >V; vector<double>u; /* V is a vector or a matrix whose entries are again vectors having the type 'Entry', u is a vector whose entries are of double type*/
   for (int r=-M; r<=M; ++r)
    {
      double ep, si, ka;

      if (r>=0)
      {
         ep=e1;
         si=s1;
         ka=k1;
      }
      else
      {
        ep=e2;
        si=s2;
        ka=k2;
      }

      double A=(ep*ka*ka*dx*dz)/9.+(ep*dz)/(3.*dx)+(ep*dx)/(3.*dz);
      double B=(ep*ka*ka*dx*dz)/18.-(ep*dz)/(3.*dx)+(ep*dx)/(6.*dz);
      double C=(ep*ka*ka*dx*dz)/18.+(ep*dz)/(6.*dx)-(ep*dx)/(3.*dz);
      double D=(ep*ka*ka*dx*dz)/36.-(ep*dz)/(6.*dx)-(ep*dx)/(6.*dz);

      for (int s=0; s<=N; ++s)
       {

         vector<Entry>V1;
         double u1=0.;
         Entry tmp;

         //1st term

        if(r!=0)
        {
          tmp.coeff=(wall(r,s)? 2. : 4.)*A;
          tmp.r=r;
          tmp.s=s;
          V1.push_back(tmp);
        }
        else
	    {
	      tmp.coeff=(wall(r,s)? 1. : 2.)*(A1+A2);
          tmp.r=r;
          tmp.s=s;
          V1.push_back(tmp);
	    }

        //2nd term

        if(inside(r-1,s)&&(r!=0))
	    {
	      tmp.coeff=(wall(r,s)? 1. : 2.)*B;
	      tmp.r=r-1;
	      tmp.s=s;
          V1.push_back(tmp);
	    }
        else if(inside(r-1,s)&&(r==0))
	    {
	      tmp.coeff=(wall(r,s)? 1. : 2.)*B2;
	      tmp.r=r-1;
	      tmp.s=s;
          V1.push_back(tmp);
	    }
        else
	      u1-=(wall(r,s)? 1. : 2.)*B*phi_low[s];

        //3rd term

        if(inside(r,s-1)&&(r!=0))
	    {
	      tmp.coeff=2.*C;
	      tmp.r=r;
	      tmp.s=s-1;
	      V1.push_back(tmp);
	    }
        else if(inside(r,s-1)&&(r==0))
        {
	      tmp.coeff=C1+C2;
	      tmp.r=r;
	      tmp.s=s-1;
	      V1.push_back(tmp);
	    }

        //4th term

        if(inside(r+1,s))
	    {
	      tmp.coeff=(wall(r,s)? 1. : 2.)*B;
	      tmp.r=r+1;
	      tmp.s=s;
	      V1.push_back(tmp);
	    }
        else
	      u1-=(wall(r,s)? 1. : 2.)*B*phi_up[s];

        //5th term

        if(inside(r,s+1)&&(r!=0))
	    {
	      tmp.coeff=2.*C;
	      tmp.r=r;
	      tmp.s=s+1;
	      V1.push_back(tmp);
	    }
        else if(inside(r,s+1)&&(r==0))
        {
	      tmp.coeff=C1+C2;
	      tmp.r=r;
	      tmp.s=s+1;
	      V1.push_back(tmp);
	    }

        //6th term

        if(inside(r-1,s-1)&&(r!=0))
	    {
	      tmp.coeff=D;
	      tmp.r=r-1;
	      tmp.s=s-1;
	      V1.push_back(tmp);
	    }
        else if(inside(r-1,s-1)&&(r==0))
	    {
	      tmp.coeff=D2;
	      tmp.r=r-1;
	      tmp.s=s-1;
	      V1.push_back(tmp);
	    }
        else
	    {
	      if(s!=0)
	        u1-=D*phi_low[s-1];
	    }

        //7th term

        if(inside(r+1,s-1))
	    {
	      tmp.coeff=D;
	      tmp.r=r+1;
	      tmp.s=s-1;
	      V1.push_back(tmp);
	    }
        else
	    {
	      if(s!=0)
	        u1-=D*phi_up[s-1];
	    }

        //8th term

        if(inside(r-1,s+1)&&(r!=0))
	    {
	      tmp.coeff=D;
	      tmp.r=r-1;
	      tmp.s=s+1;
	      V1.push_back(tmp);
	    }
        else if(inside(r-1,s+1)&&(r==0))
	    {
	      tmp.coeff=D2;
	      tmp.r=r-1;
	      tmp.s=s+1;
	      V1.push_back(tmp);
	    }
        else
	    {
	      if(s!=N)
	        u1-=D*phi_low[s+1];
	    }

        //9th term

        if(inside(r+1,s+1))
	    {
	      tmp.coeff=D;
	      tmp.r=r+1;
	      tmp.s=s+1;
	      V1.push_back(tmp);
	    }
        else
	    {
	      if(s!=N)
	        u1-=D*phi_up[s+1];
	    }

        //10th and 11th term

        u1+=theta1(r,s)+theta2(r,s);

        V.push_back(V1);
        u.push_back(u1);
        }
    }

    vector< vector<double> >phi;
    vector<double>phi_int(N+1,0.);

    for (int r=-M; r<=M; ++r)
    {
      for (int s=0; s<=N; ++s)
        phi_int[s]=phi_low[s]*(M-r)/((double)(2*M))+phi_up[s]*(r+M)/((double)(2*M));
        phi.push_back(phi_int);
    }

    double delta;

    for(int iter=0; iter<iter_max; ++iter)
    {
      delta=0.;
      for(int r=-M; r<=M; ++r)
        for(int s=0; s<=N; ++s)
        {
          int i=(r+M)*(N+1)+s;
          double phi_new=u[i];
          for(int j=1; j<V[i].size(); ++j)
	        phi_new-=V[i][j].coeff*phi[V[i][j].r+M][V[i][j].s];
            phi_new/=V[i][0].coeff;
            delta=max(delta,fabs(phi_new-phi[r+M][s]));
            phi[r+M][s]=phi_new;
        }
        if (delta<=tol)
          break;
    }
    if (delta>tol)
      cerr <<"No convergence"<< endl;

    cout.precision(numeric_limits<double>::digits10);
    cout << "# e1p = " << e1p << endl;
    cout << "# e2p = " << e2p << endl;
    cout << "# e1 = " << e1 << endl;
    cout << "# e2 = " << e2 << endl;
    cout << "# s1 = " << s1 << endl;
    cout << "# s2 = " << s2 << endl;
    cout << "# I1 = " << I1 << endl;
    cout << "# I2 = " << I2 << endl;
    cout << "# k1 = " << k1 << endl;
    cout << "# k2 = " << k2 << endl;
    cout << "# dx = " << dx << endl;
    cout << "# dz = " << dz << endl;
    cout << "# phi_b = " << phi_b << endl;
    cout << "# M = " << M << endl;
    cout << "# N = " << N << endl;
    cout << "# tol = " << tol << endl;
    cout << "# iter_max = " << iter_max << endl;
    cout << "# delta = " << delta << endl;
    cout << "#" << endl;
    cout << "# r s phi(r,s)" << endl;

    for(int r=-M; r<=M; ++r)
      for(int s=0; s<=N; ++s)
      {
        cout << r << " " << s << " " << phi[r+M][s] << endl;
      }
      return 0;

}


 
