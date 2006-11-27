#include<vector>
#include<assert.h>
using namespace std;

//#include <statmat.h>
//#include <map.h>
//#include <bandmat.h>

//----------------------------------------------------------
//--------------------KLASSENDEKLARATION--------------------
//----------------------------------------------------------

class envmatrix
  {

  private:
  vector<double> diag;                    //vector containing diagonal elements
  vector<double> env;                     //vector containing nonzero-entries
  vector<double> ldiag;                   //vector containing the diagonal elements
                                     //of the lower triangular factor of the
                                     //matrix if it is cholesky decomposed and
                                     //the elements of the matrix D if it is
                                     //rational cholesky decomposed.
  vector<double> lenv;                    //vector containing the envelope elements
                                     //of the lower triangular factor of the
                                     //matrix if it is cholesky decomposed and
                                     //the elements of the envelope of U without
                                     //if it is rational cholesky decomposed.
  vector<unsigned> xenv;             //vector of length dim+1 containing the
                                     //envelope-structure and the size of the
                                     //enevelope stored in xenv[dim]
  unsigned dim;                      //dimension of the matrix
  bool decomposed;                   //wether the matrix is cholesky
                                     //decomposed or not
  bool rational_decomposed;          //wether the matrix is rational cholesky
                                     //decomposed or not
  int bandwidth;                     //the bandwidth of the matrix. 0 indicates
                                     //a diagonal matrix, -1 indicates a matrix
                                     //with a more general envelope structure

  public:

  // DEFAULT CONSTRUCTOR
  envmatrix(void) { dim = 0; }

  // CONSTRUCTOR2
  // Initializes an envelope-matrix and stores the data provided in v in the
  // envelope of the matrix according to the storing scheme provided in xe
  // and the data provided in d in the diag-vector of the matrix
  // Therefore d.size()+1==xe.size() and v.size()==xe[d.size()]
  envmatrix(const vector<double> & d, const vector<double> & v,
            const vector<unsigned> & xe);

  // Copy CONSTRUCTOR
  envmatrix(const envmatrix & em);

  // OVERLOADED ASSIGNMENT OPERATOR
  const envmatrix & operator=(const envmatrix & em);

  //
  double operator()(const unsigned & i, const unsigned & j) const;

  // DESTRUCTOR
  ~envmatrix() {}

//-----------------------------------------------------------------------------
//--- Functions that decompose a matrix, solve systems of linear equations ----
//--------------- or compute the envelope of the inverse ----------------------
//-----------------------------------------------------------------------------

  // FUNCTION: decomp_rational
  // TASK: Computes the rational cholesky decomposition U'D^-1U of the calling
  //       matrix and stores U in lenv and D in ldiag
  void decomp_rational();

  // FUNCTION: inverse_envelope
  // TASK: computes the envelope of the inverse of the calling matrix and stores
  //       it in inv. inv is assumed to have the same envelope structure as the
  //       calling matrix.
  void inverse_envelope(envmatrix & inv);

//------------------------------------------------------------------------------
//---- Functions to assess elements or characteristics of the matrix------------
//------------------------------------------------------------------------------


  // FUNCTION: getL
  // TASK:     returns the (i,j) element of the cholesky-factor of the calling
  //           matrix
  double getL(const unsigned & i, const unsigned & j) const;

  // FUNCTION: getBandwidth
  // TASK:     returns the bandwidth of the matrix. -1 indicates a matrix
  //           without banded structure
  int getBandwidth(void) const;

  // FUNCTION: getDim
  // TASK: returns the dimension of the matrix
  unsigned getDim(void) const;

  // FUNCTION: getDiag
  // TASK: returns the ith element of diag, the vector with the diagonal
  //       elements of the calling matrix
  double getDiag(unsigned i);

  // FUNCTION: getDiag
  // TASK: returns diag, the vector with the diagonal elements
  //       of the calling matrix
  vector<double> getDiag();

  // FUNCTION: getEnv
  // TASK: returns the value of env(i)
  double getEnv(const unsigned i);

  // FUNCTION: getEnv
  // TASK: returns env, the envelope of the calling matrix
  vector<double> getEnv();

  //FUNCTION: getXenv
  //TASK: returns the value of xenv(i)
  unsigned getXenv(const unsigned i);

  // FUNCTION: getXenv
  // TASK: returns xenv, the vector containing the indices of the envelope
  //       of the calling matrix
  vector<unsigned> getXenv();

  // FUNCTION: traceOfProduct
  // TASK: returns the trace of the product A*B, where A is the calling matrix.
  double traceOfProduct(envmatrix & B);


//------------------------------------------------------------------------------
//---------- Functions to get pointers to elements of the matrix----------------
//------------------------------------------------------------------------------


  vector<double>::iterator getDiagIterator();

  vector<double>::iterator getEnvIterator();

  vector<unsigned>::iterator getXenvIterator();


//------------------------------------------------------------------------------
//------------- Functions for changing elements of the matrix-------------------
//------------------------------------------------------------------------------


  // FUNCTION setDecomposed
  // TASK: changes decomposed to the specified value. For decomposed=true the
  //       elements of ldiag and lenv are treated as the cholesky factor of the
  //       matrix
  void setDecomposed(const bool &t);

  // FUNCTION setRational_decomposed
  // TASK: changes rational_decomposed to the specified value. For
  //       rational_decomposed=true the elements of ldiag and lenv are treated
  //       as the rational cholesky factorisation of the matrix
  void setRational_decomposed(const bool &t);

}; //ende der klassendeklaration


//----------------------------------------------------------
//--------------------METHODENDEFINITION--------------------
//----------------------------------------------------------


envmatrix::envmatrix(const vector<double> & d, const vector<double> & v,
                        const vector<unsigned> & xe)
  {
  assert(d.size()+1==xe.size());
  assert(v.size()==xe[d.size()]);
  xenv=xe;
  diag=d;
  dim=diag.size();
  env=v;
  ldiag=vector<double>(dim,0);
  lenv=vector<double>(env.size(),0);
  decomposed=false;
  rational_decomposed=false;
  bandwidth=-1;
  }

double envmatrix::operator()(const unsigned & i, const unsigned & j) const
  {
  assert(i<dim);
  assert(j<dim);

  unsigned ih, jh;
  int kl, ku, zeroes;

  if(i>j)
    {
    ih = i;
    jh = j;
    }
  else if(i<j)
    {
    ih = j;
    jh = i;
    }
  else
    return diag[i];

  kl=xenv[ih];
  ku=xenv[ih+1];
  zeroes=ih-ku+kl;

  if(jh<zeroes)
    return double(0);
  else
    return env[kl+jh-zeroes];
  }

envmatrix::envmatrix(const envmatrix & em)
  {
  diag = em.diag;
  env = em.env;
  ldiag = em.ldiag;
  lenv = em.lenv;
  xenv = em.xenv;
  dim = em.dim;
  decomposed = em.decomposed;
  rational_decomposed = em.rational_decomposed;
  bandwidth = em.bandwidth;
  }

  // OVERLOADED ASSIGNMENT OPERATOR
const envmatrix & envmatrix::operator=(const envmatrix & em)
  {
  if (this == &em)
    return *this;

  diag = em.diag;
  env = em.env;
  ldiag = em.ldiag;
  lenv = em.lenv;
  xenv = em.xenv;
  dim = em.dim;
  decomposed = em.decomposed;
  rational_decomposed = em.rational_decomposed;
  bandwidth = em.bandwidth;

  return *this;
  }


//-----------------------------------------------------------------------------
//--- Functions that decompose a matrix or solve systems of linear equations---
//-----------------------------------------------------------------------------

void envmatrix::decomp_rational()
  {
  if(!rational_decomposed)
    {
    if(bandwidth==0)
      {
//      unsigned i;
      vector<double>::iterator ld = ldiag.begin();
      vector<double>::iterator d = diag.begin();
//      for(i=0; i<dim; i++)
      for(; d!=diag.end();++ld, ++d)
        {
//        ldiag[i]=1/diag[i];
        *ld = 1/ *d;
        }
      }
    else if(bandwidth==1)
      {
      unsigned i;
      vector<double>::iterator ld = ldiag.begin();
      vector<double>::iterator d = diag.begin();
      vector<double>::iterator le = lenv.begin();
      vector<double>::iterator e = env.begin();

//      ldiag[0]=1/diag[0];
      *ld=1/ *d;
//      lenv[0]=env[0]*ldiag[0];
      *le=*e * *ld;
      ++ld; ++d;
      for(i=1; i<dim-1; i++)
        {
//        ldiag[i]=1/(diag[i]-env[i-1]*lenv[i-1]);
        *ld=1/(*d- *e* *le);
        ++e; ++le;
//        lenv[i]=env[i]*ldiag[i];
        *le=*e* *ld;
        ++d; ++ld;
        }
//      ldiag[dim-1]=1/(diag[dim-1]-env[dim-2]*lenv[dim-2]);
      *ld=1/(*d-*e* *le);
      }
    else if(bandwidth==2)
      {
//      unsigned i, h;
      unsigned i;
      vector<double>::iterator ld = ldiag.begin();
      vector<double>::iterator d = diag.begin();
      vector<double>::iterator le = lenv.begin();
      vector<double>::iterator e = env.begin();

//      ldiag[0]=1/diag[0];
      *ld=1/ *d;
//      lenv[0]=env[0]*ldiag[0];
      *le=*e* *ld;
      ++ld; ++d;
//      ldiag[1]=1/(diag[1]-env[0]*lenv[0]);
      *ld=1/(*d-*e* *le);
      ++e; ++le;
//      for(i=2, h=1; i<dim; i++)
      for(i=2; i<dim; i++)
        {
//      lenv[h]=env[h]*ldiag[i-2];
        *le=*e* *(ld-1);
        ++le; ++e;
//        h++;
//      lenv[h]=(env[h]-lenv[h-1]*lenv[h-2]/ldiag[i-2])*ldiag[i-1];
        *le=(*e-*(le-1)**(le-2)/ *(ld-1))**ld;
        ++ld; ++d;
//        ldiag[i]=1/(diag[i]-lenv[h]*lenv[h]/ldiag[i-1]
//               -lenv[h-1]*lenv[h-1]/ldiag[i-2]);
        *ld=1/(*d-*le**le/ *(ld-1)
                 -*(le-1)**(le-1)/ *(ld-2));
        ++le; ++e;
//        h++;
        }
      }
    else if(bandwidth>2)
      {
//      unsigned i, j, k, h, l, m;
      unsigned i, k, h, l;
      h=0;
      vector<double>::iterator ld = ldiag.begin();
      vector<double>::iterator ldk;
      vector<double>::iterator ldm;
      vector<double>::iterator d = diag.begin();
      vector<double>::iterator le = lenv.begin();
      vector<double>::iterator lej;
      vector<double>::iterator lel;
      vector<double>::iterator e = env.begin();
      vector<unsigned>::iterator xe;

      for(i=0; i<(unsigned)bandwidth; i++)
        {
//        ldiag[i]=diag[i];
        *ld=*d;
        xe=xenv.begin();
        ldk=ldiag.begin();
        for(k=0; k<i; k++, h++)
          {
//          j=xenv[k];
//          lenv[h]=env[h];
          *le=*e;
//          m=0;
          ldm=ldiag.begin();
          lej=lenv.begin()+*xe;
          ++xe;
          lel=lenv.begin()+h-k;
//          for(l=h-k; l<h; l++, j++,m++)
          for(l=h-k; l<h; l++)
            {
//            lenv[h] -= lenv[j]*lenv[l]/ldiag[m];
            *le -= *lej**lel/ *ldm;
            ldm++; ++lej; ++lel;
            }
//          lenv[h] *= ldiag[k];
          *le *= *ldk;
//          ldiag[i] -= lenv[h]*lenv[h]/ldiag[k];
          *ld -= *le* *le/ *ldk;
          ++le; ++e; ++ldk;
          }
//        ldiag[i] = 1/ldiag[i];
        *ld = 1/ *ld;
        ++ld; ++d;
        }

      for(i=bandwidth; i<dim; i++)
        {
//        ldiag[i]=diag[i];
        *ld=*d;
        xe=xenv.begin()+i-bandwidth+1;
        ldk=ldiag.begin()+i-bandwidth;
        for(k=0; k<(unsigned)bandwidth; k++, h++)
          {
//          lenv[h]=env[h];
          *le=*e;
//          j=xenv[i-bandwidth+k+1]-k;
//          m=i-bandwidth;
          ldm=ldiag.begin()+i-bandwidth;
          lel=lenv.begin()+h-k;
          lej=lenv.begin()+*xe-k;
          ++xe;
//          for(l=h-k; l<h; l++, j++, m++)
          for(l=h-k; l<h; l++)
            {
//            lenv[h] -= lenv[j]*lenv[l]/ldiag[m];
            *le -= *lej* *lel/ *ldm;
            ++ldm; ++lel; ++lej;
            }
//          lenv[h]*=ldiag[i-bandwidth+k];
          *le*=*ldk;
//          ldiag[i] -= lenv[h]*lenv[h]/ldiag[m];
          *ld -= *le**le/ *ldm;
          ++le; ++e; ++ldk;
          }
//        ldiag[i] = 1/ldiag[i];
        *ld = 1/ *ld;
        ++ld; ++d;
        }
      }
    else
      {
      int i, k, j, h, jstart;
      h=0;

      vector<double>::iterator d = diag.begin();
      vector<double>::iterator ld = ldiag.begin();
      vector<double>::iterator ldk;
      vector<double>::iterator ldj;
      vector<double>::iterator e = env.begin();
      vector<double>::iterator le = lenv.begin();
      vector<unsigned>::iterator xe = xenv.begin();
      vector<unsigned>::iterator xek;

      for(i=0; i<(int)dim; i++)
        {
//        ldiag[i]=diag[i];
        *ld=*d;
//        for(k=i-xenv[i+1]+xenv[i]; k<i; k++, h++)
        xek=xenv.begin()+i-*(xe+1)+*xe;
        ldk=ldiag.begin()+i-*(xe+1)+*xe;
        for(k=i-*(xe+1)+*xe; k<i; k++, h++)
          {
//          lenv[h] = env[h];
          *le = *e;
//          if((xenv[i+1]-xenv[i])>(xenv[k+1]-xenv[k]))
          if((*(xe+1)-*xe)>(*(xek+1)-*xek))
            {
            jstart = i-*(xe+1)+*xe;
            }
          else
            {
            jstart = i-*(xek+1)+*xek;
            }
          ldj=ldiag.begin()+jstart;
          for(j=jstart; j<k; j++)
            {
//            lenv[h] -= getL(k,j)*getL(i,j)/getL(j,j);
            *le -= getL(k,j)*getL(i,j)/ *ldj;
            ++ldj;
            }
//          lenv[h] *= getL(k,k);
          *le *= *ldk;
//          ldiag[i] -= lenv[h]*lenv[h]/ldiag[k];
          *ld -= *le**le/ *ldk;
          ++le; ++e; ++xek; ++ldk;
          }
//        ldiag[i] = 1/ldiag[i];
        *ld = 1/ *ld;
        ++d; ++ld; ++xe;
        }
      }
    }
  rational_decomposed=true;
  decomposed=false;
  } //ende decomp_rational

void envmatrix::inverse_envelope(envmatrix & inv)
  {
  assert(xenv==inv.getXenv());
  assert(bandwidth==inv.getBandwidth());
  assert(dim==inv.getDim());
/*  if(bandwidth==0)
    {
    unsigned i;
    vector<double>::iterator invdiag=inv.getDiagIterator();
    vector<double>::iterator d=diag.begin();
    for(i=0; i<dim; i++, invdiag++, ++d)
      {
      *invdiag = 1/ *d;
      }
    }
   else if(bandwidth==1)
    {
    int i;
    vector<double>::iterator invdiag=inv.getDiagIterator()+dim-1;
    vector<double>::iterator ld=ldiag.end()-1;
    vector<double>::iterator invenv=inv.getEnvIterator()+env.size()-1;
    vector<double>::iterator le=lenv.end()-1;

    decomp_rational();

//    *invdiag = ldiag[dim-1];
    *invdiag = *ld;

    for(i=dim-2; i>=0; i--)
      {
//      *invenv = -lenv[i]* *invdiag;
      *invenv = -*le* *invdiag;
      invdiag--; ld--;
//      *invdiag = ldiag[i]-lenv[i]* *invenv;
      *invdiag = *ld-*le* *invenv;
      invenv--; le--;
      }
    }
  else if(bandwidth==2)
    {
    int i;
    vector<double>::iterator invdiag=inv.getDiagIterator()+dim-1;
    vector<double>::iterator ld=ldiag.end()-1;
    vector<double>::iterator invenv=inv.getEnvIterator()+env.size()-1;
    vector<double>::iterator lenvit=lenv.end()-1;

    decomp_rational();

//    *invdiag = ldiag[dim-1];
    *invdiag = *ld;
    *invenv = -*invdiag* *lenvit;
    invenv--; lenvit--;
//    *invenv = -inv(dim-1,dim-2)*getL(dim-2,dim-3)
//              -inv(dim-1,dim-1)*getL(dim-1,dim-3);
    *invenv = -*(invenv+1)**(lenvit-1)
              -*invdiag**lenvit;

    invdiag--; ld--;
//    *invdiag = ldiag[dim-2]-inv(dim-1,dim-2)*getL(dim-1,dim-2);
    *invdiag = *ld-*(invenv+1)**(lenvit+1);
    invenv--; lenvit--;

    for(i=dim-2; i>1; i--)
      {
//      *invenv = -inv(i,i)*getL(i,i-1)-inv(i,i+1)*getL(i+1,i-1);
      *invenv = -*invdiag**lenvit-*(invenv+2)**(lenvit+1);
      invenv--; lenvit--;
//      *invenv = -inv(i,i-1)*getL(i-1,i-2)-inv(i,i)*getL(i,i-2);
      *invenv = -*(invenv+1)**(lenvit-1)-*invdiag**lenvit;

      invdiag--; ld--;
//      *invdiag = ldiag[i-1]-inv(i,i-1)*getL(i,i-1)-inv(i+1,i-1)*getL(i+1,i-1);
      *invdiag = *ld-*(invenv+1)**(lenvit+1)-*(invenv+2)**(lenvit+2);
      invenv--; lenvit--;
      }

//    *invenv = -inv(1,1)*getL(1,0)-inv(1,2)*getL(2,0);
    *invenv = -*invdiag**lenvit-*(invenv+2)**(lenvit+1);
    invdiag--; ld--;
//    *invdiag = ldiag[0]-inv(1,0)*getL(1,0)-inv(2,0)*getL(2,0);
    *invdiag = *ld-*invenv**lenvit-*(invenv+1)**(lenvit+1);

    }
  else if(bandwidth>2)
    {
    int i, k, l, upperk, upperl;
    vector<double>::iterator invdiag=inv.getDiagIterator()+dim-1;
    vector<double>::iterator ld=ldiag.end()-1;
    vector<double>::iterator invenv=inv.getEnvIterator()+env.size()-1;

    decomp_rational();

    for(i=dim-1; i>0; i--)
      {
      if(dim-i>(unsigned)bandwidth)
        {
        upperk=bandwidth;
        }
      else
        {
        upperk=dim-i-1;
        }
//      *invdiag = ldiag[i];
      *invdiag = *ld;
      for(k=0; k<upperk; k++)
        {
        *invdiag -= inv(i+k+1,i)*getL(i+k+1,i);
        }
      invdiag--; ld--;

      if(i>bandwidth)
        {
        upperk=bandwidth;
        }
      else
        {
        upperk=i;
        }
      for(k=0; k<upperk; k++)
        {
//        if(i==(int)dim-1)
        if(dim-i+k<bandwidth)
          {
//          upperl = k+1;
          upperl = dim-i+k;
          }
        else
          {
          upperl = bandwidth;
          }
        *invenv=0;
        for(l=0; l<upperl; l++)
          {
          *invenv -= inv(i,i-k+l)*getL(i-k+l,i-k-1);
          }
        invenv--;
        }
      }
//    *invdiag = ldiag[0];
    *invdiag = *ld;
    for(k=0; k<(int)bandwidth; k++)
      {
      *invdiag -= inv(k+1,0)*getL(k+1,0);
      }
    }
  else */
// AB HIER WIRD GERECHNET !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    {
    int i, k, l, upperk, upperl;
    vector<double>::iterator invdiag=inv.getDiagIterator()+dim-1;
    vector<double>::iterator invenv=inv.getEnvIterator()+env.size()-1;

    vector<double>::iterator invenv1;
// fehler in der schleife bei diagonalmatrizen!!!!
    for(invenv1=inv.getEnvIterator(); invenv1<invenv; invenv1++)
      {
      *invenv1=0;
      }
// fehler in der schleife

    int maxbw=0;
    for(i=0; i<(int)dim; i++)
      {
      if((int)xenv[i+1]-(int)xenv[i]>maxbw)
        {
        maxbw=xenv[i+1]-xenv[i];
        }
      }

    decomp_rational();

    for(i=dim-1; i>=0; i--)// outerFor
      {
      *invdiag = ldiag[i];

      if(i+maxbw>=(int)dim)
        {
        upperk=dim;
        }
      else
        {
        upperk=i+maxbw+1;
        }

      for(k=i+1; k<upperk; k++)
        {
        *invdiag -= inv(k,i)*getL(k,i);
        }
      invdiag--;

      for(k=i-1; k>=(i-(int)xenv[i+1]+(int)xenv[i]); k--) // innerFor
        {

        if(k+maxbw>=(int)dim)
          {
          upperl=dim;
          }
        else
          {
          upperl=k+maxbw+1;
          }

        for(l=k+1; l<upperl; l++)
          {
          *invenv -= inv(i,l)*getL(l,k);
          }
        invenv--;
        } //end innerFor
      } //end outerFor
    } // ende bandwidth=-1
  inv.setDecomposed(false);
  inv.setRational_decomposed(false);
  } //ende inverse_envelope

//------------------------------------------------------------------------------
//------ Functions to assess elements or characteristics of the matrix----------
//------------------------------------------------------------------------------

double envmatrix::getL(const unsigned & i, const unsigned & j) const
  {
  int kl, ku, zeroes;
  unsigned ih, jh;
  assert(i<dim);
  assert(j<dim);
  if(i>j)
    {
    ih = i;
    jh = j;
    }
  else if(i<j)
    {
    ih = j;
    jh = i;
    }
  else
    return ldiag[i];
  kl=xenv[ih];
  ku=xenv[ih+1];
  zeroes=ih-ku+kl;
  if(jh<zeroes)
    return double(0);
  else
    return lenv[kl+jh-zeroes];
  }

int envmatrix::getBandwidth(void) const
  { return bandwidth; }

unsigned envmatrix::getDim(void) const
  { return dim; }


double envmatrix::getDiag(unsigned i)
  { return diag[i]; }

vector<double> envmatrix::getDiag()
  { return diag; }

double envmatrix::getEnv(const unsigned i)
  { return env[i]; }

vector<double> envmatrix::getEnv()
  { return env; }

unsigned envmatrix::getXenv(const unsigned i)
  { return xenv[i]; }

vector<unsigned> envmatrix::getXenv()
  { return xenv; }



//------------------------------------------------------------------------------
//---------- Functions to get pointers to elements of the matrix----------------
//------------------------------------------------------------------------------

vector<double>::iterator envmatrix::getDiagIterator()
  { return diag.begin(); }

vector<double>::iterator envmatrix::getEnvIterator()
  { return env.begin(); }

vector<unsigned>::iterator envmatrix::getXenvIterator()
  { return xenv.begin(); }


//------------------------------------------------------------------------------
//------------- Functions for changing elements of the matrix-------------------
//------------------------------------------------------------------------------

void envmatrix::setDecomposed(const bool &t)
  {
  decomposed=t;
  rational_decomposed=false;
  }

void envmatrix::setRational_decomposed(const bool &t)
  {
  decomposed=false;
  rational_decomposed=t;
  }

//--------------------------------
//----------------TRACEOFPRODUCT--
//--------------------------------
double envmatrix::traceOfProduct(envmatrix & B)
  {
  double trace=0;
  if((bandwidth==0)||(B.getBandwidth()==0))
    {
    vector<double>::iterator d1=diag.begin();
    vector<double>::iterator d2=B.getDiagIterator();
    unsigned i;

    for(i=0; i<dim; i++, d1++, d2++)
      {
      trace += *d1* *d2;
      }
    }
  else if(bandwidth>0 && B.getBandwidth()>0)
    {
    if(bandwidth==B.getBandwidth())
      {
      vector<double>::iterator d1=diag.begin();
      vector<double>::iterator d2=B.getDiagIterator();
      vector<double>::iterator env1=env.begin();
      vector<double>::iterator env2=B.getEnvIterator();
      for(; d1<diag.end(); d1++, d2++)
        {
        trace += *d1* *d2;
        }
      for(; env1<env.end(); env1++, env2++)
        {
        trace += 2* *env1* *env2;
        }
      }
    else if(bandwidth<B.getBandwidth())
      {
      vector<double>::iterator d1=diag.begin();
      vector<double>::iterator d2=B.getDiagIterator();
      vector<double>::iterator env1=env.begin();
      vector<double>::iterator env2=B.getEnvIterator();
      int i, k, diff,bbw;
      bbw=B.getBandwidth();
      diff=bbw-bandwidth;

      for(; d1<diag.end(); d1++, d2++)
        {
        trace += *d1* *d2;
        }
      for(i=0; i<bandwidth; i++)
        {
        for(k=0; (k<i) && (env1<env.end()) ; k++, env1++, env2++)
          {
          trace += 2* *env1* *env2;
          }
        }
      for(i=bandwidth; i<bbw; i++, env2+=(i-bandwidth))
        {
        for(k=0; (k<bandwidth) && (env1<env.end()) ; k++, env1++, env2++)
          {
          trace += 2* *env1* *env2;
          }
        }
      for(i=bandwidth; i<(int)dim; i++, env2+=diff)
        {
        for(k=0; (k<bandwidth) && (env1<env.end()) ; k++, env1++, env2++)
          {
          trace += 2* *env1* *env2;
          }
        }
      }
    else
      {
      vector<double>::iterator d1=diag.begin();
      vector<double>::iterator d2=B.getDiagIterator();
      vector<double>::iterator env1=env.begin();
      vector<double>::iterator env2=B.getEnvIterator();
      vector<double>::iterator env2end=B.getEnvIterator()+B.getXenv(dim);

      int i, k, diff, bbw;
      bbw=B.getBandwidth();
      diff=bandwidth-bbw;

      for(; d1<diag.end(); d1++, d2++)
        {
        trace += *d1* *d2;
        }
      for(i=0; i<bbw; i++)
        {
        for(k=0; (k<i) && (env2<env2end) ; k++, env1++, env2++)
          {
          trace += 2* *env1* *env2;
          }
        }
      for(i=bbw; i<bandwidth; i++, env1+=(i-bbw))
        {
        for(k=0; (k<bbw) && (env2<env2end) ; k++, env1++, env2++)
          {
          trace += 2* *env1* *env2;
          }
        }
      for(i=bandwidth; i<(int)dim; i++, env1+=diff)
        {
        for(k=0; (k<bbw) && (env2<env2end) ; k++, env1++, env2++)
          {
          trace += 2* *env1* *env2;
          }
        }
      }
    }
  else
    {
    vector<double>::iterator d1=diag.begin();
    vector<double>::iterator d2=B.getDiagIterator();
    vector<double>::iterator env1=env.begin();
    vector<double>::iterator env2=B.getEnvIterator();
    vector<unsigned>::iterator xenv1=xenv.begin();
    vector<unsigned>::iterator xenv2=B.getXenvIterator();

    unsigned i, k;

    for(; d1<diag.end(); d1++, d2++)
      {
      trace += *d1* *d2;
      }

    for(i=0; i<dim; i++, ++xenv1, ++xenv2)
      {
      if(*(xenv1+1)-*xenv1>=*(xenv2+1)-*xenv2)
        {
        for(k=0; k<*(xenv1+1)-*xenv1-*(xenv2+1)+*xenv2; k++)
          {
          ++env1;
          }
        for(k=0; k<*(xenv2+1)-*xenv2; k++)
          {
          trace += 2* *env1 * *env2;
          ++env1; ++env2;
          }
        }
      else
        {
        for(k=0; k<*(xenv2+1)-*xenv2-*(xenv1+1)+*xenv1; k++)
          {
          ++env2;
          }
        for(k=0; k<*(xenv1+1)-*xenv1; k++)
          {
          trace += 2* *env1 * *env2;
          ++env1; ++env2;
          }
        }
      }
    }
  return trace;
  }



extern "C" {
void myInversion(double* d, double* v, int* xe, int* n1, int* n2, int* n3){

	// vektoren für envFormat initialisieren
	vector<double> Vd(*n1,0);
	vector<double> Vv(*n2,0);
	vector<unsigned> Vxe(*n3,0);

	// cArrays auf vektoren übertragen
	for (size_t i=0; i<Vd.size(); ++i) Vd[i]=d[i];
	for (size_t i=0; i<Vxe.size(); ++i) Vxe[i]=xe[i];

	// target-objekt instanziieren
	// wichtig: der envelopeVektor mus mit null initialisiert sein
	envmatrix target(Vd, Vv, Vxe);

	// source-objekt instanziieren
	// der envelopeVektor des source-objekts muss beschrieben werden
	for (size_t i=0; i<Vv.size(); ++i) Vv[i]=v[i];
	envmatrix source(Vd, Vv, Vxe);

	// inversion aufrufen
	source.inverse_envelope(target);

	// vektoren aufcArrays übertragen
	for (size_t i=0; i<Vd.size(); ++i) d[i]=target.getDiag(i);
	for (size_t i=0; i<Vv.size(); ++i) v[i]=target.getEnv(i);
	for (size_t i=0; i<Vxe.size(); ++i) xe[i]=target.getXenv(i);

//	*n1=target.;
//	*n2=;

	return;

}
}


extern "C" {
void myTrace(double* Sd, double* Sv, int* Sxe, int* Sn1, int* Sn2, int* Sn3, double* BBd, double* BBv, int* BBxe, int* BBn1, int* BBn2, int* BBn3, double* result){

	// vektoren für envFormat initialisieren
	vector<double> SVd(*Sn1,0);
	vector<double> SVv(*Sn2,0);
	vector<unsigned> SVxe(*Sn3,0);

	vector<double> BBVd(*BBn1,0);
	vector<double> BBVv(*BBn2,0);
	vector<unsigned> BBVxe(*BBn3,0);

	// cArrays auf vektoren übertragen
	for (size_t i=0; i<SVd.size(); ++i) SVd[i]=Sd[i];
	for (size_t i=0; i<SVv.size(); ++i) SVv[i]=Sv[i];
	for (size_t i=0; i<SVxe.size(); ++i) SVxe[i]=Sxe[i];

	for (size_t i=0; i<BBVd.size(); ++i) BBVd[i]=BBd[i];
	for (size_t i=0; i<BBVv.size(); ++i) BBVv[i]=BBv[i];
	for (size_t i=0; i<BBVxe.size(); ++i) BBVxe[i]=BBxe[i];

	// matrix-objekte instanziieren
	envmatrix mS(SVd, SVv, SVxe);
	envmatrix mBB(BBVd, BBVv, BBVxe);

	// traceOfProduct aufrufen
	*result=mS.traceOfProduct(mBB);

	return;

}
}
