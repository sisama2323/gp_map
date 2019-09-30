#ifndef __SYMPERSYMMAT_
#define __SYMPERSYMMAT_

#include <cmath>
#include <vector>
#include <sstream>

namespace ifom
{

const double sqrt_pi = std::sqrt(M_PI);
const double sqrt_2 = std::sqrt(2);

// Storage class for a square symmetric and persymmetric matrix
class SymPersymMat
{
  size_t n_, n_half_;
  std::vector<std::vector<double>> contents_;
  
public:
  SymPersymMat(size_t n)
    : n_(n), n_half_( static_cast<size_t>(std::floor((n+1)/2.0)) ),
      contents_(n_half_,std::vector<double>())
  {
    for( size_t i = 0; i < n_half_; ++i )
      contents_[i].resize(n_-2*i);
  }
  
  size_t n() const
  { return n_; }

  size_t n_half() const
  {return n_half_;}
  
  SymPersymMat square() const
  {
    SymPersymMat M(n_);
    for( size_t i = 0; i < n_half_; ++i )
      for( size_t j = i; j < n_-i; ++j)
      {
        M(i,j) = operator()(i,0) * operator()(0,j);
        for( size_t k = 1; k < n_; ++k)
          M(i,j) += operator()(i,k) * operator()(k,j);
      }
    return M;
  }  
  
  double& operator()(size_t i, size_t j)
  {
    // top triangle
    if( 0 <= i && i < n_half_ && i <= j && j < n_-i )
      return contents_[i][j-i];
    // right triangle
    if( n_-j <= i && i <= j && n_half_ <= j && j < n_ )
      return contents_[n_-j-1][j-i];    
    // left triangle
    if( 0<= j && j < n_half_ && j + 1 <= i && i < n_-j)
      return contents_[j][i-j];
    // bottom triangle
    if( n_half_ <= i && i < n_ && n_-i <= j && j < i)
      return contents_[n_-i-1][i-j];
  }
  
  // https://stackoverflow.com/questions/5762042/const-overloaded-operator-function-and-its-invocation
  double operator()(size_t i, size_t j) const
  {
    // top triangle
    if( 0 <= i && i < n_half_ && i <= j && j < n_-i )
      return contents_[i][j-i];
    // right triangle
    if( n_-j <= i && i <= j && n_half_ <= j && j < n_ )
      return contents_[n_-j-1][j-i];    
    // left triangle
    if( 0<= j && j < n_half_ && j + 1 <= i && i < n_-j)
      return contents_[j][i-j];
    // bottom triangle
    if( n_half_ <= i && i < n_ && n_-i <= j && j < i)
      return contents_[n_-i-1][i-j];
  }
    
  std::string to_str() const
  {
    std::stringstream ss;
    for( size_t i = 0; i < n_; ++i )
    {
      for( size_t j = 0; j < n_; ++j )
        ss << operator()(i,j) << " ";
      ss << std::endl;
    }
    return ss.str();
  }
  
  friend std::ostream& operator<<(std::ostream &os, const SymPersymMat& M)
  {
    return os << M.to_str();
  }
};


inline double gaussianKernel(int i, int j, double a, double b)
{
  if(i == j) return 1.0/b;
  return std::pow(a,(i-j)*(i-j))/b;
}

inline SymPersymMat gaussianKernel( size_t n, double a, double b )
{
  SymPersymMat Sigma(n);
  for(size_t i = 0; i < Sigma.n_half(); ++i)
    for( size_t j = i; j < n-i; ++j)
      Sigma(i,j) = gaussianKernel(i,j,a,b);
  return Sigma;
}


inline SymPersymMat gaussianKernelInverse( size_t n, double a, double b )
{ 
  SymPersymMat Omega(n);
  if( n == 1 )
    Omega(0,0) = b;
  if( n <= 1)
    return Omega;
  
  std::vector<double> one_a2k(n-1);
  double a2 = a*a;
  double a_2k = a2;
  for( size_t k = 0; k < n-1; ++k )
  {
    one_a2k[k] = 1 - a_2k;
    a_2k *= a2;
  }
  
  // Compute the (0,0) element of the inverse kernel
  Omega(0,0) = 1.0/one_a2k[0];
  for( size_t k = 1; k < n-1; ++k )
    Omega(0,0) = Omega(0,0)/one_a2k[k];  
  
  // Compute the first row of the inverse kernel
  std::vector<double> Kinv1_Kinv11(n-1);
  for( size_t k = 0; k < n-1; ++k )
  {
    Omega(0,k+1) = -a*one_a2k[n-k-2]/one_a2k[k]*Omega(0,k);
    Kinv1_Kinv11[k] = Omega(0,k+1) / Omega(0,0);
  }
  
  // Initialize the recursion
  // See: M. GOVER, "PROPERTIES OF THE INVERSE OF THE GAUSSIAN MATRIX,"
  //      SIAM J. MATRIX ANAL. APPL. Vol. 12, No. 3, pp. 541-548, July 1991
  std::vector<double> X(n-2);
  for(size_t j = 0; j < n-2; ++j)
    X[j] = one_a2k[n-j-2]*Omega(0,j);
  
  for(size_t i = 1; i < Omega.n_half(); ++i)
  {
    for(size_t j = 0; j < n-2*i-1; ++j)
    {
      double Y = Kinv1_Kinv11[i-1]*Omega(0,j+i);
      Omega(i,j+i) = X[j] + Y;
      X[j] += one_a2k[n-j-2*i-2]*Y; // note: when j = n-2*i-1, the update for X is garbage
    }
    Omega(i,n-i-1) = X[n-2*i-1] + Kinv1_Kinv11[i-1]*Omega(0,n-i-1); 
  }
  
  // Scale the Kernel
  for( size_t i = 0; i < Omega.n_half(); ++i )
    for( size_t j = i; j < n-i; ++j)
      Omega(i,j) *= b;
  return Omega;
}

/*
inline double kroneckerProductIJ( size_t i, size_t j, const std::vector<size_t>& dims, const std::vector<SymPersymMat>& A )
{
  std::vector<size_t> subvi = ind2subv_rowmajor(i, dims.cbegin(), dims.cend());
  std::vector<size_t> subvj = ind2subv_rowmajor(j, dims.cbegin(), dims.cend());
  double Aij = 1.0;
  for(size_t d = 0; d < dims.size(); ++d)
    Aij *= A[d](subvi[d], subvj[d]);
  return Aij;
}
*/

// multiplication of a vector and a matrix built up from a kronecker product
inline std::vector<double> vecKronProd( const std::vector<double>& p,
                                        const std::vector<size_t>& dims,
                                        const std::vector<SymPersymMat>& A )
{
  std::vector<double> q1(p);
  std::vector<double> q2(p.size());
  
  size_t i_left = 1;
  size_t i_right = 1;
  for( size_t d = 1; d<dims.size(); ++d )
    i_right *= dims[d];
  
  for( size_t d = 0; d<dims.size(); ++d )
  {
    size_t base_i = 0;
    size_t base_j = 0;
    for( size_t il = 0; il<i_left; ++il )
    {
      for( size_t ir = 0; ir<i_right; ++ir )
      {
        size_t index_j = base_j + ir;
        for( size_t row = 0; row < dims[d]; ++row )
        {
          size_t index_i = base_i + ir;
          q2[index_j] = A[d](row,0) * q1[index_i];
          index_i += i_right;
          for( size_t col = 1; col < dims[d]; ++col )
          {
            q2[index_j] += A[d](row,col) * q1[index_i];
            index_i += i_right;
          }
          index_j += i_right;
        }
      }
      base_i += dims[d]*i_right;
      base_j += dims[d]*i_right;
    }
    q1 = q2;
    if( d == dims.size()-1 ) break;
    i_left = i_left * dims[d];
    i_right = i_right / dims[d+1];
  }
  return q1; 
}
                                        

// conjugate gradient descent
inline std::vector<double> conjugateGradSolver( const std::vector<size_t>& dims, 
                                                const std::vector<ifom::SymPersymMat> &A,
                                                const std::vector<double> &b,
                                                size_t num_iter )
{
  size_t D = b.size();
  std::vector<double> r(b);
  std::vector<double> x(D,0.0);
  std::vector<double> p,q;
  double alpha, beta, rho, rho_1, pq;
  
  for (size_t k = 0; k < num_iter; ++k)
  {
    rho = 0.0;
    for(size_t i = 0; i < D; ++i)
      rho += r[i]*r[i];
    
    if( k == 0)
      p = r;
    else
    {
      beta = rho / rho_1;
      for(size_t i = 0; i < D; ++i)
        p[i] = r[i] + beta*p[i];
    }
    
    q = vecKronProd( p, dims, A );

    pq = 0.0;
    for(size_t i = 0; i < D; ++i)
      pq += p[i]*q[i];

    alpha = rho / pq;
    for(size_t i = 0; i < D; ++i)
    {
      x[i] += alpha * p[i];
      r[i] -= alpha * q[i];
    }
    
    rho_1 = rho;
  }
  return x;
}



}

#endif
