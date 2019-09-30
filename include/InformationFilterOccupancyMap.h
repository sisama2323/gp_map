#ifndef __INFORMATION_FILTER_OCCUPANCY_MAP_
#define __INFORMATION_FILTER_OCCUPANCY_MAP_

#include "nx_map_utils.h"
#include "SymPersymMat.h"
#include <Eigen/Sparse>
#include <chrono>

inline std::chrono::high_resolution_clock::time_point tic()
{ return std::chrono::high_resolution_clock::now(); }
inline double toc( const std::chrono::high_resolution_clock::time_point& t2)
{ return std::chrono::duration<double>(std::chrono::high_resolution_clock::now()-t2).count(); }

namespace ifom
{

inline double GaussianCDF(double x)
{ 
  // constants
  double a1 =  0.254829592;
  double a2 = -0.284496736;
  double a3 =  1.421413741;
  double a4 = -1.453152027;
  double a5 =  1.061405429;
  double p  =  0.3275911;

  // Save the sign of x
  int sign = 1;
  if (x < 0)
      sign = -1;
  x = std::fabs(x)/std::sqrt(2.0);

  // A&S formula 7.1.26
  double t = 1.0/(1.0 + p*x);
  double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*std::exp(-x*x);
  
  return 0.5*(1.0 + sign*y);
}


// calculate gaussian kernel based on two points
double GaussianKernel(const std::vector<double> kernel_sigma, const std::vector<int> &pt1, const std::vector<int> &pt2){
    int x1 = pt1[0];
    int y1 = pt1[1];
    int z1 = pt1[2];
    int x2 = pt2[0];
    int y2 = pt2[1];
    int z2 = pt2[2];

    int dis = pow((x1 - x2) / kernel_sigma[0], 2) + pow((y1 - y2) / kernel_sigma[1], 2) + pow((z1 - z2) / kernel_sigma[2], 2);
    double det_sigma = std::sqrt(kernel_sigma[0] * kernel_sigma[1] * kernel_sigma[2]);
    double gk = (1 / (15.73 * det_sigma)) * std::exp (-static_cast<double>(dis) / 2.0);
    return gk;
}


class InformationFilterOccupancyMap
{
  // Map dimensions
  size_t num_dim_; // number of dimensions (usually 3 but can be more or fewer)
  size_t half_cell_;
  std::vector<size_t> dims_;
  
  // Inverse Gaussian Kernel 
  std::vector<SymPersymMat> kinv_; // num_dim_ x halfsize_ x (num_cell_, num_cell_-2, num_cell_-4,...)
  std::vector<double> omega0_diag_; // diagonal of the initial information matrix
  double knn_;           // kernel with itself
  std::vector<double> res_;    // map resolution
  std::vector<double> devm_;    // map minimum

  // Observations
  double obs_stdev2_; // squared standard deviation of the observation noise
  double radius_;
  double ad;
  double bd;
  std::vector<SymPersymMat> kmm;
      // compute pseudo point index mask
  int x_idx_range_;
  int y_idx_range_;
  int z_idx_range_;

public:
  std::vector<double> alpha; // difference of 1 and -1 observations per cell
  std::vector<double> obs_num_; // number of observations per cell

  std::vector<std::map<int, double>> beta;

  size_t num_cell_; // total number of cells in the map

  nx::map_nd converter;

  InformationFilterOccupancyMap( const std::vector<double>& devmin,
                                 const std::vector<double>& devmax,
                                 const std::vector<double>& resolution,
                                 const std::vector<double>& kernel_stdev,
                                 double observation_stdev,
                                 double r)
   : num_dim_(devmin.size()), num_cell_(1), dims_(devmin.size()),
     obs_stdev2_(observation_stdev*observation_stdev),
     converter(devmin, devmax, resolution)
  {
    res_ = resolution;       // need resolution for later use
    devm_ = devmin;
    radius_ = r;

    // compute how many closest points to keep
    x_idx_range_ = static_cast<int>(floor(radius_ / res_[0]));
    y_idx_range_ = static_cast<int>(floor(radius_ / res_[1]));
    z_idx_range_ = static_cast<int>(floor(radius_ / res_[2]));

    for( size_t d = 0; d < num_dim_; ++d )
    {
      dims_[d] = static_cast<size_t>(converter.size()[d]);
      ad = std::exp(-1/2.0/kernel_stdev[d]/kernel_stdev[d]);
      bd = sqrt_pi*sqrt_2*kernel_stdev[d];

      num_cell_ *= dims_[d];
      kinv_.push_back(gaussianKernelInverse(dims_[d], ad, bd ));
    }
    half_cell_ = static_cast<size_t>(std::floor((num_cell_+1)/2.0));

    obs_num_.resize(num_cell_,0);
    alpha.resize(num_cell_,0);
    beta.resize(num_cell_);

    // Precompute the diagonal of the inverse Kernel and the squared inverse Kernel
    omega0_diag_.resize(num_cell_,1.0);
    for( size_t i = 0; i < num_cell_; ++i )
    {
      std::vector<int> subvi = nx::ind2subv_rowmajor(i, converter.size().cbegin(), converter.size().cend());
      for(size_t d = 0; d < num_dim_; ++d)
      {
        omega0_diag_[i] *= kinv_[d](subvi[d], subvi[d]);
      }
    }

    // compute kernel with itself
    std::vector<double> pt = {0, 0, 0};
    knn_ = GaussianKernel( kernel_stdev,  point2ind(pt),  point2ind(pt));
  }


  size_t num_cells() const
  { return num_cell_; }

  size_t half_cells() const
  { return half_cell_; }

  // convert coordinate to index
  std::vector<int> point2ind( std::vector<double>& pt_h){

    return std::vector<int>{static_cast<int>(floor((pt_h[0] - devm_[0])/res_[0] + 0.0001)),
                            static_cast<int>(floor((pt_h[1] - devm_[1])/res_[1] + 0.0001)),
                            static_cast<int>(floor((pt_h[2] - devm_[2])/res_[2] + 0.0001))};
  }


  std::vector<double> ind2point( std::vector<int>& pt_w, std::vector<double> r, std::vector<double> minx){
        return std::vector<double>{static_cast<double>(pt_w[0]) * r[0] + minx[0],
                                   static_cast<double>(pt_w[1]) * r[1] + minx[1],
                                   static_cast<double>(pt_w[2]) * r[2] + minx[2]};

  }


  int omega_sub2ind(const int row,const int col) {
      return row*num_cell_+col;
  }


  void omega_ind2sub(const int sub, int &row, int &col) {
      row=sub/num_cell_;
      col=sub%num_cell_;
  }


  double loop_x(const std::vector<std::map<int, double>>& X_1,
                const std::vector<double>& X_2, int i){
    std::map<int, double> x1 = X_1[i];
    std::map<int, double>::iterator it;

    double result = 0;
    for ( it = x1.begin(); it != x1.end(); it++ )
    {
        int index = it->first;  // string (key)
        double value = it->second;   // string's value
        result += value * X_2[index];

    }
    return result;
  }

/*
 * pt_w: new point location in meter, kernel_sigma: kernel width, f: whether this is a free space or not, debug: debug will print stuff
 */
  void addObservation(Eigen::Vector3d& pt_w, const std::vector<double> kernel_sigma, const bool f, const bool debug)
  {
    std::vector<double> pt_w_v = {pt_w(0), pt_w(1), pt_w(2)};

    // get closest pseudo point
    std::vector<double> pt_h_c = converter.meter2meter( pt_w_v );

    // convert eigen to std::vector

    // store all the p point
    std::vector< std::vector<double> > p_point(0);

    // find all closed by pseudo points based on kernel radius
    std::vector<double> pt_h;
    for (int i=-x_idx_range_; i<=x_idx_range_; ++i){
      for (int j=-y_idx_range_; j<=y_idx_range_; ++j) {
        for (int k=-z_idx_range_; k<=z_idx_range_; ++k) {
          pt_h = {static_cast<double> (pt_h_c[0]) + i * res_[0],
                  static_cast<double> (pt_h_c[1]) + j * res_[1],
                  static_cast<double> (pt_h_c[2]) + k * res_[2]};

          // check if pt_h is within the circle
          if (pow(pt_h[0]-pt_w(0), 2) + pow(pt_h[1]-pt_w(1), 2) + pow(pt_h[2]-pt_w(2), 2) <= pow(radius_, 2)) {

              // find pseudo point idx
              int idx = converter.subv2ind_rowmajor(point2ind(pt_h));

              // make sure point fall inside the grid
              if (idx < num_cell_) {

                  double kmn = GaussianKernel(kernel_sigma, point2ind(pt_h), point2ind(pt_w_v));
                  double kmm_inv = omega0_diag_[idx];

                  double z = 1 / (knn_ - kmn * kmm_inv * kmn + obs_stdev2_);

                  if (debug) {
                      std::cout << "kmn " << kmn << std::endl;
                      std::cout << "kmm_inv " << kmm_inv << std::endl;
                      std::cout << "z " << z << std::endl;
                      std::cout << "knn_ " << knn_ << std::endl;
                      std::cout << "denominator " << (knn_ - kmn * kmm_inv * kmn + obs_stdev2_) << std::endl;
                  }

                  std::map<int, double> lb = beta[idx];
                  if (lb.find(idx) == lb.end()){           // doesn't exist
                      lb.insert(std::pair<int, double>(idx, kmn * z * kmn));
                  } else {
                      lb[idx] = lb[idx] + kmn * z * kmn;
                  }

                  // now update the non diag value of beta
                  for (int j=0; j<p_point.size(); j++){
                      std::vector<double> p2 = p_point[j];
                      double z_2 = p2[3];
                      p2.pop_back();
                      double knm_i = GaussianKernel(kernel_sigma, point2ind(pt_h), point2ind(p2));
                      int idx_2 = converter.subv2ind_rowmajor(point2ind(p2));

                      if (idx_2 < num_cell_) {
                        if (lb.find(idx_2) == lb.end()){           // doesn't exist
                            lb.insert(std::pair<int, double>(idx_2, knm_i * z * knm_i));
                        } else {
                            lb[idx_2] = lb[idx_2] + knm_i * z * knm_i;
                        }
                        beta[idx] = lb;

                        std::map<int, double> lb_2 = beta[idx_2];
                        if (lb_2.find(idx) == lb.end()){           // doesn't exist
                            lb_2.insert(std::pair<int, double>(idx, knm_i * z_2 * knm_i));
                        } else {
                            lb_2[idx] = lb_2[idx] + knm_i * z_2 * knm_i;
                        }
                        beta[idx_2] = lb_2;
                      }
                  }
                  beta[idx] = lb;

                  if (f) {
                      if (alpha[idx] < 255) {
                          alpha[idx] += kmn * z;
                      }
                  } else {
                      if (alpha[idx] > -255) {
                          alpha[idx] -= kmn * z;
                      }
                  }

                  pt_h.push_back(z);
                  p_point.push_back(pt_h);

                  if (debug) {
                      std::cout << "alpha " << alpha[idx] << std::endl;
                  }
              }
          }
        }
      }
    }
    if (debug){
        std::cout << "pseudo size " << p_point.size() << std::endl;
    }
  }


  std::vector<double> latentInfoVector() const
  {
    std::vector<double> v(num_cell_);
    for( size_t i = 0; i < num_cell_; ++i )
      v[i] = static_cast<double>(alpha[i]) / obs_stdev2_;
    return v;
  }

  std::vector<double> deltaInfoMat() const
  {
    std::vector<double> d(num_cell_);
    for( size_t i = 0; i < num_cell_; ++i )
      d[i] = static_cast<double>(obs_num_[i]) / obs_stdev2_;
    return d;
  }


  std::vector<double> latentMean( size_t num_iter )
  {
    std::vector<double> x(num_cell_,0.0);

    // Initialize residual as the information mean
    std::vector<double> r(num_cell_);
    r = alpha;

    // compute kernel
    for( size_t d = 0; d < num_dim_; ++d ) {
        kmm.push_back(gaussianKernel( dims_[d], ad, bd ));
    };

    std::vector<double> p(r);
    std::vector<double> q;
    double al, be, rho, rho_1, pq;
    for (size_t k = 0; k < num_iter; ++k)
    {
      rho = 0.0;
      for(size_t i = 0; i < num_cell_; ++i)
        rho += r[i]*r[i];

      if( k > 0)
      {
        be = rho / rho_1;
        for(size_t i = 0; i < num_cell_; ++i)
          p[i] = r[i] + be*p[i];
      }

      // Main computation
      q = vecKronProd( p, dims_, kmm );
      for(size_t i = 0; i < num_cell_; ++i)
        q[i] +=  loop_x(beta, p, i);

      pq = 0.0;
      for(size_t i = 0; i < num_cell_; ++i)
        pq += p[i]*q[i];

      al = rho / pq;
      for(size_t i = 0; i < num_cell_; ++i)
      {
        x[i] += al * p[i];
        r[i] -= al * q[i];
      }

      rho_1 = rho;

      double sum_of_elems = 0;
      for (auto& n : r)
          sum_of_elems += n;

      std::cout << "iter " << k << " residue: " << sum_of_elems << std::endl;

    }
    return x;
  }


  std::vector<double> latentCov() const
  {
    std::vector<double> CovDiag(num_cell_);
    for( size_t i = 0; i < num_cell_; ++i )  {

      std::map<int, double> x1 = beta[i];
      std::map<int, double>::iterator it;
      double denom = 0;
      for ( it = x1.begin(); it != x1.end(); it++ ) {
          int index = it->first;  // string (key)
          double value = it->second;   // string's value

          double kmm_i = 1;
          std::vector<int> subvi_col = nx::ind2subv_rowmajor(index, converter.size().cbegin(), converter.size().cend());
          std::vector<int> subvi_row = nx::ind2subv_rowmajor(i, converter.size().cbegin(), converter.size().cend());
          for(size_t d = 0; d < num_dim_; ++d)
          {
            kmm_i *= kmm[d](subvi_row[d], subvi_col[d]);
          }

          denom += value + kmm_i;
      }

      CovDiag[i] = (omega0_diag_[i] + x1[i]) / denom;

    }
    return CovDiag;
  }


  std::vector<double> latentMeanfromNum( size_t num_iter, std::vector<int> obs_num_tt, std::vector<int> obs_diff_tt ) const
  {
    std::vector<double> x(num_cell_,0.0);

    // Initialize residual as the information mean
    std::vector<double> r(num_cell_);
    for( size_t i = 0; i < num_cell_; ++i )
      r[i] = static_cast<double>(obs_diff_tt[i]) / obs_stdev2_;

    // Compute the matrix containing the number of observations
    std::vector<double> D(num_cell_);
    for( size_t i = 0; i < num_cell_; ++i )
      D[i] = static_cast<double>(obs_num_tt[i]) / obs_stdev2_;

    std::vector<double> p(r);
    std::vector<double> q;
    double al, be, rho, rho_1, pq;
    for (size_t k = 0; k < num_iter; ++k)
    {
      rho = 0.0;
      for(size_t i = 0; i < num_cell_; ++i)
        rho += r[i]*r[i];

      if( k > 0)
      {
        be = rho / rho_1;
        for(size_t i = 0; i < num_cell_; ++i)
          p[i] = r[i] + be*p[i];
      }

      // Main computation
      q = vecKronProd( p, dims_, kinv_ );
      for(size_t i = 0; i < num_cell_; ++i)
        q[i] +=  D[i]*p[i];

      pq = 0.0;
      for(size_t i = 0; i < num_cell_; ++i)
        pq += p[i]*q[i];

      al = rho / pq;
      for(size_t i = 0; i < num_cell_; ++i)
      {
        x[i] += al * p[i];
        r[i] -= al * q[i];
      }
      
      rho_1 = rho;
    }
    return x;    
  }

    
  std::vector<double> mapPdf( const std::vector<double>& Mu,
                              const std::vector<double>& CovDiag )
  {
    std::vector<double> mpdf(num_cell_);
    for( size_t i = 0; i < num_cell_; ++i )
      mpdf[i] = GaussianCDF(Mu[i] / std::sqrt(CovDiag[i]+obs_stdev2_));

    return mpdf;
  }
};

}

#endif
