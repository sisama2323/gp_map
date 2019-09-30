#ifndef _OWENS_HPP_INCLUDED
#define _OWENS_HPP_INCLUDED

namespace ifom
{

void bivariate_normal_cdf_values ( int &n_data, double &x, double &y,
  double &r, double &fxy );
double bivnor ( double h, double k, double r );
double gauss ( double t );
void normal_01_cdf_values ( int &n_data, double &x, double &fx );
void owen_values ( int &n_data, double &h, double &a, double &t );
double q ( double h, double ah );
double r8_abs ( double x );
double r8_max ( double x, double y );
double r8_min ( double x, double y );
double t_compute ( double h, double a );
double tfun ( double h, double a, double ah );
void timestamp ( void );
double znorm1 ( double z );
double znorm2 ( double z );

}

#endif
