
#include <Rcpp.h>
using namespace Rcpp;


Rcpp::NumericVector john_powercpp(Rcpp::NumericVector x){
  Rcpp::LogicalVector between = (x != 1.0) & ( x != 0.0);
  Rcpp::NumericVector out(x.size());
  std::transform(x.begin(), x.end(), between.begin(), out.begin(), ::pow);
  return out;
}


// grid / mouse position in mux/y/z
// data in x/y/z
double normal_kernel_3d_indicatorcpp(
    Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::NumericVector z,  
    double mux, double muy, double muz,
    double sd ) {
  Rcpp::NumericVector temp_x, temp_y, temp_z;
  temp_x  = john_powercpp(Rcpp::dnorm(x, mux, sd));
  temp_y  = john_powercpp(Rcpp::dnorm(y, muy, sd));
  temp_z  = john_powercpp(Rcpp::dnorm(z, muz, sd));
  double temp = mean(temp_x* temp_y* temp_z);
  return temp;
}



// function to loop through vector for `x` and `mu` in dnorm
// i.e. calculate normal_kernel_3d_indicatorcpp in a grid
//' Calculates 3D Gaussian intensity over a grid.
//'
//' @param N An integer that defines where the intensity will be calculated in the unit cube into. If N=10 then 10^3 cubes are formed.
//' @param x,y,z A vector of positions on one dimension to calculate the kernel on.
//' @param sd The standard deviation which is assumed to be the same for each dimension.
//' @return Returns a vector giving the intensity at each grid point.
//' @examples
//' set.seed(1)
//' N = 2
//' x = rnorm(5); mux=1
//' y = rnorm(5); muy=1
//' z = rnorm(5); muz=1
//' sd=1
//' normal_kernel_3d_indicator_vectorcpp(N, x, y, z, sd)
//' @export
//' @useDynLib jqmpp
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::NumericVector normal_kernel_3d_indicator_vectorcpp(
    int N,
    Rcpp::NumericVector x, Rcpp::NumericVector y, Rcpp::NumericVector z,  
    double sd=0.1 ) {
  N = N + 1;
  double reps = pow(N, 3);
  Rcpp::NumericVector intensity(reps);

  Rcpp::IntegerVector iv = Rcpp::seq(0, N-1);
  Rcpp::NumericVector v = as<Rcpp::NumericVector>(iv)/(N-1);
  
  for(int i=0; i < N; i++){
    for(int j=0; j < N; j++){
      for(int k=0; k < N; k++){
        intensity[(i * N + j) * N + k] = normal_kernel_3d_indicatorcpp(x, y, z, v(i), v(j), v(k), sd);
      }
    }
  }
  return intensity;
}




//' Calculates distance between two matrices with the same number of columns.
//'
//' @param m1,m2 Two matrices with the same number of columns.
//' @return Returns a matrix of distances with dimension nrow(m1) by nrow(m2).
//' @details The distance is calculated between each row of the two input matrices.
//' @examples
//' my_dist(diag(2), diag(1, 3, 2))
//' @export
//' @useDynLib jqmpp
//' @importFrom Rcpp sourceCpp
// [[Rcpp::export]]
Rcpp::NumericMatrix my_dist(Rcpp::NumericMatrix m1,  Rcpp::NumericMatrix m2) {
  
  int nr = m1.nrow();
  int nc = m2.nrow();
  Rcpp::NumericMatrix dmat(nr, nc);
  
  for(int i=0; i < nr; i++){
    for(int j=0; j < nc; j++){
      dmat(i, j) = sqrt(sum(pow(m1.row(i) - m2.row(j), 2.0)));
    }
  }
  return dmat;
}


