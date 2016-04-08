#include "libs.h"

// - Use ::vl_rodrigues to compute the Rodrigues formula and its derivative.
// - Use ::vl_irodrigues to compute the inverse Rodrigues formula and
//       its derivative.

/** @brief Compute the minimum between two values
 ** @param x value
 ** @param y value
 ** @return the minimum of @a x and @a y.
 **/
#define MIN(x,y) (((x)<(y))?(x):(y))

/** @brief Compute the maximum between two values
 ** @param x value.
 ** @param y value.
 ** @return the maximum of @a x and @a y.
 **/
#define MAX(x,y) (((x)>(y))?(x):(y))

/** @internal @brief IEEE double precision quiet NaN constant */
static union {long long int raw ; double value; };
const double nan_d = std::numeric_limits<double>::quiet_NaN();
#define NAN_D (nan_d)

class common {
public:
    void rodrigues  (double* R_pt,  double* dR_pt, const double* om_pt);
    void irodrigues (double* om_pt, double* dom_pt, const double* R_pt);
    void print_matrix(double* A);
    void print_matrix_sparse(Eigen::SparseMatrix<double> A);
    void write_vectorXd_tofile(VecXd v);
};
