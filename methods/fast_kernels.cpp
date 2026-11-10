#include <Rcpp.h>
using namespace Rcpp;

// Sum active rows, processing four columns together without changing row order.
static NumericVector colsums_times(const NumericMatrix& Z, const std::vector<double>& v) {
    const int n = Z.nrow(), d = Z.ncol();
    std::vector<int> act;
    act.reserve(n);
    for (int i = 0; i < n; ++i)
        if (v[i] != 0.0)
            act.push_back(i);
    const int m = act.size();
    NumericVector out(d);
    const double* z = Z.begin();
    int j = 0;
    for (; j + 3 < d; j += 4) {
        const double* z0 = z + (size_t)j * n;
        const double* z1 = z0 + n;
        const double* z2 = z1 + n;
        const double* z3 = z2 + n;
        long double s0 = 0.0L, s1 = 0.0L, s2 = 0.0L, s3 = 0.0L;
        for (int k = 0; k < m; ++k) {
            const int i = act[k];
            const double vi = v[i];
            double p0 = z0[i] * vi, p1 = z1[i] * vi, p2 = z2[i] * vi, p3 = z3[i] * vi;
            s0 += p0;
            s1 += p1;
            s2 += p2;
            s3 += p3;
        }
        out[j] = (double)s0;
        out[j + 1] = (double)s1;
        out[j + 2] = (double)s2;
        out[j + 3] = (double)s3;
    }
    for (; j < d; ++j) {
        const double* zj = z + (size_t)j * n;
        long double s = 0.0L;
        for (int k = 0; k < m; ++k) {
            const int i = act[k];
            double p = zj[i] * v[i];
            s += p;
        }
        out[j] = (double)s;
    }
    return out;
}

// No fast-math: preserve product order and long-double accumulation.
// [[Rcpp::export]]
NumericVector itr_gradient(NumericMatrix Z, NumericVector weight,
                           NumericVector label, NumericVector margin, bool concave) {
    int n = Z.nrow();
    NumericVector multiplier(n);
    for (int i = 0; i < n; ++i) {
        double u = margin[i];
        double g = concave ? (u >= 0 ? 0 : (u >= -1 ? 2 * u : -2))
                           : (u >= 1 ? 0 : (u >= 0 ? -2 * (1 - u) : -2));
        multiplier[i] = (weight[i] * g) * label[i];
    }
    std::vector<double> values(multiplier.begin(), multiplier.end());
    return colsums_times(Z, values);
}

// [[Rcpp::export]]
double itr_loss(NumericVector weight, NumericVector margin, bool ramp) {
    long double s = 0;
    for (int i = 0; i < margin.size(); ++i) {
        double u = margin[i];
        double a = u >= 1 ? 0 : (u >= 0 ? (1 - u) * (1 - u) : 1 - 2 * u);
        if (ramp) {
            double b = u >= 0 ? 0 : (u >= -1 ? u * u : -2 * u - 1);
            a = 0.5 * (a - b);
        }
        double product = weight[i] * a;
        s += product;
    }
    return static_cast<double>(s);
}
