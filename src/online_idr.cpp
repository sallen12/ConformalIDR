#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

// pd_ds headers
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;

// useful definitions and functions --------------------------------------------
// ordered set
typedef tree<double, null_type, less<double>, rb_tree_tag,
             tree_order_statistics_node_update>
  ordered_set;

// computation with sequential update ------------------------------------------

// [[Rcpp::export]]
List cidr_sequential(
    NumericVector x_r,
    NumericVector x_out,
    List w,
    List W,
    List w_out,
    List pos_x,
    NumericVector y_unique_r,
    NumericVector y_out,
    int n_thr,
    int n_x) {
  int n = x_r.size();
  int n_out = x_out.size();
  int n_y = y_unique_r.size();
  int n_prt = 3 * (pow(n_x, 2.0 / 3.0) + 1) + 1;
  int n_psx = pos_x.length();


  vector<double> x;
  x.reserve(n);
  for (int i = 0; i < n; i++) {
    x.push_back(x_r[i]);
  }
  ordered_set y_unique;
  for (int i = 0; i < n_y; i++) {
    y_unique.insert(y_unique_r[i]);
  }


  int j0, s0, a0, b0, c, n_px, lwr, upr;
  int pos_new_x = 1;
  int pos_new_y;
  double w_lwr, w_upr;
  bool copy, update_pos;
  vector<int>::iterator pos_int;
  vector<double>::iterator pos_dbl;
  IntegerVector psx;
  List W_tmp, W_update;
  NumericVector w_tmp, w_out_tmp, W_tmp_entries;
  NumericVector w_update, w_out_update, W_update_entries;
  NumericVector z(n_x);
  std::vector<int> pp;
  pp.reserve(n_prt);
  std::vector<double> ww;
  ww.reserve(n_prt);
  std::vector<double> mm;
  mm.reserve(n_prt);
  std::vector<int> rem_pp;
  rem_pp.reserve(n_prt);
  std::vector<double> rem_ww;
  rem_ww.reserve(n_prt);
  std::vector<double> rem_mm;
  rem_mm.reserve(n_prt);

  pp.push_back(-1);
  ww.push_back(1.0);
  mm.push_back(10.0);

  NumericMatrix cdf(n_thr + 2, n_x);
  NumericMatrix cdf_lwr(n_thr + 2, n_out);
  NumericMatrix cdf_upr(n_thr + 2, n_out);
  NumericMatrix cdf_oos(n_thr + 2, n_out);
  NumericMatrix cdf_is(n_thr + 2, n_out);
  NumericMatrix cdf_lcnf(n_thr + 2, n_out);
  NumericMatrix cdf_ucnf(n_thr + 2, n_out);
  NumericMatrix points(n_thr + 2, n_out);

  for (int i = 0; i < n_out; i++) {
    cdf_lwr(0, i) = 0.0;
    cdf_upr(0, i) = 0.0;
    cdf_oos(0, i) = 0.0;
    cdf_is(0, i) = 0.0;
    points(0, i) = R_NegInf;
    for (int j = 1; j < n_thr + 2; j++) {
      points(j, i) = R_PosInf;
      cdf_lwr(j, i) = 1.0;
      cdf_upr(j, i) = 1.0;
      cdf_oos(j, i) = 1.0;
      cdf_is(j, i) = 1.0;
      points(j, i) = R_PosInf;
    }
  }
  for (int i = 0; i < n_x; i++) {
    cdf(0, i) = 0.0;
    for (int j = 1; j < n_thr + 2; j++) {
      cdf(j, i) = 1.0;
    }
  }

  for (int i = 0; i < n_out; i++) {
    W_tmp = W[i];
    w_tmp = w[i];
    w_out_tmp = w_out[i];
    c = 1;
    for (ordered_set::iterator it = y_unique.begin(); it != y_unique.end(); it++) {
      points(c, i) = *it;
      c = c + 1;
    }

    update_pos = false;
    n = x.size();
    pp.erase(pp.begin() + 1, pp.end());
    ww.erase(ww.begin() + 1, ww.end());
    mm.erase(mm.begin() + 1, mm.end());
    pp.push_back(n - 1);
    ww.push_back(sum(w_tmp));
    mm.push_back(0.0);
    z.fill(0.0);
    n_psx = pos_x.size();

    for (int k = 0; k < n_psx; k++) {
      psx = pos_x[k];
      W_tmp_entries = W_tmp[k];
      n_px = psx.length();
      for (int j = 0; j < n_px; j++) {
        copy = false;
        j0 = psx[j];

        z[j0] = z[j0] + W_tmp_entries[j] / w_tmp[j0];
        pos_int = std::lower_bound(pp.begin(), pp.end(), j0);
        s0 = std::distance(pp.begin(), pos_int);
        a0 = pp[s0 - 1] + 1;
        b0 = pp[s0];

        std::vector<int> rem_pp;
        std::vector<double> rem_ww;
        std::vector<double> rem_mm;

        if (s0 < pp.size() - 1) {
          copy = true;
          for (int l = s0 + 1; l < pp.size(); l++) {
            rem_pp.push_back(pp[l]);
            rem_ww.push_back(ww[l]);
            rem_mm.push_back(mm[l]);
          }
        }

        pp.erase(pp.begin() + s0, pp.end());
        ww.erase(ww.begin() + s0, ww.end());
        mm.erase(mm.begin() + s0, mm.end());
        pp.push_back(j0);
        ww.push_back(sum(w_tmp[Range(a0, j0)]));
        mm.push_back(sum(w_tmp[Range(a0, j0)] * z[Range(a0, j0)]) / *(ww.end() - 1));

        while (*(mm.end() - 1) >= *(mm.end() - 2)) {
          pp[pp.size() - 2] = *(pp.end() - 1);
          mm[mm.size() - 2] = *(ww.end() - 1) * (*(mm.end() - 1)) +
            *(ww.end() - 2) * (*(mm.end() - 2));
            ww[ww.size() - 2] = *(ww.end() - 1) + *(ww.end() - 2);
            pp.pop_back();
            ww.pop_back();
            mm.pop_back();
            mm[mm.size() - 1] = *(mm.end() - 1) / *(ww.end() - 1);
        }

        if (j0 < b0) {
          for (int l = j0 + 1; l <= b0; l++) {
            pp.push_back(l);
            ww.push_back(w_tmp[l]);
            mm.push_back(z[l]);
            while (*(mm.end() - 1) >= *(mm.end() - 2)) {
              pp[pp.size() - 2] = l;
              mm[mm.size() - 2] = *(ww.end() - 1) * (*(mm.end() - 1)) +
                *(ww.end() - 2) * (*(mm.end() - 2));
                ww[ww.size() - 2] = *(ww.end() - 1) + *(ww.end() - 2);
                pp.pop_back();
                ww.pop_back();
                mm.pop_back();
                mm[mm.size() - 1] = *(mm.end() - 1) / *(ww.end() - 1);
            }
          }
        }

        if (copy) {
          pp.insert(pp.end(), rem_pp.begin(), rem_pp.end());
          ww.insert(ww.end(), rem_ww.begin(), rem_ww.end());
          mm.insert(mm.end(), rem_mm.begin(), rem_mm.end());
        }
        rem_pp.clear();
        rem_ww.clear();
        rem_mm.clear();
      }

      for (int l = 1; l < pp.size(); l++) {
        a0 = pp[l - 1] + 1;
        b0 = pp[l];
        for (int r = a0; r <= b0; r++) {
          cdf(k + 1, r) = mm[l];
        }
      }
    }

    pos_dbl = std::lower_bound(x.begin(), x.end(), x_out[i]);
    pos_new_x = std::distance(x.begin(), pos_dbl);
    if (pos_new_x == x.size()) {
      w_lwr = 1.0;
      w_upr = 0.0;
      lwr = x.size() - 1;
      upr = x.size() - 1;
      x.push_back(x_out[i]);
      w_tmp.push_back(w_out_tmp[i]);
      for (int r = i + 1; r < n_out; r++) {
        w_update = w[r];
        w_out_update = w_out[r];
        w_update.push_back(w_out_update[i]);
        w[r] = w_update;
      }
      update_pos = false;
    } else if (*pos_dbl == x_out[i]) {
      lwr = pos_new_x;
      upr = pos_new_x;
      w_lwr = 1.0;
      w_upr = 0.0;
      update_pos = false;
      w_tmp[pos_new_x] = w_tmp[pos_new_x] + w_out_tmp[i];
      for (int r = i + 1; r < n_out; r++) {
        w_update = w[r];
        w_out_update = w_out[r];
        w_update[pos_new_x] = w_update[pos_new_x] + w_out_update[i];
        w[r] = w_update;
      }
    } else if (pos_new_x == 0) {
      upr = 0;
      lwr = 0;
      w_upr = 1.0;
      w_lwr = 0.0;
      x.insert(x.begin(), x_out[i]);
      w_tmp.insert(w_tmp.begin(), w_out_tmp[i]);
      update_pos = true;
      for (int r = i + 1; r < n_out; r++) {
        w_update = w[r];
        w_out_update = w_out[r];
        w_update.insert(w_update.begin(), w_out_update[i]);
        w[r] = w_update;
      }
    } else {
      upr = pos_new_x;
      lwr = upr - 1;
      w_lwr = (x[upr] - x_out[i]) / (x[upr] - x[lwr]);
      w_upr = 1.0 - w_lwr;
      x.insert(pos_dbl, x_out[i]);
      w_tmp.insert(w_tmp.begin() + pos_new_x, w_out_tmp[i]);
      update_pos = true;
      for (int r = i + 1; r < n_out; r++) {
        w_update = w[r];
        w_out_update = w_out[r];
        w_update.insert(w_update.begin() + pos_new_x, w_out_tmp[i]);
        w[r] = w_update;
      }
    }

    for (int l = 0; l < n_thr + 2; l++) {
      cdf_lwr(l, i) = cdf(l, lwr);
      cdf_upr(l, i) = cdf(l, upr);
      cdf_oos(l, i) = w_lwr * cdf_lwr(l, i) + w_upr * cdf_upr(l, i);
    }

    if (update_pos) {
      for (int l = 0; l < pos_x.size(); l++) {
        psx = pos_x[l];
        for (int r = 0; r < psx.size(); r++) {
          if (pos_new_x <= psx[r]) {
            psx[r] = psx[r] + 1;
          }
        }
        pos_x[l] = psx;
      }
    }

    n = x.size();
    pp.erase(pp.begin() + 1, pp.end());
    ww.erase(ww.begin() + 1, ww.end());
    mm.erase(mm.begin() + 1, mm.end());
    z.fill(0.0);

    z[pos_new_x] = w_out_tmp[i] / w_tmp[pos_new_x];
    pp.push_back(pos_new_x);
    ww.push_back(sum(w_tmp[Range(0, pos_new_x)]));
    mm.push_back(w_out_tmp[i] / *(ww.end() - 1));
    for (int l = 0; l < pos_new_x + 1; l++) {
      cdf(0, l) = *(mm.end() - 1);
    }
    if (pos_new_x < n - 1) {
      pp.push_back(n - 1);
      ww.push_back(sum(w_tmp[Range(pos_new_x + 1, n - 1)]));
      mm.push_back(0.0);
    }

    for (int k = 0; k < n_psx; k++) {
      psx = pos_x[k];
      W_tmp_entries = W_tmp[k];
      n_px = psx.length();
      for (int j = 0; j < n_px; j++) {
        copy = false;
        j0 = psx[j];

        z[j0] = z[j0] + W_tmp_entries[j] / w_tmp[j0];
        pos_int = std::lower_bound(pp.begin(), pp.end(), j0);
        s0 = std::distance(pp.begin(), pos_int);
        a0 = pp[s0 - 1] + 1;
        b0 = pp[s0];

        std::vector<int> rem_pp;
        std::vector<double> rem_ww;
        std::vector<double> rem_mm;

        if (s0 < pp.size() - 1) {
          copy = true;
          for (int l = s0 + 1; l < pp.size(); l++) {
            rem_pp.push_back(pp[l]);
            rem_ww.push_back(ww[l]);
            rem_mm.push_back(mm[l]);
          }
        }

        pp.erase(pp.begin() + s0, pp.end());
        ww.erase(ww.begin() + s0, ww.end());
        mm.erase(mm.begin() + s0, mm.end());
        pp.push_back(j0);
        ww.push_back(sum(w_tmp[Range(a0, j0)]));
        mm.push_back(sum(w_tmp[Range(a0, j0)] * z[Range(a0, j0)]) / *(ww.end() - 1));

        while (*(mm.end() - 1) >= *(mm.end() - 2)) {
          pp[pp.size() - 2] = *(pp.end() - 1);
          mm[mm.size() - 2] = *(ww.end() - 1) * (*(mm.end() - 1)) +
            *(ww.end() - 2) * (*(mm.end() - 2));
            ww[ww.size() - 2] = *(ww.end() - 1) + *(ww.end() - 2);
            pp.pop_back();
            ww.pop_back();
            mm.pop_back();
            mm[mm.size() - 1] = *(mm.end() - 1) / *(ww.end() - 1);
        }

        if (j0 < b0) {
          for (int l = j0 + 1; l <= b0; l++) {
            pp.push_back(l);
            ww.push_back(w_tmp[l]);
            mm.push_back(z[l]);
            while (*(mm.end() - 1) >= *(mm.end() - 2)) {
              pp[pp.size() - 2] = l;
              mm[mm.size() - 2] = *(ww.end() - 1) * (*(mm.end() - 1)) +
                *(ww.end() - 2) * (*(mm.end() - 2));
                ww[ww.size() - 2] = *(ww.end() - 1) + *(ww.end() - 2);
                pp.pop_back();
                ww.pop_back();
                mm.pop_back();
                mm[mm.size() - 1] = *(mm.end() - 1) / *(ww.end() - 1);
            }
          }
        }

        if (copy) {
          pp.insert(pp.end(), rem_pp.begin(), rem_pp.end());
          ww.insert(ww.end(), rem_ww.begin(), rem_ww.end());
          mm.insert(mm.end(), rem_mm.begin(), rem_mm.end());
        }
        rem_pp.clear();
        rem_ww.clear();
        rem_mm.clear();
      }

      for (int l = 1; l < pp.size(); l++) {
        a0 = pp[l - 1] + 1;
        b0 = pp[l];
        for (int r = a0; r <= b0; r++) {
          cdf(k + 1, r) = mm[l];
        }
      }
    }

    for (int l = 0; l < n_thr + 2; l++) {
      cdf_lcnf(l, i) = cdf(l, pos_new_x);
    }

    for (int l = 0; l < pos_new_x + 1; l++) {
      cdf(0, l) = 0.0;
    }

    pp.erase(pp.begin() + 1, pp.end());
    ww.erase(ww.begin() + 1, ww.end());
    mm.erase(mm.begin() + 1, mm.end());
    pp.push_back(n - 1);
    ww.push_back(sum(w_tmp));
    mm.push_back(0.0);
    z.fill(0.0);

    for (int k = 0; k < n_psx; k++) {
      psx = pos_x[k];
      W_tmp_entries = W_tmp[k];
      n_px = psx.length();
      for (int j = 0; j < n_px; j++) {
        copy = false;
        j0 = psx[j];

        z[j0] = z[j0] + W_tmp_entries[j] / w_tmp[j0];
        pos_int = std::lower_bound(pp.begin(), pp.end(), j0);
        s0 = std::distance(pp.begin(), pos_int);
        a0 = pp[s0 - 1] + 1;
        b0 = pp[s0];

        std::vector<int> rem_pp;
        std::vector<double> rem_ww;
        std::vector<double> rem_mm;

        if (s0 < pp.size() - 1) {
          copy = true;
          for (int l = s0 + 1; l < pp.size(); l++) {
            rem_pp.push_back(pp[l]);
            rem_ww.push_back(ww[l]);
            rem_mm.push_back(mm[l]);
          }
        }

        pp.erase(pp.begin() + s0, pp.end());
        ww.erase(ww.begin() + s0, ww.end());
        mm.erase(mm.begin() + s0, mm.end());
        pp.push_back(j0);
        ww.push_back(sum(w_tmp[Range(a0, j0)]));
        mm.push_back(sum(w_tmp[Range(a0, j0)] * z[Range(a0, j0)]) / *(ww.end() - 1));

        while (*(mm.end() - 1) >= *(mm.end() - 2)) {
          pp[pp.size() - 2] = *(pp.end() - 1);
          mm[mm.size() - 2] = *(ww.end() - 1) * (*(mm.end() - 1)) +
            *(ww.end() - 2) * (*(mm.end() - 2));
            ww[ww.size() - 2] = *(ww.end() - 1) + *(ww.end() - 2);
            pp.pop_back();
            ww.pop_back();
            mm.pop_back();
            mm[mm.size() - 1] = *(mm.end() - 1) / *(ww.end() - 1);
        }

        if (j0 < b0) {
          for (int l = j0 + 1; l <= b0; l++) {
            pp.push_back(l);
            ww.push_back(w_tmp[l]);
            mm.push_back(z[l]);
            while (*(mm.end() - 1) >= *(mm.end() - 2)) {
              pp[pp.size() - 2] = l;
              mm[mm.size() - 2] = *(ww.end() - 1) * (*(mm.end() - 1)) +
                *(ww.end() - 2) * (*(mm.end() - 2));
                ww[ww.size() - 2] = *(ww.end() - 1) + *(ww.end() - 2);
                pp.pop_back();
                ww.pop_back();
                mm.pop_back();
                mm[mm.size() - 1] = *(mm.end() - 1) / *(ww.end() - 1);
            }
          }
        }

        if (copy) {
          pp.insert(pp.end(), rem_pp.begin(), rem_pp.end());
          ww.insert(ww.end(), rem_ww.begin(), rem_ww.end());
          mm.insert(mm.end(), rem_mm.begin(), rem_mm.end());
        }
        rem_pp.clear();
        rem_ww.clear();
        rem_mm.clear();
      }

      for (int l = 1; l < pp.size(); l++) {
        a0 = pp[l - 1] + 1;
        b0 = pp[l];
        for (int r = a0; r <= b0; r++) {
          cdf(k + 1, r) = mm[l];
        }
      }
    }

    for (int l = 0; l < n_thr + 2; l++) {
      cdf_ucnf(l, i) = cdf(l, pos_new_x);
    }

    pos_new_y = y_unique.order_of_key(y_out[i]);
    if (y_unique.find(y_out[i]) == y_unique.end()) {
      y_unique.insert(y_out[i]);
      psx = {pos_new_x};
      pos_x.insert(pos_x.begin() + pos_new_y, psx);
      for (int r = i + 1; r < n_out; r++) {
        W_update = W[r];
        w_out_update = w_out[r];
        W_update.insert(W_update.begin() + pos_new_y, w_out_update[i]);
        W[r] = W_update;
      }
    } else {
      psx = pos_x[pos_new_y];
      psx.push_back(pos_new_x);
      pos_x[pos_new_y] = psx;
      for (int r = i + 1; r < n_out; r++) {
        W_update = W[r];
        W_update_entries = W_update[pos_new_y];
        w_out_update = w_out[r];
        W_update_entries.push_back(w_out_update[i]);
        W[r] = W_update;
      }
    }
    Rcpp::checkUserInterrupt();
  }

  return List::create(
    _["points"] = points,
    _["cdf_lwr"] = cdf_lwr,
    _["cdf_upr"] = cdf_upr,
    _["cdf_oos"] = cdf_oos,
    _["cdf_lcnf"] = cdf_lcnf,
    _["cdf_ucnf"] = cdf_ucnf
  );
}

// computation with static training data set -----------------------------------


// [[Rcpp::export]]
List cidr_static(
    NumericVector x_r,
    NumericVector x_out,
    NumericVector w,
    List W,
    double w_out,
    List pos_x,
    NumericVector y_unique_r,
    NumericVector y_out,
    int n_thr,
    int n_x) {
  int n = x_r.size();
  int n_out = x_out.size();
  int n_y = y_unique_r.size();
  int n_prt = 3 * (pow(n_x, 2.0 / 3.0) + 1) + 1;
  int n_psx = pos_x.length();

  vector<double> x;
  x.reserve(n);
  for (int i = 0; i < n; i++) {
    x.push_back(x_r[i]);
  }
  ordered_set y_unique;
  for (int i = 0; i < n_y; i++) {
    y_unique.insert(y_unique_r[i]);
  }

  int j0, s0, a0, b0, c, n_px, lwr, upr;
  double w_lwr, w_upr;
  bool copy;
  vector<int> pos_new_x;
  vector<int>::iterator pos_int;
  vector<double>::iterator pos_dbl;
  IntegerVector psx;
  NumericVector W_tmp;
  NumericVector z(n_x + 1);
  std::vector<int> pp;
  pp.reserve(n_prt);
  std::vector<double> ww;
  ww.reserve(n_prt);
  std::vector<double> mm;
  mm.reserve(n_prt);
  std::vector<int> rem_pp;
  rem_pp.reserve(n_prt);
  std::vector<double> rem_ww;
  rem_ww.reserve(n_prt);
  std::vector<double> rem_mm;
  rem_mm.reserve(n_prt);

  pp.push_back(-1);
  ww.push_back(1.0);
  mm.push_back(10.0);

  NumericMatrix cdf(n_thr + 2, n_x + 1);
  NumericMatrix cdf_lwr(n_thr + 2, n_out);
  NumericMatrix cdf_upr(n_thr + 2, n_out);
  NumericMatrix cdf_oos(n_thr + 2, n_out);
  NumericMatrix cdf_is(n_thr + 2, n_out);
  NumericMatrix cdf_lcnf(n_thr + 2, n_out);
  NumericMatrix cdf_ucnf(n_thr + 2, n_out);
  NumericVector points = y_unique_r;
  points.insert(points.begin(), R_NegInf);
  points.push_back(R_PosInf);
  for (int i = 0; i < n_out; i++) {
    cdf_lwr(0, i) = 0.0;
    cdf_upr(0, i) = 0.0;
    cdf_oos(0, i) = 0.0;
    cdf_is(0, i) = 0.0;
    for (int j = 1; j < n_thr + 2; j++) {
      cdf_lwr(j, i) = 1.0;
      cdf_upr(j, i) = 1.0;
      cdf_oos(j, i) = 1.0;
      cdf_is(j, i) = 1.0;
    }
  }
  for (int i = 0; i < n_x + 1; i++) {
    cdf(0, i) = 0.0;
    for (int j = 1; j < n_thr + 2; j++) {
      cdf(j, i) = 1.0;
    }
  }

  n = x.size();
  pp.erase(pp.begin() + 1, pp.end());
  ww.erase(ww.begin() + 1, ww.end());
  mm.erase(mm.begin() + 1, mm.end());
  pp.push_back(n - 1);
  ww.push_back(sum(w));
  mm.push_back(0.0);
  z.fill(0.0);
  n_psx = pos_x.size();

  for (int k = 0; k < n_psx; k++) {
    psx = pos_x[k];
    W_tmp = W[k];
    n_px = psx.length();
    for (int j = 0; j < n_px; j++) {
      copy = false;
      j0 = psx[j];

      z[j0] = z[j0] + W_tmp[j] / w[j0];
      pos_int = std::lower_bound(pp.begin(), pp.end(), j0);
      s0 = std::distance(pp.begin(), pos_int);
      a0 = pp[s0 - 1] + 1;
      b0 = pp[s0];

      std::vector<int> rem_pp;
      std::vector<double> rem_ww;
      std::vector<double> rem_mm;

      if (s0 < pp.size() - 1) {
        copy = true;
        for (int l = s0 + 1; l < pp.size(); l++) {
          rem_pp.push_back(pp[l]);
          rem_ww.push_back(ww[l]);
          rem_mm.push_back(mm[l]);
        }
      }

      pp.erase(pp.begin() + s0, pp.end());
      ww.erase(ww.begin() + s0, ww.end());
      mm.erase(mm.begin() + s0, mm.end());
      pp.push_back(j0);
      ww.push_back(sum(w[Range(a0, j0)]));
      mm.push_back(sum(w[Range(a0, j0)] * z[Range(a0, j0)]) / *(ww.end() - 1));

      while (*(mm.end() - 1) >= *(mm.end() - 2)) {
        pp[pp.size() - 2] = *(pp.end() - 1);
        mm[mm.size() - 2] = *(ww.end() - 1) * (*(mm.end() - 1)) +
          *(ww.end() - 2) * (*(mm.end() - 2));
          ww[ww.size() - 2] = *(ww.end() - 1) + *(ww.end() - 2);
          pp.pop_back();
          ww.pop_back();
          mm.pop_back();
          mm[mm.size() - 1] = *(mm.end() - 1) / *(ww.end() - 1);
      }

      if (j0 < b0) {
        for (int l = j0 + 1; l <= b0; l++) {
          pp.push_back(l);
          ww.push_back(w[l]);
          mm.push_back(z[l]);
          while (*(mm.end() - 1) >= *(mm.end() - 2)) {
            pp[pp.size() - 2] = l;
            mm[mm.size() - 2] = *(ww.end() - 1) * (*(mm.end() - 1)) +
              *(ww.end() - 2) * (*(mm.end() - 2));
              ww[ww.size() - 2] = *(ww.end() - 1) + *(ww.end() - 2);
              pp.pop_back();
              ww.pop_back();
              mm.pop_back();
              mm[mm.size() - 1] = *(mm.end() - 1) / *(ww.end() - 1);
          }
        }
      }

      if (copy) {
        pp.insert(pp.end(), rem_pp.begin(), rem_pp.end());
        ww.insert(ww.end(), rem_ww.begin(), rem_ww.end());
        mm.insert(mm.end(), rem_mm.begin(), rem_mm.end());
      }
      rem_pp.clear();
      rem_ww.clear();
      rem_mm.clear();
    }

    for (int l = 1; l < pp.size(); l++) {
      a0 = pp[l - 1] + 1;
      b0 = pp[l];
      for (int r = a0; r <= b0; r++) {
        cdf(k + 1, r) = mm[l];
      }
    }
  }
  vector<bool> update_pos;
  vector<bool> contained;
  update_pos.reserve(n_out);
  pos_new_x.reserve(n_out);
  contained.reserve(n_out);
  NumericVector w_tmp;
  List pos_x_tmp;

  for (int i = 0; i < n_out; i++) {
    pos_dbl = std::lower_bound(x.begin(), x.end(), x_out[i]);
    pos_new_x.push_back(std::distance(x.begin(), pos_dbl));
    if (*(pos_new_x.end() - 1) == x.size()) {
      w_lwr = 1.0;
      w_upr = 0.0;
      lwr = x.size() - 1;
      upr = x.size() - 1;
      contained.push_back(false);
      update_pos.push_back(false);
    } else if (*pos_dbl == x_out[i]) {
      lwr = *(pos_new_x.end() - 1);
      upr = *(pos_new_x.end() - 1);
      w_lwr = 1.0;
      w_upr = 0.0;
      contained.push_back(true);
      update_pos.push_back(false);
    } else if (*(pos_new_x.end() - 1) == 0) {
      upr = 0;
      lwr = 0;
      w_upr = 1.0;
      w_lwr = 0.0;
      contained.push_back(false);
      update_pos.push_back(true);
    } else {
      upr = *(pos_new_x.end() - 1);
      lwr = upr - 1;
      w_lwr = (x[upr] - x_out[i]) / (x[upr] - x[lwr]);
      w_upr = 1.0 - w_lwr;
      contained.push_back(false);
      update_pos.push_back(true);
    }

    for (int l = 0; l < n_thr + 2; l++) {
      cdf_lwr(l, i) = cdf(l, lwr);
      cdf_upr(l, i) = cdf(l, upr);
      cdf_oos(l, i) = w_lwr * cdf_lwr(l, i) + w_upr * cdf_upr(l, i);
    }
  }

  for (int i = 0; i < n_out; i++) {
    pos_x_tmp = clone(pos_x);
    w_tmp = clone(w);
    if (update_pos[i]) {
      for (int l = 0; l < pos_x_tmp.size(); l++) {
        psx = pos_x_tmp[l];
        for (int r = 0; r < psx.size(); r++) {
          if (pos_new_x[i] <= psx[r]) {
            psx[r] = psx[r] + 1;
          }
        }
        pos_x_tmp[l] = psx;
      }
    }

    if (contained[i]) {
      w_tmp[pos_new_x[i]] = w_tmp[pos_new_x[i]] + w_out;
    } else {
      w_tmp.insert(w_tmp.begin() + pos_new_x[i], w_out);
    }

    n = w_tmp.size();
    pp.erase(pp.begin() + 1, pp.end());
    ww.erase(ww.begin() + 1, ww.end());
    mm.erase(mm.begin() + 1, mm.end());
    z.fill(0.0);

    z[pos_new_x[i]] = w_out / w_tmp[pos_new_x[i]];
    pp.push_back(pos_new_x[i]);
    ww.push_back(sum(w_tmp[Range(0, pos_new_x[i])]));
    mm.push_back(w_out / *(ww.end() - 1));
    for (int l = 0; l < pos_new_x[i] + 1; l++) {
      cdf(0, l) = *(mm.end() - 1);
    }
    if (pos_new_x[i] < n - 1) {
      pp.push_back(n - 1);
      ww.push_back(sum(w_tmp[Range(pos_new_x[i] + 1, n - 1)]));
      mm.push_back(0.0);
    }

    for (int k = 0; k < n_psx; k++) {
      psx = pos_x_tmp[k];
      W_tmp = W[k];
      n_px = psx.length();
      for (int j = 0; j < n_px; j++) {
        copy = false;
        j0 = psx[j];

        z[j0] = z[j0] + W_tmp[j] / w_tmp[j0];
        pos_int = std::lower_bound(pp.begin(), pp.end(), j0);
        s0 = std::distance(pp.begin(), pos_int);
        a0 = pp[s0 - 1] + 1;
        b0 = pp[s0];

        std::vector<int> rem_pp;
        std::vector<double> rem_ww;
        std::vector<double> rem_mm;

        if (s0 < pp.size() - 1) {
          copy = true;
          for (int l = s0 + 1; l < pp.size(); l++) {
            rem_pp.push_back(pp[l]);
            rem_ww.push_back(ww[l]);
            rem_mm.push_back(mm[l]);
          }
        }

        pp.erase(pp.begin() + s0, pp.end());
        ww.erase(ww.begin() + s0, ww.end());
        mm.erase(mm.begin() + s0, mm.end());
        pp.push_back(j0);
        ww.push_back(sum(w_tmp[Range(a0, j0)]));
        mm.push_back(sum(w_tmp[Range(a0, j0)] * z[Range(a0, j0)]) / *(ww.end() - 1));

        while (*(mm.end() - 1) >= *(mm.end() - 2)) {
          pp[pp.size() - 2] = *(pp.end() - 1);
          mm[mm.size() - 2] = *(ww.end() - 1) * (*(mm.end() - 1)) +
            *(ww.end() - 2) * (*(mm.end() - 2));
            ww[ww.size() - 2] = *(ww.end() - 1) + *(ww.end() - 2);
            pp.pop_back();
            ww.pop_back();
            mm.pop_back();
            mm[mm.size() - 1] = *(mm.end() - 1) / *(ww.end() - 1);
        }

        if (j0 < b0) {
          for (int l = j0 + 1; l <= b0; l++) {
            pp.push_back(l);
            ww.push_back(w_tmp[l]);
            mm.push_back(z[l]);
            while (*(mm.end() - 1) >= *(mm.end() - 2)) {
              pp[pp.size() - 2] = l;
              mm[mm.size() - 2] = *(ww.end() - 1) * (*(mm.end() - 1)) +
                *(ww.end() - 2) * (*(mm.end() - 2));
                ww[ww.size() - 2] = *(ww.end() - 1) + *(ww.end() - 2);
                pp.pop_back();
                ww.pop_back();
                mm.pop_back();
                mm[mm.size() - 1] = *(mm.end() - 1) / *(ww.end() - 1);
            }
          }
        }

        if (copy) {
          pp.insert(pp.end(), rem_pp.begin(), rem_pp.end());
          ww.insert(ww.end(), rem_ww.begin(), rem_ww.end());
          mm.insert(mm.end(), rem_mm.begin(), rem_mm.end());
        }
        rem_pp.clear();
        rem_ww.clear();
        rem_mm.clear();
      }

      for (int l = 1; l < pp.size(); l++) {
        a0 = pp[l - 1] + 1;
        b0 = pp[l];
        for (int r = a0; r <= b0; r++) {
          cdf(k + 1, r) = mm[l];
        }
      }
    }

    for (int l = 0; l < n_thr + 2; l++) {
      cdf_lcnf(l, i) = cdf(l, pos_new_x[i]);
    }

    for (int l = 0; l < pos_new_x[i] + 1; l++) {
      cdf(0, l) = 0.0;
    }

    pp.erase(pp.begin() + 1, pp.end());
    ww.erase(ww.begin() + 1, ww.end());
    mm.erase(mm.begin() + 1, mm.end());
    pp.push_back(n - 1);
    ww.push_back(sum(w));
    mm.push_back(0.0);
    z.fill(0.0);

    for (int k = 0; k < n_psx; k++) {
      psx = pos_x_tmp[k];
      W_tmp = W[k];
      n_px = psx.length();
      for (int j = 0; j < n_px; j++) {
        copy = false;
        j0 = psx[j];

        z[j0] = z[j0] + W_tmp[j] / w_tmp[j0];
        pos_int = std::lower_bound(pp.begin(), pp.end(), j0);
        s0 = std::distance(pp.begin(), pos_int);
        a0 = pp[s0 - 1] + 1;
        b0 = pp[s0];

        std::vector<int> rem_pp;
        std::vector<double> rem_ww;
        std::vector<double> rem_mm;

        if (s0 < pp.size() - 1) {
          copy = true;
          for (int l = s0 + 1; l < pp.size(); l++) {
            rem_pp.push_back(pp[l]);
            rem_ww.push_back(ww[l]);
            rem_mm.push_back(mm[l]);
          }
        }

        pp.erase(pp.begin() + s0, pp.end());
        ww.erase(ww.begin() + s0, ww.end());
        mm.erase(mm.begin() + s0, mm.end());
        pp.push_back(j0);
        ww.push_back(sum(w_tmp[Range(a0, j0)]));
        mm.push_back(sum(w_tmp[Range(a0, j0)] * z[Range(a0, j0)]) / *(ww.end() - 1));

        while (*(mm.end() - 1) >= *(mm.end() - 2)) {
          pp[pp.size() - 2] = *(pp.end() - 1);
          mm[mm.size() - 2] = *(ww.end() - 1) * (*(mm.end() - 1)) +
            *(ww.end() - 2) * (*(mm.end() - 2));
            ww[ww.size() - 2] = *(ww.end() - 1) + *(ww.end() - 2);
            pp.pop_back();
            ww.pop_back();
            mm.pop_back();
            mm[mm.size() - 1] = *(mm.end() - 1) / *(ww.end() - 1);
        }

        if (j0 < b0) {
          for (int l = j0 + 1; l <= b0; l++) {
            pp.push_back(l);
            ww.push_back(w_tmp[l]);
            mm.push_back(z[l]);
            while (*(mm.end() - 1) >= *(mm.end() - 2)) {
              pp[pp.size() - 2] = l;
              mm[mm.size() - 2] = *(ww.end() - 1) * (*(mm.end() - 1)) +
                *(ww.end() - 2) * (*(mm.end() - 2));
                ww[ww.size() - 2] = *(ww.end() - 1) + *(ww.end() - 2);
                pp.pop_back();
                ww.pop_back();
                mm.pop_back();
                mm[mm.size() - 1] = *(mm.end() - 1) / *(ww.end() - 1);
            }
          }
        }

        if (copy) {
          pp.insert(pp.end(), rem_pp.begin(), rem_pp.end());
          ww.insert(ww.end(), rem_ww.begin(), rem_ww.end());
          mm.insert(mm.end(), rem_mm.begin(), rem_mm.end());
        }
        rem_pp.clear();
        rem_ww.clear();
        rem_mm.clear();
      }

      for (int l = 1; l < pp.size(); l++) {
        a0 = pp[l - 1] + 1;
        b0 = pp[l];
        for (int r = a0; r <= b0; r++) {
          cdf(k + 1, r) = mm[l];
        }
      }
    }

    for (int l = 0; l < n_thr + 2; l++) {
      cdf_ucnf(l, i) = cdf(l, pos_new_x[i]);
    }
    Rcpp::checkUserInterrupt();
  }

  return List::create(
    _["points"] = points,
    _["cdf_lwr"] = cdf_lwr,
    _["cdf_upr"] = cdf_upr,
    _["cdf_oos"] = cdf_oos,
    _["cdf_lcnf"] = cdf_lcnf,
    _["cdf_ucnf"] = cdf_ucnf
  );
}
