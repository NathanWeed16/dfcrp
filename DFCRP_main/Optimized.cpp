// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <vector>
#include <algorithm>  // for std::shuffle
#include <numeric>    // for std::iota
#include <random>     // for std::random_device, std::mt19937
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector relabel_cpp(IntegerVector cluster_vec, IntegerVector permutation) {
  int n = cluster_vec.size();
  IntegerVector relabeled = clone(cluster_vec);
  std::map<int, int> label_map;
  int next_label = 1;
  for (int i = 0; i < permutation.size(); i++) {
    int j = permutation[i] -1 ; // R is 1-indexed
    int old_label = cluster_vec[j];
    // If this old label hasnt been seen, assign a new one
    if (label_map.find(old_label) == label_map.end()) {
      label_map[old_label] = next_label;
      next_label++;
    }
  }

  // Relabel all entries based on the map built from permutation
  for (int i = 0; i < n; i++) {
    int old_label = cluster_vec[i];
    if (label_map.find(old_label) != label_map.end()) {
      relabeled[i] = label_map[old_label];
    }
  }
  return relabeled;
}

// assume relabel_cpp exists and is linked
// IntegerVector relabel_cpp(const IntegerVector &cluster_vec, const IntegerVector &perm);

// [[Rcpp::export]]
double log_dfcrp_pmf_cond_cpp(const IntegerVector &family_vec,
                                           const IntegerVector &cluster_vec,
                                           double alpha = 1.0,
                                           Nullable<IntegerVector> permutation = R_NilValue,
                                           int start_val = 1) {
  int n = cluster_vec.size();
  if (n == 0 || start_val > n) return 0.0;
  if (alpha <= 0.0) Rcpp::stop("alpha must be > 0");
  
  // --- build perm_input (1-based) for relabel_cpp ---
  IntegerVector perm_input(n);
  if (permutation.isNull()) {
    for (int i = 0; i < n; ++i) perm_input[i] = i + 1;
  } else {
    IntegerVector p = permutation.get();
    if ((int)p.size() != n) Rcpp::stop("permutation length mismatch");
    for (int i = 0; i < n; ++i) perm_input[i] = p[i];
  }
  
  IntegerVector cl = relabel_cpp(cluster_vec, perm_input); // 1-based cluster labels
  
  // --- maxima and capacities ---
  int max_cluster_label = 0;
  for (int i = 0; i < cl.size(); ++i) if (cl[i] > max_cluster_label) max_cluster_label = cl[i];
  int max_family = 0;
  for (int i = 0; i < family_vec.size(); ++i) if (family_vec[i] > max_family) max_family = family_vec[i];
  
  int Kcap = max_cluster_label + n + 3; // safe headroom
  int Fcap = max_family + 1; // families 1..max_family
  
  // --- core state ---
  std::vector<int> sizes(Kcap, 0);                // cluster sizes (1..)
  std::vector<double> log_sizes(Kcap, -INFINITY); // cache for log(size)
  std::vector< std::vector<unsigned char> > family_in_cluster(Fcap, std::vector<unsigned char>(Kcap, 1));
  // denom_per_family[f] = alpha + sum_{k allowed for f} sizes[k]
  std::vector<double> denom_per_family(Fcap, alpha);
  
  // For convenience: per-cluster list of families that were initially allowed
  // (We will iterate this list on updates; it contains all families 1..max_family)
  std::vector< std::vector<int> > families_by_cluster(Kcap);
  for (int k = 1; k < Kcap; ++k) {
    families_by_cluster[k].reserve(Fcap);
    for (int f = 1; f < Fcap; ++f) families_by_cluster[k].push_back(f);
  }
  
  // --- prefill for items before start_val (1-based start_val) ---
  for (int ii = 0; ii < start_val - 1; ++ii) {
    int j = perm_input[ii] - 1;   // 0-based index into family_vec and cl
    int xj = family_vec[j];       // family (1-based)
    int k_actual = cl[j];         // cluster label (1-based)
    
    if (k_actual < 1) Rcpp::stop("cluster labels must be >= 1");
    
    // update sizes/log_sizes
    sizes[k_actual] += 1;
    log_sizes[k_actual] = std::log((double)sizes[k_actual]);
    
    // mark family used in that cluster
    if (family_in_cluster[xj][k_actual]) {
      // currently family was allowed; remove it and update denom_per_family[xj]
      // denom_per_family[xj] currently contains alpha + sum of sizes (before removal),
      // and the cluster's size has already been increased above; but we must maintain
      // invariants consistently. To keep it simple, we will compute denom_per_family
      // after the prefill loop below. (Keeps prefill simple and correct.)
      family_in_cluster[xj][k_actual] = 0;
    }
  } // end prefill
  
  // --- initialize denom_per_family after prefill (safe and correct) ---
  for (int f = 1; f < Fcap; ++f) {
    double s = alpha;
    for (int k = 1; k <= max_cluster_label + n; ++k) {
      if (family_in_cluster[f][k] && sizes[k] > 0) s += (double)sizes[k];
    }
    denom_per_family[f] = s;
  }
  
  // After prefill, ensure next_empty and current_max_cluster computed
  int current_max_cluster = max_cluster_label;
  // if some cluster labels > max_cluster_label were created in cl (shouldn't), adjust
  for (int i = 0; i < cl.size(); ++i) if (cl[i] > current_max_cluster) current_max_cluster = cl[i];
  
  int next_empty = 1;
  while (next_empty <= current_max_cluster && sizes[next_empty] > 0) ++next_empty;
  
  double log_joint = 0.0;
  
  // --- main loop for items start_val..n ---
  for (int ii = start_val - 1; ii < n; ++ii) {
    int j = perm_input[ii] - 1;
    int xj = family_vec[j];
    int k_actual = cl[j];
    if (k_actual < 1) Rcpp::stop("cluster labels must be >= 1");
    
    if (k_actual > current_max_cluster) current_max_cluster = k_actual;
    while (next_empty <= current_max_cluster && sizes[next_empty] > 0) ++next_empty;
    int first_empty = (next_empty <= current_max_cluster ? next_empty : current_max_cluster + 1);
    
    // numerator: alpha if assigning to new cluster, else sizes[k_actual]
    double numer = (k_actual == first_empty ? alpha : (double)sizes[k_actual]);
    
    // Use cached denom for family xj (this is O(1))
    double denom = denom_per_family[xj];
    // guard
    if (denom <= 0.0) denom = std::numeric_limits<double>::min();
    if (numer <= 0.0) return -INFINITY;
    
    log_joint += std::log(numer) - std::log(denom);
    
    // --- Now update state (must use "before" family_in_cluster to determine updates) ---
    // We need to:
    // 1) For every family f that is still allowed in cluster k_actual (family_in_cluster[f][k_actual] == 1),
    //    add +1 to denom_per_family[f] (because sizes[k_actual] will increase by 1).
    // 2) If xj was previously allowed in k_actual, then after assignment family_in_cluster[xj][k_actual]=0,
    //    so denom_per_family[xj] must be reduced by the OLD sizes[k_actual] (since that cluster no longer
    //    contributes sizes[k_actual] to xj's denom); then we should NOT add the +1 for xj (because it's removed).
    //
    // To implement correctly, we:
    // - iterate families_by_cluster[k_actual] (which contains 1..max_family) and check family_in_cluster[f][k_actual].
    //   This costs O(F) per assignment (good when F < K).
    // - apply updates accordingly.
    
    int old_size = sizes[k_actual];
    
    // iterate families known for this cluster (we use families_by_cluster to avoid constructing the loop over 1..F each time)
    const std::vector<int> &famlist = families_by_cluster[k_actual];
    for (int idxf = 0; idxf < (int)famlist.size(); ++idxf) {
      int f = famlist[idxf];
      if (!family_in_cluster[f][k_actual]) continue; // already disallowed -> no change
      if (f == xj) {
        // For the current family: removing cluster entirely from that family's allowed set
        // subtract the old_size (if old_size>0) from denom_per_family[xj]
        if (old_size > 0) denom_per_family[f] -= (double)old_size;
        // mark disallowed (cluster will not be part of xj's denom anymore)
        family_in_cluster[f][k_actual] = 0;
      } else {
        // For other families still allowed in this cluster: they see sizes[k_actual] increase by 1
        denom_per_family[f] += 1.0;
      }
    }
    
    // --- Finally increment the cluster size and update caches for sizes/log_sizes ---
    sizes[k_actual] = old_size + 1;
    log_sizes[k_actual] = std::log((double)sizes[k_actual]);
    
    // If we used the first_empty slot, advance next_empty
    if (k_actual == first_empty) {
      ++next_empty;
      while (next_empty <= current_max_cluster && sizes[next_empty] > 0) ++next_empty;
    }
  } // end main loop
  
  return log_joint;
}


