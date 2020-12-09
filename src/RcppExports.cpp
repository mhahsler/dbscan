// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// JP_int
IntegerVector JP_int(IntegerMatrix nn, unsigned int kt);
RcppExport SEXP _dbscan_JP_int(SEXP nnSEXP, SEXP ktSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< unsigned int >::type kt(ktSEXP);
    rcpp_result_gen = Rcpp::wrap(JP_int(nn, kt));
    return rcpp_result_gen;
END_RCPP
}
// SNN_sim_int
IntegerMatrix SNN_sim_int(IntegerMatrix nn, LogicalVector jp);
RcppExport SEXP _dbscan_SNN_sim_int(SEXP nnSEXP, SEXP jpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type nn(nnSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type jp(jpSEXP);
    rcpp_result_gen = Rcpp::wrap(SNN_sim_int(nn, jp));
    return rcpp_result_gen;
END_RCPP
}
// ANN_cleanup
void ANN_cleanup();
RcppExport SEXP _dbscan_ANN_cleanup() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    ANN_cleanup();
    return R_NilValue;
END_RCPP
}
// dbscan_int
IntegerVector dbscan_int(NumericMatrix data, double eps, int minPts, NumericVector weights, int borderPoints, int type, int bucketSize, int splitRule, double approx, List frNN);
RcppExport SEXP _dbscan_dbscan_int(SEXP dataSEXP, SEXP epsSEXP, SEXP minPtsSEXP, SEXP weightsSEXP, SEXP borderPointsSEXP, SEXP typeSEXP, SEXP bucketSizeSEXP, SEXP splitRuleSEXP, SEXP approxSEXP, SEXP frNNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type minPts(minPtsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< int >::type borderPoints(borderPointsSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type splitRule(splitRuleSEXP);
    Rcpp::traits::input_parameter< double >::type approx(approxSEXP);
    Rcpp::traits::input_parameter< List >::type frNN(frNNSEXP);
    rcpp_result_gen = Rcpp::wrap(dbscan_int(data, eps, minPts, weights, borderPoints, type, bucketSize, splitRule, approx, frNN));
    return rcpp_result_gen;
END_RCPP
}
// dbscan_density_int
IntegerVector dbscan_density_int(NumericMatrix data, double eps, int type, int bucketSize, int splitRule, double approx);
RcppExport SEXP _dbscan_dbscan_density_int(SEXP dataSEXP, SEXP epsSEXP, SEXP typeSEXP, SEXP bucketSizeSEXP, SEXP splitRuleSEXP, SEXP approxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type splitRule(splitRuleSEXP);
    Rcpp::traits::input_parameter< double >::type approx(approxSEXP);
    rcpp_result_gen = Rcpp::wrap(dbscan_density_int(data, eps, type, bucketSize, splitRule, approx));
    return rcpp_result_gen;
END_RCPP
}
// frNN_int
List frNN_int(NumericMatrix data, double eps, int type, int bucketSize, int splitRule, double approx);
RcppExport SEXP _dbscan_frNN_int(SEXP dataSEXP, SEXP epsSEXP, SEXP typeSEXP, SEXP bucketSizeSEXP, SEXP splitRuleSEXP, SEXP approxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type splitRule(splitRuleSEXP);
    Rcpp::traits::input_parameter< double >::type approx(approxSEXP);
    rcpp_result_gen = Rcpp::wrap(frNN_int(data, eps, type, bucketSize, splitRule, approx));
    return rcpp_result_gen;
END_RCPP
}
// frNN_query_int
List frNN_query_int(NumericMatrix data, NumericMatrix query, double eps, int type, int bucketSize, int splitRule, double approx);
RcppExport SEXP _dbscan_frNN_query_int(SEXP dataSEXP, SEXP querySEXP, SEXP epsSEXP, SEXP typeSEXP, SEXP bucketSizeSEXP, SEXP splitRuleSEXP, SEXP approxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type query(querySEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type splitRule(splitRuleSEXP);
    Rcpp::traits::input_parameter< double >::type approx(approxSEXP);
    rcpp_result_gen = Rcpp::wrap(frNN_query_int(data, query, eps, type, bucketSize, splitRule, approx));
    return rcpp_result_gen;
END_RCPP
}
// kNN_int
List kNN_int(NumericMatrix data, int k, int type, int bucketSize, int splitRule, double approx);
RcppExport SEXP _dbscan_kNN_int(SEXP dataSEXP, SEXP kSEXP, SEXP typeSEXP, SEXP bucketSizeSEXP, SEXP splitRuleSEXP, SEXP approxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type splitRule(splitRuleSEXP);
    Rcpp::traits::input_parameter< double >::type approx(approxSEXP);
    rcpp_result_gen = Rcpp::wrap(kNN_int(data, k, type, bucketSize, splitRule, approx));
    return rcpp_result_gen;
END_RCPP
}
// kNN_query_int
List kNN_query_int(NumericMatrix data, NumericMatrix query, int k, int type, int bucketSize, int splitRule, double approx);
RcppExport SEXP _dbscan_kNN_query_int(SEXP dataSEXP, SEXP querySEXP, SEXP kSEXP, SEXP typeSEXP, SEXP bucketSizeSEXP, SEXP splitRuleSEXP, SEXP approxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type query(querySEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type splitRule(splitRuleSEXP);
    Rcpp::traits::input_parameter< double >::type approx(approxSEXP);
    rcpp_result_gen = Rcpp::wrap(kNN_query_int(data, query, k, type, bucketSize, splitRule, approx));
    return rcpp_result_gen;
END_RCPP
}
// optics_int
List optics_int(NumericMatrix data, double eps, int minPts, int type, int bucketSize, int splitRule, double approx, List frNN);
RcppExport SEXP _dbscan_optics_int(SEXP dataSEXP, SEXP epsSEXP, SEXP minPtsSEXP, SEXP typeSEXP, SEXP bucketSizeSEXP, SEXP splitRuleSEXP, SEXP approxSEXP, SEXP frNNSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type data(dataSEXP);
    Rcpp::traits::input_parameter< double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< int >::type minPts(minPtsSEXP);
    Rcpp::traits::input_parameter< int >::type type(typeSEXP);
    Rcpp::traits::input_parameter< int >::type bucketSize(bucketSizeSEXP);
    Rcpp::traits::input_parameter< int >::type splitRule(splitRuleSEXP);
    Rcpp::traits::input_parameter< double >::type approx(approxSEXP);
    Rcpp::traits::input_parameter< List >::type frNN(frNNSEXP);
    rcpp_result_gen = Rcpp::wrap(optics_int(data, eps, minPts, type, bucketSize, splitRule, approx, frNN));
    return rcpp_result_gen;
END_RCPP
}
// distToAdjacency
List distToAdjacency(IntegerVector constraints, const int N);
RcppExport SEXP _dbscan_distToAdjacency(SEXP constraintsSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type constraints(constraintsSEXP);
    Rcpp::traits::input_parameter< const int >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(distToAdjacency(constraints, N));
    return rcpp_result_gen;
END_RCPP
}
// buildDendrogram
List buildDendrogram(List hcl);
RcppExport SEXP _dbscan_buildDendrogram(SEXP hclSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type hcl(hclSEXP);
    rcpp_result_gen = Rcpp::wrap(buildDendrogram(hcl));
    return rcpp_result_gen;
END_RCPP
}
// all_children
IntegerVector all_children(List hier, int key, bool leaves_only);
RcppExport SEXP _dbscan_all_children(SEXP hierSEXP, SEXP keySEXP, SEXP leaves_onlySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type hier(hierSEXP);
    Rcpp::traits::input_parameter< int >::type key(keySEXP);
    Rcpp::traits::input_parameter< bool >::type leaves_only(leaves_onlySEXP);
    rcpp_result_gen = Rcpp::wrap(all_children(hier, key, leaves_only));
    return rcpp_result_gen;
END_RCPP
}
// node_xy
NumericMatrix node_xy(List cl_tree, List cl_hierarchy, int cid);
RcppExport SEXP _dbscan_node_xy(SEXP cl_treeSEXP, SEXP cl_hierarchySEXP, SEXP cidSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type cl_tree(cl_treeSEXP);
    Rcpp::traits::input_parameter< List >::type cl_hierarchy(cl_hierarchySEXP);
    Rcpp::traits::input_parameter< int >::type cid(cidSEXP);
    rcpp_result_gen = Rcpp::wrap(node_xy(cl_tree, cl_hierarchy, cid));
    return rcpp_result_gen;
END_RCPP
}
// simplifiedTree
List simplifiedTree(List cl_tree);
RcppExport SEXP _dbscan_simplifiedTree(SEXP cl_treeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type cl_tree(cl_treeSEXP);
    rcpp_result_gen = Rcpp::wrap(simplifiedTree(cl_tree));
    return rcpp_result_gen;
END_RCPP
}
// computeStability
List computeStability(const List hcl, const int minPts, bool compute_glosh);
RcppExport SEXP _dbscan_computeStability(SEXP hclSEXP, SEXP minPtsSEXP, SEXP compute_gloshSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List >::type hcl(hclSEXP);
    Rcpp::traits::input_parameter< const int >::type minPts(minPtsSEXP);
    Rcpp::traits::input_parameter< bool >::type compute_glosh(compute_gloshSEXP);
    rcpp_result_gen = Rcpp::wrap(computeStability(hcl, minPts, compute_glosh));
    return rcpp_result_gen;
END_RCPP
}
// validateConstraintList
List validateConstraintList(List& constraints, int n);
RcppExport SEXP _dbscan_validateConstraintList(SEXP constraintsSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List& >::type constraints(constraintsSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(validateConstraintList(constraints, n));
    return rcpp_result_gen;
END_RCPP
}
// computeVirtualNode
double computeVirtualNode(IntegerVector noise, List constraints);
RcppExport SEXP _dbscan_computeVirtualNode(SEXP noiseSEXP, SEXP constraintsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type noise(noiseSEXP);
    Rcpp::traits::input_parameter< List >::type constraints(constraintsSEXP);
    rcpp_result_gen = Rcpp::wrap(computeVirtualNode(noise, constraints));
    return rcpp_result_gen;
END_RCPP
}
// fosc
NumericVector fosc(List cl_tree, std::string cid, std::list<int>& sc, List cl_hierarchy, bool prune_unstable_leaves, const double alpha, bool useVirtual, const int n_constraints, List constraints);
RcppExport SEXP _dbscan_fosc(SEXP cl_treeSEXP, SEXP cidSEXP, SEXP scSEXP, SEXP cl_hierarchySEXP, SEXP prune_unstable_leavesSEXP, SEXP alphaSEXP, SEXP useVirtualSEXP, SEXP n_constraintsSEXP, SEXP constraintsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type cl_tree(cl_treeSEXP);
    Rcpp::traits::input_parameter< std::string >::type cid(cidSEXP);
    Rcpp::traits::input_parameter< std::list<int>& >::type sc(scSEXP);
    Rcpp::traits::input_parameter< List >::type cl_hierarchy(cl_hierarchySEXP);
    Rcpp::traits::input_parameter< bool >::type prune_unstable_leaves(prune_unstable_leavesSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type useVirtual(useVirtualSEXP);
    Rcpp::traits::input_parameter< const int >::type n_constraints(n_constraintsSEXP);
    Rcpp::traits::input_parameter< List >::type constraints(constraintsSEXP);
    rcpp_result_gen = Rcpp::wrap(fosc(cl_tree, cid, sc, cl_hierarchy, prune_unstable_leaves, alpha, useVirtual, n_constraints, constraints));
    return rcpp_result_gen;
END_RCPP
}
// extractUnsupervised
List extractUnsupervised(List cl_tree, bool prune_unstable);
RcppExport SEXP _dbscan_extractUnsupervised(SEXP cl_treeSEXP, SEXP prune_unstableSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type cl_tree(cl_treeSEXP);
    Rcpp::traits::input_parameter< bool >::type prune_unstable(prune_unstableSEXP);
    rcpp_result_gen = Rcpp::wrap(extractUnsupervised(cl_tree, prune_unstable));
    return rcpp_result_gen;
END_RCPP
}
// extractSemiSupervised
List extractSemiSupervised(List cl_tree, List constraints, float alpha, bool prune_unstable_leaves);
RcppExport SEXP _dbscan_extractSemiSupervised(SEXP cl_treeSEXP, SEXP constraintsSEXP, SEXP alphaSEXP, SEXP prune_unstable_leavesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type cl_tree(cl_treeSEXP);
    Rcpp::traits::input_parameter< List >::type constraints(constraintsSEXP);
    Rcpp::traits::input_parameter< float >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type prune_unstable_leaves(prune_unstable_leavesSEXP);
    rcpp_result_gen = Rcpp::wrap(extractSemiSupervised(cl_tree, constraints, alpha, prune_unstable_leaves));
    return rcpp_result_gen;
END_RCPP
}
// reach_to_dendrogram
List reach_to_dendrogram(const Rcpp::List reachability, const NumericVector pl_order);
RcppExport SEXP _dbscan_reach_to_dendrogram(SEXP reachabilitySEXP, SEXP pl_orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type reachability(reachabilitySEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type pl_order(pl_orderSEXP);
    rcpp_result_gen = Rcpp::wrap(reach_to_dendrogram(reachability, pl_order));
    return rcpp_result_gen;
END_RCPP
}
// dendrogram_to_reach
List dendrogram_to_reach(const Rcpp::List x);
RcppExport SEXP _dbscan_dendrogram_to_reach(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(dendrogram_to_reach(x));
    return rcpp_result_gen;
END_RCPP
}
// mst_to_dendrogram
List mst_to_dendrogram(const NumericMatrix mst);
RcppExport SEXP _dbscan_mst_to_dendrogram(SEXP mstSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix >::type mst(mstSEXP);
    rcpp_result_gen = Rcpp::wrap(mst_to_dendrogram(mst));
    return rcpp_result_gen;
END_RCPP
}
// mrd
NumericVector mrd(NumericVector dm, NumericVector cd);
RcppExport SEXP _dbscan_mrd(SEXP dmSEXP, SEXP cdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type dm(dmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cd(cdSEXP);
    rcpp_result_gen = Rcpp::wrap(mrd(dm, cd));
    return rcpp_result_gen;
END_RCPP
}
// mrd_m
NumericMatrix mrd_m(NumericMatrix dm, NumericVector cd);
RcppExport SEXP _dbscan_mrd_m(SEXP dmSEXP, SEXP cdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type dm(dmSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type cd(cdSEXP);
    rcpp_result_gen = Rcpp::wrap(mrd_m(dm, cd));
    return rcpp_result_gen;
END_RCPP
}
// coreFromDist
NumericVector coreFromDist(const NumericVector dist, const int n, const int minPts);
RcppExport SEXP _dbscan_coreFromDist(SEXP distSEXP, SEXP nSEXP, SEXP minPtsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type dist(distSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    Rcpp::traits::input_parameter< const int >::type minPts(minPtsSEXP);
    rcpp_result_gen = Rcpp::wrap(coreFromDist(dist, n, minPts));
    return rcpp_result_gen;
END_RCPP
}
// prims
NumericMatrix prims(const NumericVector x_dist, const int n);
RcppExport SEXP _dbscan_prims(SEXP x_distSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x_dist(x_distSEXP);
    Rcpp::traits::input_parameter< const int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(prims(x_dist, n));
    return rcpp_result_gen;
END_RCPP
}
// order_
IntegerVector order_(NumericVector x);
RcppExport SEXP _dbscan_order_(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(order_(x));
    return rcpp_result_gen;
END_RCPP
}
// hclustMergeOrder
List hclustMergeOrder(NumericMatrix mst, IntegerVector o);
RcppExport SEXP _dbscan_hclustMergeOrder(SEXP mstSEXP, SEXP oSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type mst(mstSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type o(oSEXP);
    rcpp_result_gen = Rcpp::wrap(hclustMergeOrder(mst, o));
    return rcpp_result_gen;
END_RCPP
}
// lowerTri
IntegerVector lowerTri(IntegerMatrix m);
RcppExport SEXP _dbscan_lowerTri(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(lowerTri(m));
    return rcpp_result_gen;
END_RCPP
}
// combine
NumericVector combine(const NumericVector& t1, const NumericVector& t2);
RcppExport SEXP _dbscan_combine(SEXP t1SEXP, SEXP t2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type t1(t1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type t2(t2SEXP);
    rcpp_result_gen = Rcpp::wrap(combine(t1, t2));
    return rcpp_result_gen;
END_RCPP
}
// concat_int
IntegerVector concat_int(List const& container);
RcppExport SEXP _dbscan_concat_int(SEXP containerSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List const& >::type container(containerSEXP);
    rcpp_result_gen = Rcpp::wrap(concat_int(container));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_dbscan_JP_int", (DL_FUNC) &_dbscan_JP_int, 2},
    {"_dbscan_SNN_sim_int", (DL_FUNC) &_dbscan_SNN_sim_int, 2},
    {"_dbscan_ANN_cleanup", (DL_FUNC) &_dbscan_ANN_cleanup, 0},
    {"_dbscan_dbscan_int", (DL_FUNC) &_dbscan_dbscan_int, 10},
    {"_dbscan_dbscan_density_int", (DL_FUNC) &_dbscan_dbscan_density_int, 6},
    {"_dbscan_frNN_int", (DL_FUNC) &_dbscan_frNN_int, 6},
    {"_dbscan_frNN_query_int", (DL_FUNC) &_dbscan_frNN_query_int, 7},
    {"_dbscan_kNN_int", (DL_FUNC) &_dbscan_kNN_int, 6},
    {"_dbscan_kNN_query_int", (DL_FUNC) &_dbscan_kNN_query_int, 7},
    {"_dbscan_optics_int", (DL_FUNC) &_dbscan_optics_int, 8},
    {"_dbscan_distToAdjacency", (DL_FUNC) &_dbscan_distToAdjacency, 2},
    {"_dbscan_buildDendrogram", (DL_FUNC) &_dbscan_buildDendrogram, 1},
    {"_dbscan_all_children", (DL_FUNC) &_dbscan_all_children, 3},
    {"_dbscan_node_xy", (DL_FUNC) &_dbscan_node_xy, 3},
    {"_dbscan_simplifiedTree", (DL_FUNC) &_dbscan_simplifiedTree, 1},
    {"_dbscan_computeStability", (DL_FUNC) &_dbscan_computeStability, 3},
    {"_dbscan_validateConstraintList", (DL_FUNC) &_dbscan_validateConstraintList, 2},
    {"_dbscan_computeVirtualNode", (DL_FUNC) &_dbscan_computeVirtualNode, 2},
    {"_dbscan_fosc", (DL_FUNC) &_dbscan_fosc, 9},
    {"_dbscan_extractUnsupervised", (DL_FUNC) &_dbscan_extractUnsupervised, 2},
    {"_dbscan_extractSemiSupervised", (DL_FUNC) &_dbscan_extractSemiSupervised, 4},
    {"_dbscan_reach_to_dendrogram", (DL_FUNC) &_dbscan_reach_to_dendrogram, 2},
    {"_dbscan_dendrogram_to_reach", (DL_FUNC) &_dbscan_dendrogram_to_reach, 1},
    {"_dbscan_mst_to_dendrogram", (DL_FUNC) &_dbscan_mst_to_dendrogram, 1},
    {"_dbscan_mrd", (DL_FUNC) &_dbscan_mrd, 2},
    {"_dbscan_mrd_m", (DL_FUNC) &_dbscan_mrd_m, 2},
    {"_dbscan_coreFromDist", (DL_FUNC) &_dbscan_coreFromDist, 3},
    {"_dbscan_prims", (DL_FUNC) &_dbscan_prims, 2},
    {"_dbscan_order_", (DL_FUNC) &_dbscan_order_, 1},
    {"_dbscan_hclustMergeOrder", (DL_FUNC) &_dbscan_hclustMergeOrder, 2},
    {"_dbscan_lowerTri", (DL_FUNC) &_dbscan_lowerTri, 1},
    {"_dbscan_combine", (DL_FUNC) &_dbscan_combine, 2},
    {"_dbscan_concat_int", (DL_FUNC) &_dbscan_concat_int, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_dbscan(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
