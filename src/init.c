#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _dbscan_all_children(SEXP, SEXP, SEXP);
extern SEXP _dbscan_buildDendrogram(SEXP);
extern SEXP _dbscan_combine(SEXP, SEXP);
extern SEXP _dbscan_computeStability(SEXP, SEXP, SEXP);
extern SEXP _dbscan_computeVirtualNode(SEXP, SEXP);
extern SEXP _dbscan_concat_int(SEXP);
extern SEXP _dbscan_coreFromDist(SEXP, SEXP, SEXP);
extern SEXP _dbscan_dbscan_density_int(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dbscan_dbscan_int(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dbscan_dendrogram_to_reach(SEXP);
extern SEXP _dbscan_distToAdjacency(SEXP, SEXP);
extern SEXP _dbscan_extractSemiSupervised(SEXP, SEXP, SEXP, SEXP);
extern SEXP _dbscan_extractUnsupervised(SEXP, SEXP);
extern SEXP _dbscan_fosc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dbscan_frNN_int(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dbscan_hclustMergeOrder(SEXP, SEXP);
extern SEXP _dbscan_JP_int(SEXP, SEXP);
extern SEXP _dbscan_kNN_int(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dbscan_lowerTri(SEXP);
extern SEXP _dbscan_mrd(SEXP, SEXP);
extern SEXP _dbscan_mrd_m(SEXP, SEXP);
extern SEXP _dbscan_mst_to_dendrogram(SEXP);
extern SEXP _dbscan_node_xy(SEXP, SEXP, SEXP);
extern SEXP _dbscan_optics_int(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _dbscan_order_(SEXP);
extern SEXP _dbscan_prims(SEXP, SEXP);
extern SEXP _dbscan_reach_to_dendrogram(SEXP, SEXP);
extern SEXP _dbscan_simplifiedTree(SEXP);
extern SEXP _dbscan_SNN_sim_int(SEXP);
extern SEXP _dbscan_validateConstraintList(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_dbscan_all_children",           (DL_FUNC) &_dbscan_all_children,            3},
  {"_dbscan_buildDendrogram",        (DL_FUNC) &_dbscan_buildDendrogram,         1},
  {"_dbscan_combine",                (DL_FUNC) &_dbscan_combine,                 2},
  {"_dbscan_computeStability",       (DL_FUNC) &_dbscan_computeStability,        3},
  {"_dbscan_computeVirtualNode",     (DL_FUNC) &_dbscan_computeVirtualNode,      2},
  {"_dbscan_concat_int",             (DL_FUNC) &_dbscan_concat_int,              1},
  {"_dbscan_coreFromDist",           (DL_FUNC) &_dbscan_coreFromDist,            3},
  {"_dbscan_dbscan_density_int",     (DL_FUNC) &_dbscan_dbscan_density_int,      6},
  {"_dbscan_dbscan_int",             (DL_FUNC) &_dbscan_dbscan_int,             10},
  {"_dbscan_dendrogram_to_reach",    (DL_FUNC) &_dbscan_dendrogram_to_reach,     1},
  {"_dbscan_distToAdjacency",        (DL_FUNC) &_dbscan_distToAdjacency,         2},
  {"_dbscan_extractSemiSupervised",  (DL_FUNC) &_dbscan_extractSemiSupervised,   4},
  {"_dbscan_extractUnsupervised",    (DL_FUNC) &_dbscan_extractUnsupervised,     2},
  {"_dbscan_fosc",                   (DL_FUNC) &_dbscan_fosc,                    9},
  {"_dbscan_frNN_int",               (DL_FUNC) &_dbscan_frNN_int,                6},
  {"_dbscan_hclustMergeOrder",       (DL_FUNC) &_dbscan_hclustMergeOrder,        2},
  {"_dbscan_JP_int",                 (DL_FUNC) &_dbscan_JP_int,                  2},
  {"_dbscan_kNN_int",                (DL_FUNC) &_dbscan_kNN_int,                 6},
  {"_dbscan_lowerTri",               (DL_FUNC) &_dbscan_lowerTri,                1},
  {"_dbscan_mrd",                    (DL_FUNC) &_dbscan_mrd,                     2},
  {"_dbscan_mrd_m",                  (DL_FUNC) &_dbscan_mrd_m,                   2},
  {"_dbscan_mst_to_dendrogram",      (DL_FUNC) &_dbscan_mst_to_dendrogram,       1},
  {"_dbscan_node_xy",                (DL_FUNC) &_dbscan_node_xy,                 3},
  {"_dbscan_optics_int",             (DL_FUNC) &_dbscan_optics_int,              8},
  {"_dbscan_order_",                 (DL_FUNC) &_dbscan_order_,                  1},
  {"_dbscan_prims",                  (DL_FUNC) &_dbscan_prims,                   2},
  {"_dbscan_reach_to_dendrogram",    (DL_FUNC) &_dbscan_reach_to_dendrogram,     2},
  {"_dbscan_simplifiedTree",         (DL_FUNC) &_dbscan_simplifiedTree,          1},
  {"_dbscan_SNN_sim_int",            (DL_FUNC) &_dbscan_SNN_sim_int,             1},
  {"_dbscan_validateConstraintList", (DL_FUNC) &_dbscan_validateConstraintList,  2},
  {NULL, NULL, 0}
};

void R_init_dbscan(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}