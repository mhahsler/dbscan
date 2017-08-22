#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP dbscan_all_children(SEXP, SEXP, SEXP);
extern SEXP dbscan_buildDendrogram(SEXP);
extern SEXP dbscan_combine(SEXP, SEXP);
extern SEXP dbscan_computeStability(SEXP, SEXP, SEXP);
extern SEXP dbscan_computeVirtualNode(SEXP, SEXP);
extern SEXP dbscan_concat_int(SEXP);
extern SEXP dbscan_coreFromDist(SEXP, SEXP, SEXP);
extern SEXP dbscan_dbscan_density_int(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dbscan_dbscan_int(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dbscan_dendrogram_to_reach(SEXP);
extern SEXP dbscan_distToAdjacency(SEXP, SEXP);
extern SEXP dbscan_extractSemiSupervised(SEXP, SEXP, SEXP);
extern SEXP dbscan_extractUnsupervised(SEXP);
extern SEXP dbscan_fosc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dbscan_frNN_int(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dbscan_hclustMergeOrder(SEXP, SEXP);
extern SEXP dbscan_JP_int(SEXP, SEXP);
extern SEXP dbscan_kNN_int(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dbscan_lowerTri(SEXP);
extern SEXP dbscan_mrd(SEXP, SEXP);
extern SEXP dbscan_mrd_m(SEXP, SEXP);
extern SEXP dbscan_mst_to_dendrogram(SEXP);
extern SEXP dbscan_node_xy(SEXP, SEXP, SEXP);
extern SEXP dbscan_optics_int(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP dbscan_order_(SEXP);
extern SEXP dbscan_prims(SEXP, SEXP);
extern SEXP dbscan_reach_to_dendrogram(SEXP, SEXP);
extern SEXP dbscan_simplifiedTree(SEXP);
extern SEXP dbscan_SNN_sim_int(SEXP);
extern SEXP dbscan_validateConstraintList(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"dbscan_all_children",           (DL_FUNC) &dbscan_all_children,            3},
    {"dbscan_buildDendrogram",        (DL_FUNC) &dbscan_buildDendrogram,         1},
    {"dbscan_combine",                (DL_FUNC) &dbscan_combine,                 2},
    {"dbscan_computeStability",       (DL_FUNC) &dbscan_computeStability,        3},
    {"dbscan_computeVirtualNode",     (DL_FUNC) &dbscan_computeVirtualNode,      2},
    {"dbscan_concat_int",             (DL_FUNC) &dbscan_concat_int,              1},
    {"dbscan_coreFromDist",           (DL_FUNC) &dbscan_coreFromDist,            3},
    {"dbscan_dbscan_density_int",     (DL_FUNC) &dbscan_dbscan_density_int,      6},
    {"dbscan_dbscan_int",             (DL_FUNC) &dbscan_dbscan_int,             10},
    {"dbscan_dendrogram_to_reach",    (DL_FUNC) &dbscan_dendrogram_to_reach,     1},
    {"dbscan_distToAdjacency",        (DL_FUNC) &dbscan_distToAdjacency,         2},
    {"dbscan_extractSemiSupervised",  (DL_FUNC) &dbscan_extractSemiSupervised,   3},
    {"dbscan_extractUnsupervised",    (DL_FUNC) &dbscan_extractUnsupervised,     1},
    {"dbscan_fosc",                   (DL_FUNC) &dbscan_fosc,                    7},
    {"dbscan_frNN_int",               (DL_FUNC) &dbscan_frNN_int,                6},
    {"dbscan_hclustMergeOrder",       (DL_FUNC) &dbscan_hclustMergeOrder,        2},
    {"dbscan_JP_int",                 (DL_FUNC) &dbscan_JP_int,                  2},
    {"dbscan_kNN_int",                (DL_FUNC) &dbscan_kNN_int,                 6},
    {"dbscan_lowerTri",               (DL_FUNC) &dbscan_lowerTri,                1},
    {"dbscan_mrd",                    (DL_FUNC) &dbscan_mrd,                     2},
    {"dbscan_mrd_m",                  (DL_FUNC) &dbscan_mrd_m,                   2},
    {"dbscan_mst_to_dendrogram",      (DL_FUNC) &dbscan_mst_to_dendrogram,       1},
    {"dbscan_node_xy",                (DL_FUNC) &dbscan_node_xy,                 3},
    {"dbscan_optics_int",             (DL_FUNC) &dbscan_optics_int,              8},
    {"dbscan_order_",                 (DL_FUNC) &dbscan_order_,                  1},
    {"dbscan_prims",                  (DL_FUNC) &dbscan_prims,                   2},
    {"dbscan_reach_to_dendrogram",    (DL_FUNC) &dbscan_reach_to_dendrogram,     2},
    {"dbscan_simplifiedTree",         (DL_FUNC) &dbscan_simplifiedTree,          1},
    {"dbscan_SNN_sim_int",            (DL_FUNC) &dbscan_SNN_sim_int,             1},
    {"dbscan_validateConstraintList", (DL_FUNC) &dbscan_validateConstraintList,  2},
    {NULL, NULL, 0}
};

void R_init_dbscan(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
