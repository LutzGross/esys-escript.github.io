/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef ML_VIZ_VTK_H
#define ML_VIZ_VTK_H

#include "ml_include.h"
#include "ml_viz_stats.h"

#ifndef ML_CPP
#ifdef __cplusplus
extern "C" {
#endif
#endif

int ML_Aggregate_VisualizeVTK( ML_Aggregate_Viz_Stats info,
               char base_filename[], ML_Comm *comm, double * vector);

int ML_PlotVTK(int Npoints, double* x, double* y, double* z,
               char base_filename[], ML_Comm *comm, double * vector);

#ifndef ML_CPP
#ifdef __cplusplus
}
#endif
#endif

#endif /* #ifndef ML_VIZ_VTK_H */


#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

