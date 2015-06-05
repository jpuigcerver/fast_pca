#!/bin/bash
set -e;

SDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd );
DATA="${SDIR}/../../examples/gauss2d/data.octave.mat";
PCA_REF="${SDIR}/pca.reference.mat";
DATA_PROJ_REF="${SDIR}/data.proj.reference.mat";
DATA_PROJ_NORM_REF="${SDIR}/data.proj.norm.reference.mat";
FAST_PCA_CMD="$1";

## Compute PCA
"${FAST_PCA_CMD}" -C -f octave "${DATA}" > pca.octave.sp.mat;
"${FAST_PCA_CMD}" -C -d -f octave "${DATA}" > pca.octave.dp.mat;
## Project data
"${FAST_PCA_CMD}" -P -f octave "${DATA}" -m pca.octave.sp.mat > proj.octave.sp.mat;
"${FAST_PCA_CMD}" -P -d -f octave "${DATA}" -m pca.octave.dp.mat > proj.octave.dp.mat;
## Project normalized data
"${FAST_PCA_CMD}" -P -f octave -n "${DATA}" -m pca.octave.sp.mat > proj.octave.norm.sp.mat;
"${FAST_PCA_CMD}" -P -d -f octave -n "${DATA}" -m pca.octave.dp.mat > proj.octave.norm.dp.mat;
## Compute PCA & Project data in a single pass
"${FAST_PCA_CMD}" -C -P -d -f octave -m pca.octave.sp.2.mat "${DATA}" > proj.octave.sp.2.mat;
"${FAST_PCA_CMD}" -C -P -d -f octave -m pca.octave.dp.2.mat "${DATA}" > proj.octave.dp.2.mat;
## Compute PCA & Project normalized data in a single pass
"${FAST_PCA_CMD}" -C -P -d -n -f octave -m pca.octave.sp.2.mat "${DATA}" > proj.octave.norm.sp.2.mat;
"${FAST_PCA_CMD}" -C -P -d -n -f octave -m pca.octave.dp.2.mat "${DATA}" > proj.octave.norm.dp.2.mat;


## Check PCA
"${SDIR}/../check_pca.sh" "${PCA_REF}" pca.octave.sp.mat 1E-5;
"${SDIR}/../check_pca.sh" "${PCA_REF}" pca.octave.dp.mat 1E-8;
## Check PCA 2
"${SDIR}/../check_pca.sh" "${PCA_REF}" pca.octave.sp.2.mat 1E-5;
"${SDIR}/../check_pca.sh" "${PCA_REF}" pca.octave.dp.2.mat 1E-8;
## Check data projections
"${SDIR}/../check_proj_octave.sh" "${DATA_PROJ_REF}" proj.octave.sp.mat 1E-2;
"${SDIR}/../check_proj_octave.sh" "${DATA_PROJ_REF}" proj.octave.dp.mat 1E-8;
## Check data projections 2
"${SDIR}/../check_proj_octave.sh" "${DATA_PROJ_REF}" proj.octave.sp.2.mat 1E-2;
"${SDIR}/../check_proj_octave.sh" "${DATA_PROJ_REF}" proj.octave.dp.2.mat 1E-8;
## Check normalized data projections
"${SDIR}/../check_proj_octave.sh" "${DATA_PROJ_NORM_REF}" proj.octave.norm.sp.mat 1E-2;
"${SDIR}/../check_proj_octave.sh" "${DATA_PROJ_NORM_REF}" proj.octave.norm.dp.mat 1E-8;
## Check normalized data projections 2
"${SDIR}/../check_proj_octave.sh" "${DATA_PROJ_NORM_REF}" proj.octave.norm.sp.2.mat 1E-2;
"${SDIR}/../check_proj_octave.sh" "${DATA_PROJ_NORM_REF}" proj.octave.norm.dp.2.mat 1E-8;


exit 0;
