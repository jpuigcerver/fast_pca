#!/bin/bash
set -e;

SDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd );
DATA="${SDIR}/../../examples/gauss2d/data.ascii.mat";
PCA_REF="${SDIR}/pca.reference.mat";
DATA_PROJ_REF="${SDIR}/data.proj.reference.mat";
DATA_PROJ_NORM_REF="${SDIR}/data.proj.norm.reference.mat";
FAST_PCA_CMD="$1";

## Compute PCA
"${FAST_PCA_CMD}" -C -f ascii -p 2 "${DATA}" > pca.ascii.sp.mat;
"${FAST_PCA_CMD}" -C -d -f ascii -p 2 "${DATA}" > pca.ascii.dp.mat;
## Project data
"${FAST_PCA_CMD}" -P -f ascii -p 2 "${DATA}" -m pca.ascii.sp.mat \
    > proj.ascii.sp.mat;
"${FAST_PCA_CMD}" -P -d -f ascii -p 2 "${DATA}" -m pca.ascii.dp.mat \
    > proj.ascii.dp.mat;
## Project normalized data
"${FAST_PCA_CMD}" -P -f ascii -n -p 2 "${DATA}" -m pca.ascii.sp.mat \
    > proj.ascii.norm.sp.mat;
"${FAST_PCA_CMD}" -P -d -f ascii -n -p 2 "${DATA}" -m pca.ascii.dp.mat \
    > proj.ascii.norm.dp.mat;
## Compute PCA & Project data in a single pass
"${FAST_PCA_CMD}" -C -P -d -f ascii -p 2 -m pca.ascii.sp.2.mat "${DATA}" \
    > proj.ascii.sp.2.mat;
"${FAST_PCA_CMD}" -C -P -d -f ascii -p 2 -m pca.ascii.dp.2.mat "${DATA}" \
    > proj.ascii.dp.2.mat;
## Compute PCA & Project normalized data in a single pass
"${FAST_PCA_CMD}" -C -P -d -n -f ascii -p 2 -m pca.ascii.sp.2.mat "${DATA}" \
    > proj.ascii.norm.sp.2.mat;
"${FAST_PCA_CMD}" -C -P -d -n -f ascii -p 2 -m pca.ascii.dp.2.mat "${DATA}" \
    > proj.ascii.norm.dp.2.mat;

## Check PCA
"${SDIR}/../check_pca.sh" "${PCA_REF}" pca.ascii.sp.mat 1E-5;
"${SDIR}/../check_pca.sh" "${PCA_REF}" pca.ascii.dp.mat 1E-4;
## Check PCA 2
"${SDIR}/../check_pca.sh" "${PCA_REF}" pca.ascii.sp.2.mat 1E-5;
"${SDIR}/../check_pca.sh" "${PCA_REF}" pca.ascii.dp.2.mat 1E-4;
## Check data projections
"${SDIR}/../check_proj_ascii.sh" "${DATA_PROJ_REF}" proj.ascii.sp.mat 1E-2;
"${SDIR}/../check_proj_ascii.sh" "${DATA_PROJ_REF}" proj.ascii.dp.mat 1E-4;
## Check data projections 2
"${SDIR}/../check_proj_ascii.sh" "${DATA_PROJ_REF}" proj.ascii.sp.2.mat 1E-2;
"${SDIR}/../check_proj_ascii.sh" "${DATA_PROJ_REF}" proj.ascii.dp.2.mat 1E-4;
## Check normalized data projections
"${SDIR}/../check_proj_ascii.sh" "${DATA_PROJ_NORM_REF}" \
    proj.ascii.norm.sp.mat 1E-2;
"${SDIR}/../check_proj_ascii.sh" "${DATA_PROJ_NORM_REF}" \
    proj.ascii.norm.dp.mat 1E-4;
## Check normalized data projections 2
"${SDIR}/../check_proj_ascii.sh" "${DATA_PROJ_NORM_REF}" \
    proj.ascii.norm.sp.2.mat 1E-2;
"${SDIR}/../check_proj_ascii.sh" "${DATA_PROJ_NORM_REF}" \
    proj.ascii.norm.dp.2.mat 1E-4;

exit 0;
