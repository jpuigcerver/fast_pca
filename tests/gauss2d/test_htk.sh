#!/bin/bash
set -e;

SDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd );
DATA="${SDIR}/../../examples/gauss2d/data.htk.mat";
DATA_PROJ_REF="${SDIR}/data.proj.reference.mat";
DATA_PROJ_NORM_REF="${SDIR}/data.proj.norm.reference.mat";
FAST_PCA_CMD="$1";

## Compute PCA & Project data in a single pass
"${FAST_PCA_CMD}" -C -P -f htk -m pca.htk.sp.mat "${DATA}" \
    > proj.htk.sp.mat;
"${FAST_PCA_CMD}" -C -P -d -f htk -m pca.htk.dp.mat "${DATA}" \
    > proj.htk.dp.mat;

## Check data projections
"${SDIR}/../check_proj_htk.sh" "${DATA_PROJ_REF}" proj.htk.sp.mat 1E-2;
"${SDIR}/../check_proj_htk.sh" "${DATA_PROJ_REF}" proj.htk.dp.mat 1E-4;

exit 0;
