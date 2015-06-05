#!/bin/bash
set -e;

SDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd );
DATA="${SDIR}/../../examples/gauss2d/data.matlab.mat";
DATA_PROJ_REF="${SDIR}/data.proj.reference.mat";
DATA_PROJ_NORM_REF="${SDIR}/data.proj.norm.reference.mat";
FAST_PCA_CMD="$1";

## Compute PCA & Project data in a single pass
"${FAST_PCA_CMD}" -C -P -f mat4 -m pca.matlab.sp.mat "${DATA}" \
    > proj.matlab.sp.mat;
"${FAST_PCA_CMD}" -C -P -d -f mat4 -m pca.matlab.dp.mat "${DATA}" \
    > proj.matlab.dp.mat;

## Check data projections
"${SDIR}/../check_proj_octave.sh" "${DATA_PROJ_REF}" proj.matlab.sp.mat 1E-2;
"${SDIR}/../check_proj_octave.sh" "${DATA_PROJ_REF}" proj.matlab.dp.mat 1E-8;

exit 0;
