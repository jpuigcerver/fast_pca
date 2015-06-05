#!/bin/bash
set -e;

SDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd );
DATA="${SDIR}/../../examples/gauss2d/data.vbosch.mat";
DATA_PROJ_REF="${SDIR}/data.proj.reference.mat";
DATA_PROJ_NORM_REF="${SDIR}/data.proj.norm.reference.mat";
FAST_PCA_CMD="$1";

## Compute PCA & Project data in a single pass
"${FAST_PCA_CMD}" -C -P -f vbosch -m pca.vbosch.sp.mat "${DATA}" \
    > proj.vbosch.sp.mat;
"${FAST_PCA_CMD}" -C -P -d -f vbosch -m pca.vbosch.dp.mat "${DATA}" \
    > proj.vbosch.dp.mat;

tail -n+2 proj.vbosch.sp.mat > proj.vbosch.noheader.sp.mat;
tail -n+2 proj.vbosch.dp.mat > proj.vbosch.noheader.dp.mat;

## Check data projections
"${SDIR}/../check_proj_ascii.sh" "${DATA_PROJ_REF}" \
    proj.vbosch.noheader.sp.mat 1E-2;
"${SDIR}/../check_proj_ascii.sh" "${DATA_PROJ_REF}" \
    proj.vbosch.noheader.dp.mat 1E-8;

exit 0;
