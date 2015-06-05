#!/bin/bash
set -e;

SDIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd );
DATA_SP="${SDIR}/../../examples/gauss2d/data.binary.sp.mat";
DATA_DP="${SDIR}/../../examples/gauss2d/data.binary.dp.mat";
DATA_PROJ_REF="${SDIR}/data.proj.reference.mat";
DATA_PROJ_NORM_REF="${SDIR}/data.proj.norm.reference.mat";
FAST_PCA_CMD="$1";

## Compute PCA & Project data in a single pass
"${FAST_PCA_CMD}" -C -P -f binary -p 2 -m pca.binary.sp.mat "${DATA_SP}" > proj.binary.sp.mat;
"${FAST_PCA_CMD}" -C -P -d -f binary -p 2 -m pca.binary.dp.mat "${DATA_DP}" > proj.binary.dp.mat;

LC_NUMERIC=C od -An -t f4 -w8 proj.binary.sp.mat > proj.binary2ascii.sp.mat;
LC_NUMERIC=C od -An -t f8 -w16 proj.binary.dp.mat > proj.binary2ascii.dp.mat;

## Check data projections
"${SDIR}/../check_proj_ascii.sh" "${DATA_PROJ_REF}" proj.binary2ascii.sp.mat 1E-2;
"${SDIR}/../check_proj_ascii.sh" "${DATA_PROJ_REF}" proj.binary2ascii.dp.mat 1E-8;

exit 0;
