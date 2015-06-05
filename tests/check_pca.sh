#!/bin/bash
set -e;

[ $# -ne 3 ] && {
    echo "Usage: ${0##*/} pca_ref.mat pca_test.mat tolerance" >&2;
    exit 1;
}

octave --eval "
function check_equal(A, B, tol, msg)
  sA = size(A);
  sB = size(B);
  if sum(sA ~= sB) ~= 0
    fprintf(stderr, '%s. Sizes do not match (%d,%d) vs (%d,%d)', ...
            msg, sA(1), sA(2), sB(1), sB(2));
    exit(1);
  else
    s_a = abs(A) + abs(B);
    d_a = abs(A - B);
    s_a(s_a < tol) = 1;
    max_err = max(max(d_a ./ s_a));
    if max_err > tol
      fprintf(stderr, '%s. Maximum Relative Error: %g', msg, max_err);
      exit(1);
    endif
  end
endfunction

load '$1';
Eref=E; Rref=R; Mref=M; Sref=S; Dref=D; Vref=V;
load '$2';

check_equal(Eref, E, $3, 'E does not match the reference');
check_equal(Rref, R, $3, 'R does not match the reference');
check_equal(Mref, M, $3, 'M does not match the reference');
check_equal(Sref, S, $3, 'S does not match the reference');
check_equal(Dref, D, $3, 'D does not match the reference');
check_equal(Vref, V, $3, 'V does not match the reference');
" || { echo "File \"$2\" does not match the reference \"$1\"!" >&2; exit 1; }

exit 0;
