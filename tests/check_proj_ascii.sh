#!/bin/bash
set -e;

[ $# -ne 3 ] && {
    echo "Usage: ${0##*/} proj_ref.mat proj_test.mat tolerance" >&2;
    exit 1;
}

octave --eval "
function check_equal(A, B, tol, msg)
  sA = size(A);
  sB = size(B);
  if sum(sA ~= sB) ~= 0
    fprintf(stderr, '%s. Sizes do not match (%d,%d) vs (%d,%d)\n', ...
            msg, sA(1), sA(2), sB(1), sB(2));
    exit(1);
  else
    s_a = abs(A) + abs(B);
    d_a = abs(A - B);
    s_a(s_a < tol) = 1;
    max_err = max(max(d_a ./ s_a));
    if max_err > tol
      fprintf(stderr, '%s. Maximum Relative Error: %g\n', msg, max_err);
      exit(1);
    endif
  end
endfunction

load '$1';
Xref = X;
X = load('$2');

check_equal(Xref, X, $3, 'Projected data does not match the reference');
" || { echo "File \"$2\" does not match the reference \"$1\"!" >&2; exit 1; }

exit 0;
