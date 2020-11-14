#!/bin/csh

set minF = $1
set maxF = $2
set imgfn = $3.pdf
echo $imgfn

sac << EOF
r *sac
qdp off
bp bu co $minF $maxF n 1 p 2
sss
prs
quitsub
save $imgfn
q
EOF
