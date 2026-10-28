# EPW stage, step 1-2 for Al: collect DFPT data (EPW layout) and nscf on the full coarse k grid.
set -e
Q=$HOME/micromamba/envs/qe/bin
SRC=$(ls -d ~/rtsc_work/bench/Al/run/Al-*)/dfpt
W=~/rtsc_work/epw_Al; mkdir -p $W/save $W/tmp; cd $W
# --- collect (port of epw_pp.py; dyn files are named dyn<i> here, dvscf in tmp/_ph0)
nq=$(grep -A1 NUMBER_OF_Q_POINTS $SRC/tmp/_ph0/pwscf.phsave/control_ph.xml | tail -1 | tr -d ' ')
echo "irreducible q: $nq"
for i in $(seq 1 $nq); do
  cp $SRC/dyn$i save/pwscf.dyn_q$i
  if [ $i -eq 1 ]; then s=$SRC/tmp/_ph0; else s=$SRC/tmp/_ph0/pwscf.q_$i; fi
  ln -sf $s/pwscf.dvscf1 save/pwscf.dvscf_q$i
done
rm -rf save/pwscf.phsave; cp -r $SRC/tmp/_ph0/pwscf.phsave save/
# --- scf charge density for nscf (copy save dir without wavefunctions)
rm -rf tmp/pwscf.save; mkdir -p tmp/pwscf.save
cp $SRC/tmp/pwscf.save/*.xml $SRC/tmp/pwscf.save/charge-density* tmp/pwscf.save/ 2>/dev/null || true
cp $SRC/tmp/pwscf.save/*.upf tmp/pwscf.save/ 2>/dev/null || true
ls tmp/pwscf.save
# --- nscf on the full 8x8x8 grid (commensurate with q 8^3), 20 bands (4 semicore + valence + empty)
NK=8
python3 - <<PY > nscf.in
src = open("$SRC/scf.in").read()
src = src.replace("calculation = 'scf'", "calculation = 'nscf'")
src = src.replace("degauss = 0.02", "degauss = 0.02\n    nbnd = 20\n    nosym = .true.\n    noinv = .true.")
src = src.split("K_POINTS")[0]
n = $NK; N = n**3
lines = ["K_POINTS crystal", str(N)]
for i in range(n):
    for j in range(n):
        for k in range(n):
            lines.append(f"  {i/n:.8f} {j/n:.8f} {k/n:.8f} {1/N:.6e}")
print(src + "\n".join(lines))
PY
grep -E "calculation|nbnd|nosym" nscf.in
( time $Q/mpirun -np 8 --bind-to none $Q/pw.x -nk 8 -in nscf.in > nscf.out ) 2>&1 | tail -3
grep -E "Fermi energy|JOB DONE" nscf.out
