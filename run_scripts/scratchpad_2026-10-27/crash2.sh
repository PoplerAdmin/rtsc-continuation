Q=$HOME/micromamba/envs/qe/bin; export OMP_NUM_THREADS=1
cd ~/rtsc_work/crashchar
tail -25 al_paw_eph_simple/ph.out; echo; head -30 al_paw_eph_simple/ph.err
which gdb; file -L $Q/ph.x; nm $Q/ph.x 2>&1 | head -3; nm -D $Q/ph.x | grep -ci elph
ls $Q/../lib | grep -i -E 'libqe|phlib|libph' | head
# variants on Al PAW, gamma only
base=al_paw_eph_interp
for v in nofildvscf skipgamma twostep; do
  rm -rf v_$v; cp -r $base v_$v; cd v_$v; rm -rf tmp/_ph0 dyn* ph.out ph.err
  case $v in
   nofildvscf) sed -i "s/fildvscf='dvscf',//" ph.in ;;
   skipgamma) sed -i "s/nq1=1, nq2=1, nq3=1/nq1=2, nq2=2, nq3=2, start_q=2, last_q=2/" ph.in ;;
   twostep) sed "s/electron_phonon='interpolated', el_ph_sigma=0.005, el_ph_nsigma=10,//" ph.in > ph0.in
            nice -n 10 $Q/ph.x -in ph0.in > ph0.out 2>&1; echo "twostep phonon-only rc=$?"
            sed -i "s/trans=.true./trans=.false./" ph.in ;;
  esac
  nice -n 10 $Q/ph.x -in ph.in > ph.out 2> ph.err; echo "$v rc=$? segv=$(grep -c SIGSEGV ph.err) lambda=$(grep -c 'lambda(' ph.out)"
  tail -4 ph.out
  cd ..
done
