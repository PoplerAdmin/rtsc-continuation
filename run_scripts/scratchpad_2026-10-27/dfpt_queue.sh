# Coarse DFPT for Pb, Nb, H3S, MgB2 (inputs for the EPW convergence stage), run one after another.
. ~/rtsc_work/.venv/bin/activate
ulimit -s unlimited
for m in Pb Nb H3S MgB2; do
  cd ~/rtsc_work/bench/$m
  sed -i 's/^el_ph_sigma_Ry = 0.005.*/el_ph_sigma_Ry = 0.01      # a-priori rule: smallest smearing 0.01 Ry (0.005 Ry under-resolved on Al 32^3)/' config.toml
  grep "^el_ph_sigma_Ry" config.toml
  rm -rf run
  echo "$(date) START $m" >> ~/rtsc_work/bench/queue.log
  ( time nice -n 10 rtsc run -c config.toml --no-resume ) > run.log 2>&1
  echo "$(date) END $m: $(grep -E 'ERROR|REJECTED|UNDECIDED|CERTIFIED|EXCLUDED' run.log | awk '$2>0' | tr '\n' ' ')" >> ~/rtsc_work/bench/queue.log
done
