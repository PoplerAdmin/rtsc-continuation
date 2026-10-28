# Resume after the 21:07 WSL crash: Nb resumes (ph.x recover via rtsc's own resume), then H3S and MgB2 fresh.
. ~/rtsc_work/.venv/bin/activate
ulimit -s unlimited
cd ~/rtsc_work/bench/Nb
echo "$(date) RESUME Nb (after WSL crash)" >> ~/rtsc_work/bench/queue.log
( time nice -n 10 rtsc run -c config.toml ) > run_resume.log 2>&1
echo "$(date) END Nb: $(grep -E 'ERROR|REJECTED|UNDECIDED|CERTIFIED|EXCLUDED' run_resume.log | awk '$2>0' | tr '\n' ' ')" >> ~/rtsc_work/bench/queue.log
for m in H3S MgB2; do
  cd ~/rtsc_work/bench/$m
  sed -i 's/^el_ph_sigma_Ry = 0.005.*/el_ph_sigma_Ry = 0.01      # a-priori rule: smallest smearing 0.01 Ry (0.005 Ry under-resolved on Al 32^3)/' config.toml
  rm -rf run
  echo "$(date) START $m" >> ~/rtsc_work/bench/queue.log
  ( time nice -n 10 rtsc run -c config.toml --no-resume ) > run.log 2>&1
  echo "$(date) END $m: $(grep -E 'ERROR|REJECTED|UNDECIDED|CERTIFIED|EXCLUDED' run.log | awk '$2>0' | tr '\n' ' ')" >> ~/rtsc_work/bench/queue.log
done
