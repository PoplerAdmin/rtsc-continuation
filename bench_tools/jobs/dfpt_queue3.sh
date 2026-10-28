# Resume after the 2026-10-27 23:21 WSL VM shutdown: Nb resumes (rtsc restarts ph.x with recover=.true.),
# then H3S and MgB2 fresh. Same as run_scripts/scratchpad_2026-10-27/dfpt_queue2.sh except that an old
# run/ directory is moved aside instead of removed (evidence is never deleted without a hash record).
. ~/rtsc_work/.venv/bin/activate
ulimit -s unlimited
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd ~/rtsc_work/bench/Nb
[ -f run_resume.log ] && mv run_resume.log run_resume_crash3_$(date +%Y%m%dT%H%M%S).log
echo "$(date) RESUME Nb (after 23:21 WSL VM shutdown)" >> ~/rtsc_work/bench/queue.log
( time nice -n 10 rtsc run -c config.toml ) > run_resume.log 2>&1
echo "$(date) END Nb: $(grep -E 'ERROR|REJECTED|UNDECIDED|CERTIFIED|EXCLUDED' run_resume.log | awk '$2>0' | tr '\n' ' ')" >> ~/rtsc_work/bench/queue.log
for m in H3S MgB2; do
  cd ~/rtsc_work/bench/$m
  sed -i 's/^el_ph_sigma_Ry = 0.005.*/el_ph_sigma_Ry = 0.01      # a-priori rule: smallest smearing 0.01 Ry (0.005 Ry under-resolved on Al 32^3)/' config.toml
  [ -d run ] && mv run run_prev_$(date +%Y%m%dT%H%M%S)
  echo "$(date) START $m" >> ~/rtsc_work/bench/queue.log
  ( time nice -n 10 rtsc run -c config.toml --no-resume ) > run.log 2>&1
  echo "$(date) END $m: $(grep -E 'ERROR|REJECTED|UNDECIDED|CERTIFIED|EXCLUDED' run.log | awk '$2>0' | tr '\n' ' ')" >> ~/rtsc_work/bench/queue.log
done
