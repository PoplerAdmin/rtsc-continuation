mkdir -p ~/rtsc_work/tetra_Al/src && cd ~/rtsc_work/tetra_Al/src
for f in PHonon/PH/elphon.f90 PHonon/PH/alpha2f.f90 PHonon/PH/phq_readin.f90 PHonon/Doc/INPUT_PH.def PHonon/PH/phq_setup.f90 PHonon/PH/do_phonon.f90; do
 curl -sSfL -o $(basename $f) https://raw.githubusercontent.com/QEF/q-e/qe-7.6/$f || echo FAIL $f; done
ls -la
grep -n "^#\|^##" ~/rtsc_work/rtsc_bundle/HANDOFF.md | head -40
ls ~/rtsc_work/tetra_Al/../ ; ls /mnt/g/ClaudeWorkspace/RTSC_Project/reports/ | head -50
