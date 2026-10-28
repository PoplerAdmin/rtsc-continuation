source ~/rtsc_work/.venv/bin/activate; cd ~/rtsc_work/tetra_Al/src/ref
S=/mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27
rtsc certify aluminum.a2F.dat --units Ry --rule bins --target 300 --expected-lambda 0.3963915104345056 --out ref_v.json | head -2
python $S/dump_spectra.py aluminum.a2F.dat Ry bins error 0.10 ref_sp.json 0.3963915104345056
rtsc recheck ref_v.json --spectra ref_sp.json
