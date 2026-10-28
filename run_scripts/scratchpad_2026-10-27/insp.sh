source ~/rtsc_work/.venv/bin/activate; python /mnt/g/ClaudeWorkspace/RTSC_Project/run_scripts/scratchpad_2026-10-27/insp.py $1
grep -i "fermi\|wall\|WALL" ~/rtsc_work/epw_Al/g1/discarded_nofrozen/epw_ef_K*/epw.out | head; cat ~/rtsc_work/epw_Al/g1/pw_bands/DONE; grep -i "highest\|Fermi" ~/rtsc_work/epw_Al/g1/pw_bands/bands.out
