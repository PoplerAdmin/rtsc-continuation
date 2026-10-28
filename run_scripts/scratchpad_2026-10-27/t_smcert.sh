source ~/rtsc_work/.venv/bin/activate
mkdir -p ~/rtsc_work/tetra_Al/smear_ref && cd ~/rtsc_work/tetra_Al/smear_ref
for i in 2 4 10; do cp ~/rtsc_work/bench/Al/run/Al-fcc_0GPa-d2e010742b/dfpt/a2F.dos$i .; done
for i in 2 4 10; do for c in "300 0.15" "1.2 0.10" "1.2 0.15"; do set -- $c
 echo "== a2F.dos$i T*=$1 mu_hi=$2"; rtsc certify a2F.dos$i --format qe_a2f_dos --target $1 --mu-hi $2 --mode bracket --material-id Al-smear-dos$i --out v_dos${i}_T$1_mu$2.json 2>&1 | grep -a "^Al-smear\|rror" | cut -c1-200; done; done
