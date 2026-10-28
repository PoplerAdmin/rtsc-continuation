#!/bin/bash
# Reproducer: QE 7.6 ph.x SIGSEGV in adddvscf_ (from elphel) for electron_phonon + USPP/PAW.
# Usage: QE_BIN=/path/to/qe/bin SSSP_DIR=/path/SSSP_1.3.0_PBE_precision SG15_DIR=/path/sg15 ./run.sh [case]
#   case: al_paw (default) | al_paw_simple | al_paw_noeph | al_sg15 | mg_us | all
# Pseudopotentials (not included; SSSP 1.3.0 PBE precision / SG15 1.2):
#   Al.pbe-n-kjpaw_psl.1.0.0.UPF  md5 cfc449ca30b5f3223ec38ddd88ac046d
#   mg_pbe_v1.4.uspp.F.UPF        md5 8ffbd8f729fef71095aac5bd8316fb1f
#   Al_ONCV_PBE-1.2.upf (SG15)    md5 2d9e7bc02555220d62fad08c3c4858e7
# Serial, ~15-25 s per case.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
: "${QE_BIN:?set QE_BIN}" "${SSSP_DIR:?set SSSP_DIR}"; SG15_DIR=${SG15_DIR:-}
export OMP_NUM_THREADS=1 GFORTRAN_ERROR_BACKTRACE=1
one() { # name system phinput
  local n=$1 sys=$2 ph=$3 d=${WORK:-$PWD/work}/$1
  rm -rf "$d"; mkdir -p "$d"; cd "$d"
  for k in fit scf; do sed -e "s#SSSP_DIR#$SSSP_DIR#" -e "s#SG15_DIR#$SG15_DIR#" "$HERE/$sys.$k.in" > $k.in; done
  cp "$HERE/$ph" ph.in
  [ "$ph" != ph_noeph.in ] && { "$QE_BIN/pw.x" -in fit.in > fit.out 2>&1 || { echo "$n: fit pw.x failed"; return; }; }
  "$QE_BIN/pw.x" -in scf.in > scf.out 2>&1 || { echo "$n: scf pw.x failed"; return; }
  "$QE_BIN/ph.x" -in ph.in > ph.out 2> ph.err; rc=$?
  if grep -q SIGSEGV ph.out ph.err; then r=SEGV; elif [ $rc = 0 ] && grep -q 'JOB DONE' ph.out; then r=OK; else r="FAIL rc=$rc"; fi
  echo "$n: $r  ($(grep -m1 'Program PHONON' ph.out | awk '{print $3}'))"
}
case ${1:-al_paw} in
  al_paw)        one al_paw al_paw ph_eph.in ;;
  al_paw_simple) one al_paw_simple al_paw ph_eph_simple.in ;;
  al_paw_noeph)  one al_paw_noeph al_paw ph_noeph.in ;;
  al_sg15)       one al_sg15 al_sg15 ph_eph.in ;;
  mg_us)         one mg_us mg_us ph_eph.in ;;
  all) for c in al_paw al_paw_simple al_paw_noeph al_sg15 mg_us; do (bash "$0" $c); done ;;
esac
