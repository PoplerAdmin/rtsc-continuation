cd ~/rtsc_work/tetra_Al/src; mkdir -p ref; cd ref
B=https://raw.githubusercontent.com/QEF/q-e/qe-7.6/PHonon/examples/tetra_example/reference
for f in al.a2F.out aluminum.a2F.dat al.dyn1.elph.1 al.elph.out; do curl -sSfL -o $f $B/$f; done
cat al.a2F.out | sed -n '1,200p' | grep -v "^\s*$" | tail -40
head -5 aluminum.a2F.dat; cat al.dyn1.elph.1
grep -n "k points\|nbnd\|Kohn-Sham states\|PHONON.*WALL\|lambda(" al.elph.out | head -30
