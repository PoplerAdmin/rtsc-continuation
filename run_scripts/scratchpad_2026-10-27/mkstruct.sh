B=~/rtsc_work/bench
mkdir -p $B/{Pb,Nb,MgB2,H3S}/structures
cat > $B/Pb/structures/Pb_fcc.vasp <<P
Pb fcc primitive cell, a=4.95 A (experimental; vc-relaxed by the DFT stage)
 1.0
     0.000000000   2.475000000   2.475000000
     2.475000000   0.000000000   2.475000000
     2.475000000   2.475000000   0.000000000
 Pb
   1
Direct
  0.0 0.0 0.0
P
cat > $B/Pb/structures/manifest.toml <<P
# Benchmark input. Only facts about the input: target pressure.
["Pb_fcc"]
material_id = "Pb-fcc@0GPa"
pressure_GPa = 0.0
P
cat > $B/Nb/structures/Nb_bcc.vasp <<P
Nb bcc primitive cell, a=3.30 A (experimental; vc-relaxed by the DFT stage)
 1.0
    -1.650000000   1.650000000   1.650000000
     1.650000000  -1.650000000   1.650000000
     1.650000000   1.650000000  -1.650000000
 Nb
   1
Direct
  0.0 0.0 0.0
P
cat > $B/Nb/structures/manifest.toml <<P
# Benchmark input. Only facts about the input: target pressure.
["Nb_bcc"]
material_id = "Nb-bcc@0GPa"
pressure_GPa = 0.0
P
cat > $B/MgB2/structures/MgB2_P6mmm.vasp <<P
MgB2 P6/mmm (AlB2 type), a=3.086 c=3.524 A (experimental; vc-relaxed by the DFT stage)
 1.0
     3.086000000   0.000000000   0.000000000
    -1.543000000   2.672555332   0.000000000
     0.000000000   0.000000000   3.524000000
 Mg B
   1 2
Direct
  0.000000000 0.000000000 0.000000000
  0.333333333 0.666666667 0.500000000
  0.666666667 0.333333333 0.500000000
P
cat > $B/MgB2/structures/manifest.toml <<P
# Benchmark input. Only facts about the input: target pressure.
["MgB2_P6mmm"]
material_id = "MgB2-P6mmm@0GPa"
pressure_GPa = 0.0
P
cat > $B/H3S/structures/H3S_Im-3m.vasp <<P
H3S Im-3m PRIMITIVE bcc cell (4 atoms) from prototypes/H3S_Im-3m.vasp (a=3.09 A conventional); vc-relaxed at 200 GPa by the DFT stage
 1.0
    -1.545000000   1.545000000   1.545000000
     1.545000000  -1.545000000   1.545000000
     1.545000000   1.545000000  -1.545000000
 S H
   1 3
Cartesian
  0.000000000 0.000000000 0.000000000
  1.545000000 0.000000000 0.000000000
  0.000000000 1.545000000 0.000000000
  0.000000000 0.000000000 1.545000000
P
cat > $B/H3S/structures/manifest.toml <<P
# Benchmark input. Only facts about the input: target pressure.
# Harmonic PBE Im-3m H3S is dynamically stable only above ~180 GPa (Errea et al., Nature 532, 81 (2016));
# 200 GPa keeps the harmonic DFPT inside its stability range. Expt Tc ~200 K at ~155 GPa (Drozdov 2015).
["H3S_Im-3m"]
material_id = "H3S-Im-3m@200GPa"
pressure_GPa = 200.0
P
find $B -path $B/Al -prune -o -type f -print
cat ~/rtsc_work/conv/results.txt
