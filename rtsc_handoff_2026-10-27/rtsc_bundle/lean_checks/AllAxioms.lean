/-
  Axiom audit over EVERY declaration of the bundle's Lean libraries (not a hand-picked list).
  For each constant defined in a module whose name starts with TcR1, PipelineLogic, MultiBand or DiscoveryLogic,
  print one line `AXIOMS <decl> : [ax1, ax2, ...]` (Lean's own `collectAxioms`, the function
  behind `#print axioms`), then a summary line `AXIOM-AUDIT decls=<n> bad=<m>`.
  `verify_all.sh` accepts only propext, Classical.choice and Quot.sound.
  Run from the Lean project directory:  lake env lean <path>/AllAxioms.lean
-/
import TcR1
import PipelineLogic
import MultiBand
import DiscoveryLogic
open Lean Elab Command

#eval show CommandElabM Unit from do
  let env ← getEnv
  let roots : List Name := [`TcR1, `PipelineLogic, `MultiBand, `DiscoveryLogic]
  let allowed : List Name := [``propext, ``Classical.choice, ``Quot.sound]
  let mut n := 0
  let mut bad := 0
  let mut lines : Array String := #[]
  for (c, _) in env.constants.map₁.toList do
    let some idx := env.getModuleIdxFor? c | continue
    let some modName := env.header.moduleNames[idx.toNat]? | continue
    unless roots.any (fun r => r.isPrefixOf modName) do continue
    let axs ← collectAxioms c
    n := n + 1
    let isBad := axs.any (fun a => !(allowed.contains a))
    if isBad then bad := bad + 1
    lines := lines.push s!"AXIOMS {c} : {axs.toList}{if isBad then "   <-- NOT ALLOWED" else ""}"
  for l in lines.qsort (· < ·) do
    IO.println l
  IO.println s!"AXIOM-AUDIT decls={n} bad={bad}"
