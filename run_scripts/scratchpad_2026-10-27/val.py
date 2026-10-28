import tomllib, os
from rtsc.dfpt.workflow import EPhSettings
from rtsc.dfpt.launcher import Launcher
for m in ['Pb','Nb','MgB2','H3S']:
    c=tomllib.load(open(os.path.expanduser(f'~/rtsc_work/bench/{m}/config.toml'),'rb'))
    s=EPhSettings.from_config(c['dfpt']['settings']); L=Launcher.from_config(c['dfpt']['launcher'])
    ps=c['dfpt']['pseudos']; ok=all(os.path.isfile(os.path.join(ps['dir'],f)) for f in ps['files'].values())
    print(m, 'settings OK', s.ecutwfc_Ry, s.k_scf, s.k_fit, s.q, '| launcher', L.nprocs, L.npool, '| pseudos exist', ok)
