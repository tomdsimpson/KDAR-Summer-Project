import numpy as np
import matplotlib.pyplot as plt
import uproot
import pathlib
import particle
from particle import Particle

PATH = pathlib.Path(__file__).parent.resolve()
DET_I = 5 #SBND ?
DET_A = 4*np.pi # Find detector area

flavours = {
     14: (r"$\nu_\mu$",     "blue"),
    -14: (r"$\bar\nu_\mu$", "orange"),
     12: (r"$\nu_e$",       "green"),
    -12: (r"$\bar\nu_e$",   "red"),
}



with uproot.open(PATH / "NuBeam_production_BooNE_50m_I174000A_0.dk2nu.root") as f:
    
    tree = f[f.keys()[0]]
    meta = f[f.keys()[2]]

    nuray_E   = tree["nuray"]["nuray.E"].array()
    nuray_wgt = tree["nuray"]["nuray.wgt"].array()
    nimpwt   = tree["decay"]["decay.nimpwt"].array()
    ntype     = tree["decay"]["decay.ntype"].array()
    ptype     = tree["decay"]["decay.ptype"].array() 

    pots = meta["pots"].array()

# Select SBND
energies = nuray_E[:,DET_I]
wgts     = nuray_wgt[:,DET_I]
tot_wgts = nimpwt * wgts / (DET_A * pots[0])

# Parents
is_numu      = (ntype == 14)
parent_groups = {
    "mu" : (np.isin(ptype, [13,-13]), "blue"),
    "pi" : (np.isin(ptype, [211,-211]), "orange"),
    "K0L" : (np.isin(ptype, [130]), "red"),
    "Kpm" : (np.isin(ptype, [321, -321]), "green")
}

fig, (ax1, ax2) = plt.subplots(1,2,figsize=(12, 5))
bins = np.linspace(0, 10, 101)  # GeV


for pdg, (label, colour) in flavours.items():
    mask = (ntype == pdg)
    ax1.hist(energies[mask], bins=bins, weights=tot_wgts[mask],
            label=label, color=colour, histtype="step", linewidth=1.5)

ax1.set_xlabel("Neutrino Energy (GeV)")
ax1.set_xlim(0,10)
ax1.set_ylabel(r"$\nu$ / 4$\pi$ m$^2$ / GeV / POT")
ax1.set_title("Neutrino Flux")
ax1.legend()
ax1.set_yscale("log")

for label, (parent_mask, colour) in parent_groups.items():
    mask = is_numu & parent_mask
    ax2.hist(energies[mask], bins=bins, weights=tot_wgts[mask],
            label=label, color=colour, histtype="step", linewidth=1.5)

ax2.set_xlabel("Neutrino Energy (GeV)")
ax2.set_ylabel(r"$\nu$ / cm$^2$ / GeV / POT")
ax2.set_title(r"$\nu_\mu$ Flux by Parent Particle")
ax2.legend()
ax2.set_yscale("log")






plt.tight_layout()
plt.show()




'''
mask = np.array([321 in x for x in ancestor_pdgs])
k_energies = energies[mask]
o_energies = energies[~mask]

print(f"n KDAR: {len(k_energies)}")
print(f"n Other: {len(o_energies)}")

k_energies = k_energies[(k_energies >= 0) & (k_energies <= 1)]
o_energies = o_energies[(o_energies >= 0) & (o_energies <= 1)]

plt.hist([k_energies, o_energies], bins=150, histtype="barstacked")
plt.show()
'''


