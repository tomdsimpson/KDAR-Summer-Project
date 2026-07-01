import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pathlib

PDG_MAP = {
     14: r'$\nu_\mu$',
    -14: r'$\bar{\nu}_\mu$',
     12: r'$\nu_e$',
    -12: r'$\bar{\nu}_e$',
}

P_PDGS = {
    "mu" : ([13,-13], "blue"),
    "pi" : ([211,-211], "orange"),
    "K0L" : ( [130], "red"),
    "Kpm" : ([321, -321], "green")
    }

# --- Load Data --- #

def load_mcnu_data(filepath, energy_max, n_chunks = 35):
    print('Loading data...')

    hdr_chunks  = []
    mcnu_chunks = []
    hdr_cols = ['pot']

    with pd.HDFStore(filepath, 'r') as store:

        print(store["opflash_1"].columns)

        # Browse evt structure
        for i in range(n_chunks):
            print(f'  Chunk {i}...', flush=True)
            hdr_chunks.append(store[f'hdr_{i}'][hdr_cols].copy())
            mcnu_chunks.append(store[f'mcnu_{i}'].copy())

    hdr  = pd.concat(hdr_chunks)
    mcnu = pd.concat(mcnu_chunks)

    pot_mc_total = hdr['pot'][hdr['pot'] > 0].sum()

    print(f'Non-zero POT entries : {(hdr["pot"] > 0).sum():,}')
    print(f'Total MC POT         : {pot_mc_total:.3e}')
    print(f'Total events         : {len(mcnu):,}')

    mcnu['weight'] = 1.0 / pot_mc_total  # nu / POT / GeV

    # --- Energy Cut & Binning --- #
    mcnu_cut = mcnu[mcnu['E'] < energy_max].copy()
    return mcnu_cut



def load_evt_data(filepath, chunk_i = 0):
    
    with pd.HDFStore(filepath, 'r') as store:

        print(f'  Chunk {chunk_i}...', flush=True)
        evt = store[f"evt_{chunk_i}"].copy()

    return evt



# --- Plot --- #
if __name__ == "__main__":

    ENERGY_MAX = 1.0    # GeV
    N_BINS     = 100

    this_path = pathlib.Path(__file__).parent.resolve()
    path = this_path / '../../Desktop/mc_MCP2025B_prodgenie_corsika_proton_rockbox_lowenergydirt_sbnd_CV_caf_flat_caf_sbnd_general_dfs.df'

    mcnu_cut = load_mcnu_data(filepath=path, energy_max=ENERGY_MAX)
    bins     = np.linspace(0, ENERGY_MAX, N_BINS + 1)
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))

    # Left: flux by flavour
    ax = axs[0]
    colors = ['blue', 'red', 'green', 'orange']

    for (pdg_code, label), color in zip(PDG_MAP.items(), colors):
        mask = mcnu_cut['pdg'] == pdg_code
        ax.hist(mcnu_cut.loc[mask, 'E'], bins=bins, weights=mcnu_cut.loc[mask, 'weight'], color=color, histtype="step", label=label)

    ax.set_xlabel('True Neutrino Energy (GeV)')
    ax.set_ylabel(r'$\nu$ / POT / GeV / cm$^2$')
    ax.set_yscale('log')
    ax.set_title('Neutrino Flux by Flavour')
    ax.legend()
    ax.set_xlim(0, ENERGY_MAX)

    # Right: nu_mu CC vs NC
    ax = axs[1]
    numu_mask = mcnu_cut['pdg'] == 14

    for label, (p_pdgs, color) in P_PDGS.items():
        mask = numu_mask & (np.isin(mcnu_cut['parent_pdg'], p_pdgs))
        ax.hist(mcnu_cut.loc[mask, 'E'], bins=bins, weights=mcnu_cut.loc[mask, 'weight'], color=color, label=label, histtype="step")

    ax.set_xlabel('True Neutrino Energy (GeV)')
    ax.set_ylabel(r'$\nu$ / POT / GeV / cm$^2$')
    ax.set_yscale('log')
    ax.set_title(r'$\nu_\mu$ Flux by Parent')
    ax.legend()
    ax.set_xlim(0, ENERGY_MAX)

    plt.tight_layout()
    plt.savefig('neutrino_flux.png', dpi=150)
    plt.show()
    print('Done — saved neutrino_flux.png')
