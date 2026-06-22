import numpy as np
import uproot
import particle
import matplotlib.pyplot as plt

DET_A = 4*np.pi*1e4 #Detector area

flavours = {
     14: (r"$\nu_\mu$",     "blue"),
    -14: (r"$\bar\nu_\mu$", "orange"),
     12: (r"$\nu_e$",       "green"),
    -12: (r"$\bar\nu_e$",   "red"),
}

class Neutrinos_dk2nu():

    def __init__(self, root_path, det_i=5):

        self.bins = 51 # Define granularity

        with uproot.open("./DATA/NuBeam_production_BooNE_50m_I174000A_0.dk2nu.root") as f:

            tree = f[f.keys()[0]]
            meta = f[f.keys()[2]]

            self.energies   = tree["nuray"]["nuray.E"].array()[:,det_i]
            self.ntype     = tree["decay"]["decay.ntype"].array()
            self.ptype     = tree["decay"]["decay.ptype"].array() 
            self.pot = meta["pots"].array()[0]

            nuray_wgt = tree["nuray"]["nuray.wgt"].array()[:,det_i]
            nimpwt    = tree["decay"]["decay.nimpwt"].array()
            self.tot_wgts = nimpwt * nuray_wgt / (DET_A * self.pot)

            self.px = tree["nuray"]["nuray.px"].array()[:,det_i]
            self.py = tree["nuray"]["nuray.py"].array()[:,det_i]
            self.pz = tree["nuray"]["nuray.pz"].array()[:,det_i]

            self.decay_x = tree["decay"]["decay.vx"].array()
            self.decay_y = tree["decay"]["decay.vy"].array()
            self.decay_z = tree["decay"]["decay.vz"].array()

        # Load x-sec data 0-1 (GeV)
        self.xsecs = {
            14  : np.loadtxt("./DATA/Genie_xs/Genie_xs/numu_Ar_v3.04.00_t.AR23.20i.00.000.txt"),
            -14 : np.loadtxt("./DATA/Genie_xs/Genie_xs/numu_bar_Ar_v3.04.00_t.AR23.20i.00.000.txt"),
            12  : np.loadtxt("./DATA/Genie_xs/Genie_xs/nue_Ar_v3.04.00_t.AR23.20i.00.000.txt"),
            -12 : np.loadtxt("./DATA/Genie_xs/Genie_xs/nue_bar_Ar_v3.04.00_t.AR23.20i.00.000.txt")
        }
        self.xsec_scale = 1e-38 / 40


        # Target data
        rho = 1.38434 #g/cm^{3}
        n_a = 6.022e23 #molec/mol
        m_mol = 39.95 #g/mol
        A = 40
        volume_cm3 = 400*400*500

        self.n_targets = (rho * volume_cm3 * n_a * A) / (m_mol)

        # Useful dicts
        self.flavours = {
            14: (r"$\nu_\mu$",     "blue"),
            -14: (r"$\bar\nu_\mu$", "orange"),
            12: (r"$\nu_e$",       "green"),
            -12: (r"$\bar\nu_e$",   "red"),
        }

        self.parent_groups = {
            "mu" : ([13,-13], "blue"),
            "pi" : ([211,-211], "orange"),
            "K0L" : ( [130], "red"),
            "Kpm" : ([321, -321], "green")
        }


    def flavour_mask(self, pdgs=None):
        if pdgs is None:
            return np.ones_like(self.ntype, dtype=bool)
        if np.isscalar(pdgs):
            pdgs = [pdgs]
        return np.isin(self.ntype, pdgs)


    def parent_mask(self, p_pdgs=None):
        if p_pdgs is None:
            return np.ones_like(self.ptype, dtype=bool)
        if np.isscalar(p_pdgs):
            p_pdgs = [p_pdgs]
        return np.isin(self.ptype, p_pdgs)


    def energy_mask(self, energy_range = [0,1]):
        return (self.energies >= energy_range[0]) & (self.energies <= energy_range[1])


    def event_count(self, energy_range = [0,1], flavour=None, parent=None, extra_mask=None, pot=6.5e20):
        """
        energy_range: [min, max] array of float
        flavour: [pdg1, pdg2, ...] array of int or single int
        parents: [pdg1, pdg2, ...] array of int or single int
        """


        e_mask = self.energy_mask([energy_range[0], energy_range[1]])
        p_mask = self.parent_mask(parent)
        events = 0

        if flavour is None:
            flavour = list(self.flavours.keys())
        elif np.isscalar(flavour):
            flavour = [flavour]

        for f in flavour:
            f_mask = self.flavour_mask(f)
            mask = e_mask & f_mask & p_mask
            if extra_mask is not None:
                mask = mask & extra_mask


            values, edges = np.histogram(self.energies[mask], bins=self.bins, weights=self.tot_wgts[mask])
            bin_centres = (edges[1:] + edges[:-1]) / 2

            for centre, val in zip(bin_centres, values):
                xsection = np.interp(centre, self.xsecs[f][:,0], self.xsecs[f][:,1] * self.xsec_scale)
                evt = val * pot * xsection * self.n_targets
                events += evt

        return events
    
    

if __name__ == "__main__":

    data = Neutrinos_dk2nu("./DATA/NuBeam_production_BooNE_50m_I174000A_0.dk2nu.root")
    
    print(data.event_count(14))

