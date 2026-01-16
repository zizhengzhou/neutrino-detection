# conda activate mkid_sim
import numpy as np
import matplotlib.pyplot as plt
import openmc.data
import xraylib
import os

if "OPENMC_CROSS_SECTIONS" not in os.environ:
    # os.environ['OPENMC_CROSS_SECTIONS'] = '/home/zzz/mkid_sim/cross_sections.xml'
    print("Warning: OPENMC_CROSS_SECTIONS not set. Script might fail.")

class DiffCrossSectionCalculator:
    def __init__(self):
        self.cache = {}
        self.library = openmc.data.DataLibrary.from_xml()

    def _load_isotope(self, element):
        iso_name = element
        if not iso_name:
            raise ValueError(f"No isotope mapping found for element {element}")

        if iso_name in self.cache:
            return self.cache[iso_name]
        
        try:
            path = [x['path'] for x in self.library.libraries if x['materials'][0] == iso_name][0]
            data = openmc.data.IncidentNeutron.from_hdf5(path)
            self.cache[iso_name] = data
            print(f"Loaded data for {element} ({iso_name})")
            return data
        except IndexError:
            raise FileNotFoundError(f"Data for {iso_name} not found in cross_sections.xml")

    def get_angular_pdf(self, element, E_n_eV, mu_val):
        """
        for epsicifc E_n @ mu ratio density:  p(mu)
        mu = 1 - 2E_R/E_R max
        """
        data = self._load_isotope(element)
        neutron_product = data[2].products[0]
        # neutron is the first product
        uncorrelated_dist = neutron_product.distribution[0]
        angle_container = uncorrelated_dist.angle

        if hasattr(angle_container, "energy"):
            # TODO change this part to intepotation
            idx = np.abs(angle_container.energy - E_n_eV).argmin()
            dist = angle_container.mu[idx]

            if hasattr(dist, "x") and hasattr(dist, "p"):  # Tabular
                return np.interp(mu_val, dist.x, dist.p)

            elif hasattr(dist, "coefficients"):  # Legendre
                # p(mu) = 0.5 + sum( (2l+1)/2 * al * Pl(mu) )
                coeffs = dist.coefficients
                val = 0.5  # l=0  (0.5 * 1.0)
                for l, a_l in enumerate(coeffs):
                    # openmc from l=1 
                    P_l = np.polynomial.legendre.Legendre.basis(l + 1)(mu_val)
                    val += (2 * (l + 1) + 1) / 2 * a_l * P_l
                return val

            elif isinstance(dist, openmc.data.angle_distribution.Uniform):
                return 0.5

        return 0.5  # Uniform as default

    def calculate_compound_diff_sigma(self, compound_dict, En_MeV_list, Er_keV_array):
        """
        Output: dSigma/dEr [barn/keV]
        """
        results = {}
        Er_eV = Er_keV_array * 1e3  # to eV

        for En_MeV in En_MeV_list:
            En_eV = En_MeV * 1e6
            total_ds_dEr = np.zeros_like(Er_eV)

            for element, count in compound_dict.items():
                data = self._load_isotope(element)
                A = data.atomic_weight_ratio

                # E_max (eV)
                r = 4 * A / ((1 + A) ** 2)
                E_max_eV = r * En_eV

                # total cross section: Noting that 294K is the lowest temperature data in database
                sigma_elastic = data[2].xs["294K"](En_eV)

                for i, E_rec in enumerate(Er_eV):
                    if E_rec >= E_max_eV:
                        continue

                    # E_rec = E_max * (1 - mu) / 2  => mu = 1 - 2*E_rec/E_max
                    mu = 1.0 - 2.0 * E_rec / E_max_eV
                    p_mu = self.get_angular_pdf(element, En_eV, mu)

                    # barn/eV -> 1e3 barn/keV
                    # dSigma/dEr = (Sigma_tot / E_max) * (p(mu) / 0.5)
                    #            = (Sigma_tot / E_max) * 2 * p(mu)
                    # when p(mu)=0.5 it becomes Sigma/E_max
                    val = (sigma_elastic / E_max_eV) * 2 * p_mu * 1e3

                    total_ds_dEr[i] += count * val

            results[En_MeV] = total_ds_dEr

        return results


def main():
    calc = DiffCrossSectionCalculator()

    materials = {
        "Si": {"Si28": 1},
        "Ge": {"Ge74": 1},
        "CaWO4": {"Ca40": 1, "W184": 1, "O16": 4},
        "Li2MoO4": {"Li7": 2, "Mo98": 1, "O16": 4},
        "PbWO4": {"Pb208": 1, "W184": 1, "O16": 4},
    }

    Er_array = np.logspace(-3, 4, 1000)
    En_list = [0.01, 0.1, 0.5, 1.0, 5.0]

    if not os.path.exists("figures_xsec"):
        os.makedirs("figures_xsec")

    for mat_name, formula in materials.items():
        try:
            print(f"Processing {mat_name}...")
            results = calc.calculate_compound_diff_sigma(formula, En_list, Er_array)
            plt.figure(figsize=(8, 6))
            colors = ["#1f8ab4", "#ff7f0e", "#2ca02c", "#e4ff19", "#ff00e1"]

            for idx, En in enumerate(En_list):
                y_vals = results[En]
                mask = y_vals > 0

                plt.plot(
                    Er_array[mask],
                    y_vals[mask],
                    label=f"$E_n = {En}$ MeV",
                    color=colors[idx],
                    linewidth=2,
                )

            plt.xscale("log")
            plt.yscale("log")

            plt.title(f"Neutron Differential Cross-section on {mat_name}", fontsize=14)
            plt.xlabel("Nuclear Recoil Energy $E_R$ [keV]", fontsize=12)
            plt.ylabel(r"$\frac{d\sigma}{dE_R}$ [barn / keV]", fontsize=12)

            plt.grid(True, which="both", ls="--", alpha=0.5)
            plt.legend(fontsize=10)
            plt.ylim(bottom=1e-4)

            filename = f"figures_xsec/xsec_{mat_name}.png"
            plt.savefig(filename, dpi=300)
            print(f"Saved {filename}")
            plt.close()

        except Exception as e:
            print(f"Skipped {mat_name} due to error: {e}")


if __name__ == "__main__":
    main()
