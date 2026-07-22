from __future__ import annotations

import math
from typing import Mapping

import numpy as np
import pyneb as pn
import pandas as pd

df = pd.read_csv('/home/nicolas/Documents/Research/PhD/JWST-Data/Medium_Resolution/Output_data/fluxes_for_abundance_all_filters.tsv', sep='\t')

ID_array = df["ID"]
z_array = df["redshift"]

#Maps for line labels and variable names
line_label_map = {
    "O3727": "[OII] 3727",
    "O4363": "[OIII] 4363",
    "H4861": "Hbeta 4861",
    "O4959": "[OIII] 4959",
    "O5007": "[OIII] 5007",
    "N5755": "[NII] 5755",
    "N6548": "[NII] 6548",
    "H6563": "Halpha 6563",
    "N6583": "[NII] 6583",
    "S6716": "[SII] 6716",
    "S6731": "[SII] 6731"
}




# ---------------------------------------------------------------------
# Observed integrated fluxes.
# Any consistent flux units are acceptable.
# ---------------------------------------------------------------------
flux: dict[int, float] = {
    3726: 27/2,   # [O II]
    3729: 27/2,   # [O II]
    4363: 32.3,    # [O III] auroral
    4861: 90.6,  # H beta
    4959: 191.9,  # [O III]
    5007: 507.2,  # [O III]
    5755: 6.36,   # [N II] auroral
    6548: 6.2,    # [N II]
    6563: 275.6,  # H alpha
    6584: 8,   # [N II]
    6716: 14.0,   # [S II]
    6731: 13.54,   # [S II]
}

flux_error: dict[int, float] = {
    3726: 3.63,  # replace with measured 1-sigma uncertainty
    3729: 3.63,
    4363: 1,
    4861: 3.8,
    4959: 2.1,
    5007: 1.9,
    5755: 4.06,
    6548: 5.97,
    6563: 5.94,
    6584: 5.8,
    6716: 4.8,
    6731: 4.84,
}

"""
def require_lines(data: Mapping[int, float], waves: list[int]) -> None:
    Check that all required lines are present and positive.
    missing = [wave for wave in waves if wave not in data]
    if missing:
        raise KeyError(f"Missing required wavelengths: {missing}")

    invalid = [wave for wave in waves if not np.isfinite(data[wave]) or data[wave] <= 0]
    if invalid:
        raise ValueError(f"Non-positive or invalid fluxes at: {invalid}")


require_lines(
    flux,
    [3726, 3729, 4363, 4861, 4959, 5007,
     5755, 6548, 6563, 6584, 6716, 6731],
)
"""

# ---------------------------------------------------------------------
# PyNeb atoms
#
# Spectroscopic notation:
# O2 means O+, O3 means O++, N2 means N+, etc.
# ---------------------------------------------------------------------
H1 = pn.RecAtom("H", 1)
O2 = pn.Atom("O", 2)
O3 = pn.Atom("O", 3)
N2 = pn.Atom("N", 2)
S2 = pn.Atom("S", 2)


def deredden_and_normalize(
    observed: Mapping[int, float],
    te: float,
    ne: float,
    law: str = "G03 LMC",
    r_v: float = 3.1,
) -> tuple[dict[int, float], pn.RedCorr, float]:
    """
    Deredden the fluxes using Halpha/Hbeta and normalize to Hbeta=100.

    The intrinsic Case-B Halpha/Hbeta ratio is calculated at the
    current Te and ne rather than fixed permanently at 2.86.
    """
    ha_hb_theory = (
        H1.getEmissivity(te, ne, label="3_2")
        / H1.getEmissivity(te, ne, label="4_2")
    )

    observed_over_theoretical = (
        observed[6563] / observed[4861]
    ) / ha_hb_theory

    rc = pn.RedCorr(law=law, R_V=r_v)
    rc.setCorr(
        obs_over_theo=observed_over_theoretical,
        wave1=6563.0,
        wave2=4861.0,
    )

    corrected = {
        wave: value * rc.getCorr(wave, rel_wave=4861.0)
        for wave, value in observed.items()
    }

    hb_scale = 100.0 / corrected[4861]
    normalized = {
        wave: value * hb_scale
        for wave, value in corrected.items()
    }

    return normalized, rc, float(ha_hb_theory)


# ---------------------------------------------------------------------
# Joint iteration because:
# - the Balmer decrement depends weakly on Te and ne;
# - [S II] density depends weakly on Te;
# - [O III] temperature depends weakly on ne.
# ---------------------------------------------------------------------
te_high = 10000.0
ne = 100.0

for iteration in range(20):
    intensity, reddening, intrinsic_ha_hb = deredden_and_normalize(
        flux,
        te=te_high,
        ne=ne,
        law="G03 LMC",
        r_v=2.1,
    )

    sii_ratio = intensity[6716] / intensity[6731]

    ne_new = S2.getTemDen(
        int_ratio=sii_ratio,
        tem=te_high,
        wave1=6716,
        wave2=6731,
    )

    oiii_ratio = (
        intensity[4959] + intensity[5007]
    ) / intensity[4363]

    te_new = O3.getTemDen(
        int_ratio=oiii_ratio,
        den=ne_new,
        to_eval="(L(4959) + L(5007)) / L(4363)",
    )

    if not np.isfinite(te_new) or not np.isfinite(ne_new):
        raise RuntimeError(
            "PyNeb could not invert one of the diagnostic ratios. "
            "Check the line fluxes, ratio orientation and diagnostic limits."
        )

    if abs(te_new - te_high) < 1.0 and abs(ne_new - ne) < 0.1:
        te_high = float(te_new)
        ne = float(ne_new)
        break

    te_high = float(te_new)
    ne = float(ne_new)

else:
    raise RuntimeError("Te-ne iteration did not converge.")


# Recompute the final reddening correction at the converged conditions.
intensity, reddening, intrinsic_ha_hb = deredden_and_normalize(
    flux,
    te=te_high,
    ne=ne,
)


# ---------------------------------------------------------------------
# Low-ionization-zone temperature from [N II].
# ---------------------------------------------------------------------
nii_ratio = (
    intensity[6548] + intensity[6584]
) / intensity[5755]

te_low = N2.getTemDen(
    int_ratio=nii_ratio,
    den=ne,
    to_eval="(L(6548) + L(6584)) / L(5755)",
)


# ---------------------------------------------------------------------
# Ionic abundances relative to H+.
# Input intensities are normalized to Hbeta=100.
# ---------------------------------------------------------------------
o_plus = O2.getIonAbundance(
    int_ratio=intensity[3726] + intensity[3729],
    tem=te_low,
    den=ne,
    to_eval="L(3726) + L(3729)",
    Hbeta=100.0,
)

o_double_plus = O3.getIonAbundance(
    int_ratio=intensity[4959] + intensity[5007],
    tem=te_high,
    den=ne,
    to_eval="L(4959) + L(5007)",
    Hbeta=100.0,
)

n_plus = N2.getIonAbundance(
    int_ratio=intensity[6548] + intensity[6584],
    tem=te_low,
    den=ne,
    to_eval="L(6548) + L(6584)",
    Hbeta=100.0,
)


# For an ordinary stellar H II region without significant He II:
oxygen_hydrogen = o_plus + o_double_plus
metallicity_oxygen = 12.0 + math.log10(oxygen_hydrogen)

# Standard approximation for H II regions:
nitrogen_oxygen = n_plus / o_plus
log_nitrogen_oxygen = math.log10(nitrogen_oxygen)


print(f"PyNeb version           = {pn.__version__}")
print(f"E(B-V)                  = {reddening.E_BV:.4f}")
print(f"c(Hbeta)                = {reddening.cHbeta:.4f}")
print(f"Intrinsic Halpha/Hbeta  = {intrinsic_ha_hb:.4f}")
print(f"Te([O III])             = {te_high:.0f} K")
print(f"Te([N II])              = {te_low:.0f} K")
print(f"ne([S II])              = {ne:.1f} cm^-3")
print(f"O+/H+                   = {o_plus:.4e}")
print(f"O++/H+                  = {o_double_plus:.4e}")
print(f"O/H                     = {oxygen_hydrogen:.4e}")
print(f"12 + log10(O/H)         = {metallicity_oxygen:.3f}")
print(f"log10(N/O)              = {log_nitrogen_oxygen:.3f}")
