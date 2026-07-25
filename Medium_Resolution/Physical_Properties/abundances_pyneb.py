from __future__ import annotations
import numpy as np
import pyneb as pn
import pandas as pd
from typing import Mapping





# Input convention: wavelength in Angstrom -> integrated observed flux.
# 3727 means the unresolved sum [O II] 3726+3729.
PYNEB_LABEL = {
    3727: "O2_3727A+",
    4363: "O3_4363A",
    4861: "H1_4861A",
    4959: "O3_4959A",
    5007: "O3_5007A",
    6548: "N2_6548A",
    6563: "H1_6563A",
    6584: "N2_6584A",
    6716: "S2_6716A",
    6731: "S2_6731A",
}

H1 = pn.RecAtom("H", 1)
O2 = pn.Atom("O", 2)
O3 = pn.Atom("O", 3)
N2 = pn.Atom("N", 2)
S2 = pn.Atom("S", 2)




df = pd.read_csv("/home/nicolas/Documents/Research/PhD/JWST-Data/Medium_Resolution/Output_data/fluxes_for_abundance_all_filters.tsv", sep='\t')
ID_array = df["ID"]
z_array = df["redshift"]


def main():



    N_available_IDs = (4404 , 5759, 7892, 8083, 8113, 9422, 16625, 18846, 19519, 19606, 22251, 10009506, 10013618, 10013704, 10016186, 10035295, 10056849)

    for ID in ID_array:
        if ID in N_available_IDs:
            index = np.where(ID_array == ID)[0][0]
            
            flux = flux_from_df(index, df=df, source='JADES')
            print("======== Results for object with ID:", ID, "========================")
            try:
                result = analyze_pyneb_auto(flux)
            except RuntimeError:
            
                print("PyNeb analysis failed for this object.")
                continue
            except ValueError:
                print("Fluxes must be finite and positive, fails for this objects")
                continue

            
            print_pyneb_result(result)




def _check_input(flux: Mapping[int, float]) -> None:
    missing = sorted(set(PYNEB_LABEL) - set(flux))
    if missing:
        raise KeyError(f"Missing required wavelengths: {missing}")

    bad = [
        wave for wave in PYNEB_LABEL
        if not np.isfinite(flux[wave]) or flux[wave] <= 0
    ]
    if bad:
        raise ValueError(f"Fluxes must be finite and positive: {bad}")
        


def _deredden(
    flux: Mapping[int, float],
    te: float,
    ne: float,
    *,
    extinction_law: str,
    r_v: float,
    clip_negative_ebv: bool,
) -> tuple[dict[int, float], float, float, float, float]:
    """Return dereddened intensities normalized to Hbeta=100."""
    intrinsic_ha_hb = float( H1.getEmissivity(te, ne, label="3_2") / H1.getEmissivity(te, ne, label="4_2") )

    observed_over_theoretical = (flux[6563] / flux[4861] ) / intrinsic_ha_hb

    raw_rc = pn.RedCorr(law=extinction_law, R_V=r_v)
    raw_rc.setCorr( observed_over_theoretical, 6563.0, 4861.0)
    ebv_raw = float(raw_rc.E_BV)

    ebv_used = max(0.0, ebv_raw) if clip_negative_ebv else ebv_raw

    rc = pn.RedCorr(
        E_BV=ebv_used,
        law=extinction_law,
        R_V=r_v,
    )

    corrected = {
        wave: float(value * rc.getCorr(wave, rel_wave=4861.0))
        for wave, value in flux.items()
    }
    
    scale = 100.0 / corrected[4861]

    intensity =   {
        wave: value * scale
        for wave, value in corrected.items()
    }

    return (
        intensity,
        intrinsic_ha_hb,
        ebv_raw,
        ebv_used,
        float(rc.cHbeta),
    )


def _make_observation(intensity: Mapping[int, float]) -> pn.Observation:
    """Build the Observation used by Diagnostics.addDiagsFromObs."""
    obs = pn.Observation(corrected=True)

    # H I is used for reddening, not for the collisionally excited diagnostics.
    for wave, label in PYNEB_LABEL.items():
        if wave in (4861, 6563):
            continue

        obs.addLine(
            pn.EmissionLine(
                label=label,
                obsIntens=float(intensity[wave]),
                obsError=0.0,
                corrected=True,
            )
        )

    return obs


def analyze_pyneb_auto(flux: Mapping[int, float], *, extinction_law: str = "G03 LMC", r_v: float = 2.77, 
                       clip_negative_ebv: bool = True,
                       input_density: float = 100.0,
                       max_iterations: int = 8,
):
    """
    Automatic PyNeb analysis for the requested line set.

    Automatic parts
    ---------------
    1. Discover available PyNeb diagnostics from an Observation.
    2. Prefer [O III] 4363/(4959+5007) for Te.
    3. Prefer [S II] 6731/6716 for ne.
    4. Use Diagnostics.getCrossTemDen when both ratios are invertible.

    Necessary scientific fallbacks
    ------------------------------
    - If [S II] is not invertible, adopt fallback_density.
    - Since no low-zone auroral line is supplied, use
      Te(low) = 0.7 Te([O III]) + 3000 K.
    """
    _check_input(flux)

    te = 10000.0
    ne = float(input_density)
    warnings: list[str] = []

    available: tuple[str, ...] = ()
    te_diag: str | None = None
    ne_diag: str | None = None
    te_method = ""
    ne_method = ""
    sii_invalid = False

    for _ in range(max_iterations):
        intensity, _, _, _, _ = _deredden(
            flux,
            te,
            ne,
            extinction_law=extinction_law,
            r_v=r_v,
            clip_negative_ebv=clip_negative_ebv,
        )

        obs = _make_observation(intensity)
        diagnostics = pn.Diagnostics()
        diagnostics.addDiagsFromObs(obs)
        available = tuple(diagnostics.getDiagLabels())

        if "[OIII] 4363/5007+" in available:
            te_diag = "[OIII] 4363/5007+"
        elif "[OIII] 4363/5007" in available:
            te_diag = "[OIII] 4363/5007"
        else:
            raise RuntimeError(
                "No [O III] temperature diagnostic was discovered."
            )

        ne_diag = (
            "[SII] 6731/6716"
            if "[SII] 6731/6716" in available
            else None
        )

        te_new = np.nan
        ne_new = np.nan

        # Official PyNeb cross-convergence when both diagnostics are usable.
        if ne_diag is not None:
            te_cross, ne_cross = diagnostics.getCrossTemDen(
                te_diag,
                ne_diag,
                obs=obs,
                guess_tem=te,
                max_iter=20,
            )

            if np.isfinite(te_cross) and np.isfinite(ne_cross):
                te_new = float(te_cross)
                ne_new = float(ne_cross)
                te_method = f"cross-solution using {te_diag}"
                ne_method = f"cross-solution using {ne_diag}"

        # Density fallback if [S II] is missing or outside its physical range.
        if not np.isfinite(ne_new):
            if ne_diag is not None:
                sii_ratio = intensity[6731] / intensity[6716]
                ne_try = S2.getTemDen(
                    int_ratio=sii_ratio,
                    tem=te,
                    wave1=6731,
                    wave2=6716,
                    start_x=0.0,   # log10(1 cm^-3)
                    end_x=8.0,     # log10(1e8 cm^-3)
                )

                if np.isfinite(ne_try) and ne_try > 0:
                    ne_new = float(ne_try)
                    ne_method = ne_diag
                else:
                    sii_invalid = True

            if not np.isfinite(ne_new):
                ne_new = float(input_density)
                ne_method = f"assumed ne={input_density:g} cm^-3"

        # Temperature fallback: invert [O III] independently at the adopted ne.
        if not np.isfinite(te_new):
            if te_diag == "[OIII] 4363/5007+":
                observed_ratio = (
                    intensity[4363]
                    / (intensity[4959] + intensity[5007])
                )
                te_try = O3.getTemDen(
                    int_ratio=observed_ratio,
                    den=ne_new,
                    to_eval="L(4363)/(L(4959)+L(5007))",
                    start_x=np.log10(3000.0),
                    end_x=np.log10(50_000.0),
                )
            else:
                observed_ratio = intensity[4363] / intensity[5007]
                te_try = O3.getTemDen(
                    int_ratio=observed_ratio,
                    den=ne_new,
                    wave1=4363,
                    wave2=5007,
                    start_x=np.log10(3000.0),
                    end_x=np.log10(50_000.0),
                )

            if not np.isfinite(te_try) or te_try <= 0:
                raise RuntimeError(
                    f"{te_diag} cannot be inverted between 3000 and 50000 K."
                )

            te_new = float(te_try)
            te_method = te_diag

        if abs(te_new - te) < 1.0 and abs(ne_new - ne) < 0.1:
            te, ne = te_new, ne_new
            break

        te, ne = te_new, ne_new

    # Final dereddening at the converged conditions.
    (
        intensity,
        intrinsic_ha_hb,
        ebv_raw,
        ebv_used,
        c_hbeta,
    ) = _deredden(
        flux,
        te,
        ne,
        extinction_law=extinction_law,
        r_v=r_v,
        clip_negative_ebv=clip_negative_ebv,
    )

    if sii_invalid:
        ratio = intensity[6731] / intensity[6716]
        low_limit = S2.getLowDensRatio(
            tem=te,
            wave1=6731,
            wave2=6716,
        )
        high_limit = S2.getHighDensRatio(
            tem=te,
            wave1=6731,
            wave2=6716,
        )
        warnings.append(
            f"[S II] 6731/6716={ratio:.4g} is outside the approximate "
            f"PyNeb range {min(low_limit, high_limit):.4g}–"
            f"{max(low_limit, high_limit):.4g}; ne was assumed."
        )

    if clip_negative_ebv and ebv_raw < 0:
        warnings.append(
            f"The Balmer decrement gives E(B-V)={ebv_raw:.3f}; "
            "the value used was clipped to zero."
        )

    # Quality-control checks; they do not alter the input fluxes.
    o3_branching = intensity[5007] / intensity[4959]
    if not 2.7 <= o3_branching <= 3.3:
        warnings.append(
            f"[O III] 5007/4959={o3_branching:.3f}, not close to ~2.98."
        )

    n2_branching = intensity[6584] / intensity[6548]
    if not 2.7 <= n2_branching <= 3.3:
        warnings.append(
            f"[N II] 6584/6548={n2_branching:.3f}, not close to ~3; "
            "N+/H+ and N/O may be unreliable."
        )

    # No [N II] 5755 or [O II] 7320+7330 is supplied, so Te(low) is inferred.
    te_low = 0.7 * te + 3000.0

    o_plus = float(
        O2.getIonAbundance(
            int_ratio=intensity[3727],
            tem=te_low,
            den=ne,
            to_eval="L(3726)+L(3729)",
            Hbeta=100.0,
        )
    )

    o_double_plus = float(
        O3.getIonAbundance(
            int_ratio=intensity[4959] + intensity[5007],
            tem=te,
            den=ne,
            to_eval="L(4959)+L(5007)",
            Hbeta=100.0,
        )
    )

    n_plus = float(
        N2.getIonAbundance(
            int_ratio=intensity[6548] + intensity[6584],
            tem=te_low,
            den=ne,
            to_eval="L(6548)+L(6584)",
            Hbeta=100.0,
        )
    )

    oxygen = o_plus + o_double_plus

    return {
        "PyNeb_version": pn.__version__,
        "available_diagnostics": available,
        "temperature_diagnostic": te_diag,
        "density_diagnostic": ne_diag,
        "temperature_method": te_method,
        "density_method": ne_method,
        "low_temperature_method": "Te(low)=0.7 Te(OIII)+3000 K",
        "intrinsic_Ha_Hb": intrinsic_ha_hb,
        "E_BV_raw": ebv_raw,
        "E_BV_used": ebv_used,
        "cHbeta_used": c_hbeta,
        "Te_OIII_K": te,
        "Te_low_K": te_low,
        "ne_cm-3": ne,
        "O_plus_H": o_plus,
        "O_double_plus_H": o_double_plus,
        "O_H": oxygen,
        "12_log_O_H": 12.0 + np.log10(oxygen),
        "N_plus_H": n_plus,
        "log_N_O": np.log10(n_plus / o_plus),
        "warnings": tuple(warnings),
        "intensity_Hbeta_100": intensity,
    }


def print_pyneb_result(result: Mapping[str, object]):
    print(f"PyNeb version           = {result['PyNeb_version']}")
    print(f"Available diagnostics   = {list(result['available_diagnostics'])}")
    print(f"Temperature diagnostic  = {result['temperature_diagnostic']}")
    print(f"Density diagnostic      = {result['density_diagnostic']}")
    print(f"Temperature method      = {result['temperature_method']}")
    print(f"Density method          = {result['density_method']}")
    print(f"E(B-V), raw / used      = {result['E_BV_raw']:.4f} / {result['E_BV_used']:.4f}")
    print(f"c(Hbeta), used          = {result['cHbeta_used']:.8f}")
    print(f"Intrinsic Halpha/Hbeta  = {result['intrinsic_Ha_Hb']:.4f}")
    print(f"Te([O III])             = {result['Te_OIII_K']:.0f} K")
    print(f"Te(low)                 = {result['Te_low_K']:.0f} K")
    print(f"ne                       = {result['ne_cm-3']:.1f} cm^-3")
    print(f"O+/H+                   = {result['O_plus_H']:.4e}")
    print(f"O++/H+                  = {result['O_double_plus_H']:.4e}")
    print(f"O/H                     = {result['O_H']:.4e}")
    print(f"12 + log10(O/H)         = {result['12_log_O_H']:.3f}")
    print(f"N+/H+                   = {result['N_plus_H']:.4e}")
    print(f"log10(N/O)              = {result['log_N_O']:.3f}")

    warnings = result["warnings"]
    if warnings:
        print("Warnings:")
        for warning in warnings:
            print(f"  - {warning}")


def flux_from_df(index, df=df, source='MAST', printWaves=False):

    line_label_map = {
    "3727": "[OII] 3727",
    "4363": "[OIII] 4363",
    "4861": "Hbeta 4861",    
    "4959": "[OIII] 4959",
    "5007": "[OIII] 5007",
    "5755": "[NII] 5755",
    "6548": "[NII] 6548",
    "6563": "Halpha 6563",    
    "6584": "[NII] 6583",
    "6716": "[SII] 6716",
    "6731": "[SII] 6731"   
    }

    keys = ['ID', 'redshift']
    cols_g235m = [f'Flux ({label}) {source} G235M' for label in line_label_map.values()]
    cols_g395m = [f'Flux ({label}) {source} G395M' for label in line_label_map.values()]
    cols_g395h = [f'Flux ({label}) {source} G395H' for label in line_label_map.values()]
    cols_errs_g235m = [f'F_err({label}) {source} G235M' for label in line_label_map.values()]
    cols_errs_g395m = [f'F_err({label}) {source} G395M' for label in line_label_map.values()]
    cols_errs_g395h = [f'F_err({label}) {source} G395H' for label in line_label_map.values()]

    df_sub_g235m = df[keys + cols_g235m + cols_errs_g235m].copy() 
    df_sub_g395m = df[keys + cols_g395m + cols_errs_g395m].copy()
    df_sub_g395h = df[keys + cols_g395h + cols_errs_g395h].copy()


    df_merged = (
            df_sub_g235m[keys].merge(df_sub_g395m[keys], on=keys, how='inner').merge(df_sub_g395h[keys], on=keys, how='inner')
        )
    

    for line, label in line_label_map.items():
        
        df_merged[line] = df_sub_g235m[f'Flux ({label}) {source} G235M'].combine_first(df_sub_g395m[f'Flux ({label}) {source} G395M']).combine_first(df_sub_g395h[f'Flux ({label}) {source} G395H'])
        df_merged[f'F_err({line})'] = df_sub_g235m[f'F_err({label}) {source} G235M'].combine_first(df_sub_g395m[f'F_err({label}) {source} G395M']).combine_first(df_sub_g395h[f'F_err({label}) {source} G395H'])

        
    flux ={        
            3727: df_merged['3727'].iloc[index],   # [O II]
            4363: df_merged['4363'].iloc[index],    # [O III] auroral
            4861: df_merged['4861'].iloc[index], # H beta
            4959: df_merged['4959'].iloc[index], # [O III]
            5007: df_merged['5007'].iloc[index], # [O III]
            6548: df_merged['6548'].iloc[index], # [N II]
            6563: df_merged['6563'].iloc[index], # H alpha
            6584: df_merged['6584'].iloc[index], # [N II]
            6716: df_merged['6716'].iloc[index], # [S II]
            6731: df_merged['6731'].iloc[index], # [S II]
        }
    
    if printWaves:
        print("Line values for object with ID:", df_merged['ID'].iloc[index])
        for line in flux:
            print(line, '\t', flux[line])
        print("                                               ")
    
    
    return flux





main()