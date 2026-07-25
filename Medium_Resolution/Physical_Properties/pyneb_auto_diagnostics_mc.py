from __future__ import annotations

"""
PyNeb semi-direct abundance analysis with Monte-Carlo uncertainties
and explicit treatment of low-S/N lines as upper limits.

Input convention
----------------
- Keys are rest-frame wavelengths in Angstrom.
- Fluxes and 1-sigma errors use any common integrated-flux unit.
- 3727 means the unresolved sum [O II] 3726+3729.
- A line is automatically classed as a detection when flux/error >=
  snr_threshold and flux > 0, unless status_override changes it.
- For a nondetection, an explicit upper limit can be supplied. If it is
  omitted, upper_limit_sigma * error is used.

Important statistical convention
--------------------------------
Upper-limit lines are NOT randomly drawn as if they were detections.
They are propagated as one-sided bounds. Sampling a nondetection from a
half-normal or uniform distribution would impose an additional prior and
can bias the result.
"""

import math
from collections import Counter
from dataclasses import dataclass
from typing import Literal, Mapping

import numpy as np
import pyneb as pn


Wave = int
LineStatus = Literal["detection", "upper_limit", "missing"]

EXPECTED_WAVES: tuple[Wave, ...] = (
    3727, 4363, 4861, 4959, 5007,
    6548, 6563, 6584, 6716, 6731,
)

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

# Atomic branching ratios used to tie fixed-ratio doublets.
O3_5007_4959 = 2.98
N2_6584_6548 = 2.96

H1 = pn.RecAtom("H", 1)
O2 = pn.Atom("O", 2)
O3 = pn.Atom("O", 3)
N2 = pn.Atom("N", 2)
S2 = pn.Atom("S", 2)


@dataclass(frozen=True)
class Measurement:
    flux: float | None
    error: float | None
    status: LineStatus
    snr: float | None
    upper_limit: float | None


@dataclass(frozen=True)
class FixedDoubletModel:
    """Latent amplitude model: weak line=A, strong line=ratio*A."""

    weak_wave: Wave
    strong_wave: Wave
    ratio: float
    amplitude: float
    amplitude_error: float

    def central_fluxes(self) -> dict[Wave, float]:
        return {
            self.weak_wave: self.amplitude,
            self.strong_wave: self.ratio * self.amplitude,
        }

    def draw(self, rng: np.random.Generator) -> dict[Wave, float]:
        amplitude = float(rng.normal(self.amplitude, self.amplitude_error))
        return {
            self.weak_wave: amplitude,
            self.strong_wave: self.ratio * amplitude,
        }


def classify_lines(
    flux: Mapping[Wave, float],
    flux_error: Mapping[Wave, float],
    *,
    upper_limits: Mapping[Wave, float] | None = None,
    status_override: Mapping[Wave, LineStatus] | None = None,
    snr_threshold: float = 3.0,
    upper_limit_sigma: float = 3.0,
) -> dict[Wave, Measurement]:
    """Classify supplied lines as detections, upper limits, or missing."""
    if snr_threshold <= 0 or upper_limit_sigma <= 0:
        raise ValueError("snr_threshold and upper_limit_sigma must be positive.")

    upper_limits = dict(upper_limits or {})
    status_override = dict(status_override or {})
    out: dict[Wave, Measurement] = {}

    for wave in EXPECTED_WAVES:
        f_raw = flux.get(wave)
        e_raw = flux_error.get(wave)

        f = float(f_raw) if f_raw is not None and np.isfinite(f_raw) else None
        e = float(e_raw) if e_raw is not None and np.isfinite(e_raw) else None

        if e is not None and e <= 0:
            raise ValueError(f"Uncertainty for line {wave} must be positive.")

        snr = None if f is None or e is None else f / e
        override = status_override.get(wave)

        if override is not None and override not in {
            "detection", "upper_limit", "missing"
        }:
            raise ValueError(
                f"Invalid status_override for {wave}: {override!r}."
            )

        if override == "missing":
            status: LineStatus = "missing"
        elif override == "detection":
            if f is None or e is None or f <= 0:
                raise ValueError(
                    f"Line {wave} was forced to detection but has invalid flux/error."
                )
            status = "detection"
        elif override == "upper_limit":
            status = "upper_limit"
        elif f is None or e is None:
            status = "missing"
        elif f > 0 and snr is not None and snr >= snr_threshold:
            status = "detection"
        else:
            status = "upper_limit"

        limit: float | None = None
        if status == "upper_limit":
            if wave in upper_limits:
                limit = float(upper_limits[wave])
            elif e is not None:
                # Common survey convention for a nondetection. Replace this
                # with a spectrum-integrated limit whenever available.
                limit = upper_limit_sigma * e

            if limit is None or not np.isfinite(limit) or limit <= 0:
                raise ValueError(
                    f"Line {wave} is an upper limit but no valid limit is available."
                )

        out[wave] = Measurement(
            flux=f,
            error=e,
            status=status,
            snr=snr,
            upper_limit=limit,
        )

    return out


def _fit_fixed_doublet(
    measurements: Mapping[Wave, Measurement],
    weak_wave: Wave,
    strong_wave: Wave,
    ratio: float,
) -> FixedDoubletModel | None:
    """Weighted fit of a fixed-ratio doublet using detected components."""
    weak = measurements[weak_wave]
    strong = measurements[strong_wave]

    num = 0.0
    den = 0.0

    if weak.status == "detection":
        assert weak.flux is not None and weak.error is not None
        num += weak.flux / weak.error**2
        den += 1.0 / weak.error**2

    if strong.status == "detection":
        assert strong.flux is not None and strong.error is not None
        num += ratio * strong.flux / strong.error**2
        den += ratio**2 / strong.error**2

    if den == 0:
        return None

    amplitude = num / den
    amplitude_error = 1.0 / math.sqrt(den)

    return FixedDoubletModel(
        weak_wave=weak_wave,
        strong_wave=strong_wave,
        ratio=ratio,
        amplitude=float(amplitude),
        amplitude_error=float(amplitude_error),
    )


def _prepare_central_fluxes(
    measurements: Mapping[Wave, Measurement],
    *,
    tie_fixed_doublets: bool,
) -> tuple[dict[Wave, float], dict[str, FixedDoubletModel | None]]:
    """Build the detected-line flux dictionary used in the central solution."""
    working: dict[Wave, float] = {}

    # Independent lines.
    for wave in (3727, 4363, 4861, 6563, 6716, 6731):
        m = measurements[wave]
        if m.status == "detection":
            assert m.flux is not None
            working[wave] = m.flux

    models: dict[str, FixedDoubletModel | None] = {
        "O3": None,
        "N2": None,
    }

    if tie_fixed_doublets:
        models["O3"] = _fit_fixed_doublet(
            measurements, 4959, 5007, O3_5007_4959
        )
        models["N2"] = _fit_fixed_doublet(
            measurements, 6548, 6584, N2_6584_6548
        )

        for model in models.values():
            if model is not None:
                working.update(model.central_fluxes())
    else:
        for wave in (4959, 5007, 6548, 6584):
            m = measurements[wave]
            if m.status == "detection":
                assert m.flux is not None
                working[wave] = m.flux

        # PyNeb temperature and abundance expressions use both [O III] lines.
        if 5007 in working and 4959 not in working:
            working[4959] = working[5007] / O3_5007_4959
        elif 4959 in working and 5007 not in working:
            working[5007] = O3_5007_4959 * working[4959]

        if 6584 in working and 6548 not in working:
            working[6548] = working[6584] / N2_6584_6548
        elif 6548 in working and 6584 not in working:
            working[6584] = N2_6584_6548 * working[6548]

    return working, models


def _draw_detected_fluxes(
    measurements: Mapping[Wave, Measurement],
    models: Mapping[str, FixedDoubletModel | None],
    rng: np.random.Generator,
    *,
    tie_fixed_doublets: bool,
) -> dict[Wave, float]:
    """Draw one realization of detected fluxes without redrawing negatives."""
    draw: dict[Wave, float] = {}

    for wave in (3727, 4363, 4861, 6563, 6716, 6731):
        m = measurements[wave]
        if m.status == "detection":
            assert m.flux is not None and m.error is not None
            draw[wave] = float(rng.normal(m.flux, m.error))

    if tie_fixed_doublets:
        for model in models.values():
            if model is not None:
                draw.update(model.draw(rng))
    else:
        for wave in (4959, 5007, 6548, 6584):
            m = measurements[wave]
            if m.status == "detection":
                assert m.flux is not None and m.error is not None
                draw[wave] = float(rng.normal(m.flux, m.error))

        if 5007 in draw and 4959 not in draw:
            draw[4959] = draw[5007] / O3_5007_4959
        elif 4959 in draw and 5007 not in draw:
            draw[5007] = O3_5007_4959 * draw[4959]

        if 6584 in draw and 6548 not in draw:
            draw[6548] = draw[6584] / N2_6584_6548
        elif 6548 in draw and 6584 not in draw:
            draw[6584] = N2_6584_6548 * draw[6548]

    return draw


def _deredden(
    flux: Mapping[Wave, float],
    te: float,
    ne: float,
    *,
    extinction_law: str,
    r_v: float,
    clip_negative_ebv: bool,
) -> tuple[dict[Wave, float], float, float, float, float]:
    """Deredden supplied fluxes and normalize them to Hbeta=100."""
    if 4861 not in flux or 6563 not in flux:
        raise ValueError("Detected Hbeta and Halpha are required for reddening.")
    if flux[4861] <= 0 or flux[6563] <= 0:
        raise ValueError("Hbeta and Halpha must be positive in this realization.")

    intrinsic_ha_hb = float(
        H1.getEmissivity(te, ne, label="3_2")
        / H1.getEmissivity(te, ne, label="4_2")
    )

    obs_over_theo = (flux[6563] / flux[4861]) / intrinsic_ha_hb
    raw_rc = pn.RedCorr(law=extinction_law, R_V=r_v)
    raw_rc.setCorr(obs_over_theo, 6563.0, 4861.0)
    ebv_raw = float(raw_rc.E_BV)

    ebv_used = max(0.0, ebv_raw) if clip_negative_ebv else ebv_raw
    rc = pn.RedCorr(E_BV=ebv_used, law=extinction_law, R_V=r_v)

    corrected = {
        wave: float(value * rc.getCorr(wave, rel_wave=4861.0))
        for wave, value in flux.items()
    }
    scale = 100.0 / corrected[4861]
    intensity = {wave: value * scale for wave, value in corrected.items()}

    return intensity, intrinsic_ha_hb, ebv_raw, ebv_used, float(rc.cHbeta)


def _make_observation(intensity: Mapping[Wave, float]) -> pn.Observation:
    obs = pn.Observation(corrected=True)
    for wave, value in intensity.items():
        if wave in (4861, 6563) or wave not in PYNEB_LABEL:
            continue
        if not np.isfinite(value) or value <= 0:
            continue
        obs.addLine(
            pn.EmissionLine(
                label=PYNEB_LABEL[wave],
                obsIntens=float(value),
                obsError=0.0,
                corrected=True,
            )
        )
    return obs


def _solve_direct_realization(
    flux: Mapping[Wave, float],
    *,
    extinction_law: str,
    r_v: float,
    clip_negative_ebv: bool,
    fallback_density: float,
    max_iterations: int,
    te_low_scatter_offset: float = 0.0,
) -> dict[str, object]:
    """Solve one realization when [O III] 4363 is detected."""
    required = (4363, 4861, 4959, 5007, 6563)
    missing = [wave for wave in required if wave not in flux]
    if missing:
        raise ValueError(f"Direct solution is missing required detections: {missing}")
    if any(not np.isfinite(flux[w]) or flux[w] <= 0 for w in required):
        raise ValueError("A required direct-method line is non-positive.")

    te = 15_000.0
    ne = float(fallback_density)
    density_was_assumed = True
    available: tuple[str, ...] = ()
    te_diag: str | None = None
    ne_diag: str | None = None
    te_method = ""
    ne_method = f"assumed ne={fallback_density:g} cm^-3"

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
            raise RuntimeError("No [O III] temperature diagnostic was found.")

        ne_diag = (
            "[SII] 6731/6716"
            if "[SII] 6731/6716" in available
            else None
        )

        te_new = np.nan
        ne_new = np.nan

        if ne_diag is not None:
            te_cross, ne_cross = diagnostics.getCrossTemDen(
                te_diag,
                ne_diag,
                obs=obs,
                guess_tem=te,
                max_iter=20,
            )
            if (
                np.isfinite(te_cross)
                and np.isfinite(ne_cross)
                and te_cross > 0
                and ne_cross > 0
            ):
                te_new = float(te_cross)
                ne_new = float(ne_cross)
                density_was_assumed = False
                te_method = f"cross-solution using {te_diag}"
                ne_method = f"cross-solution using {ne_diag}"

        # If the coupled solution failed, independently test [S II].
        if not np.isfinite(ne_new):
            if 6716 in intensity and 6731 in intensity:
                sii_ratio = intensity[6731] / intensity[6716]
                ne_try = S2.getTemDen(
                    int_ratio=sii_ratio,
                    tem=te,
                    wave1=6731,
                    wave2=6716,
                    start_x=0.0,
                    end_x=8.0,
                )
                if np.isfinite(ne_try) and ne_try > 0:
                    ne_new = float(ne_try)
                    density_was_assumed = False
                    ne_method = "[SII] 6731/6716"

            if not np.isfinite(ne_new):
                ne_new = float(fallback_density)
                density_was_assumed = True
                ne_method = f"assumed ne={fallback_density:g} cm^-3"

        # Independently derive Te at the accepted density if needed.
        if not np.isfinite(te_new):
            ratio = intensity[4363] / (intensity[4959] + intensity[5007])
            te_try = O3.getTemDen(
                int_ratio=ratio,
                den=ne_new,
                to_eval="L(4363)/(L(4959)+L(5007))",
                start_x=np.log10(3000.0),
                end_x=np.log10(50_000.0),
            )
            if not np.isfinite(te_try) or te_try <= 0:
                raise RuntimeError(
                    "[O III] 4363/(4959+5007) is outside the invertible range."
                )
            te_new = float(te_try)
            te_method = "[OIII] 4363/(4959+5007)"

        if abs(te_new - te) < 1.0 and abs(ne_new - ne) < 0.1:
            te, ne = te_new, ne_new
            break
        te, ne = te_new, ne_new

    intensity, intrinsic, ebv_raw, ebv_used, c_hbeta = _deredden(
        flux,
        te,
        ne,
        extinction_law=extinction_law,
        r_v=r_v,
        clip_negative_ebv=clip_negative_ebv,
    )

    te_low = 0.7 * te + 3000.0 + te_low_scatter_offset
    if te_low <= 0:
        raise ValueError("The sampled low-zone temperature is non-positive.")

    result: dict[str, object] = {
        "available_diagnostics": available,
        "temperature_diagnostic": te_diag,
        "density_diagnostic": ne_diag,
        "temperature_method": te_method,
        "density_method": ne_method,
        "density_was_assumed": density_was_assumed,
        "intrinsic_Ha_Hb": intrinsic,
        "E_BV_raw": ebv_raw,
        "E_BV_used": ebv_used,
        "cHbeta_used": c_hbeta,
        "Te_OIII_K": te,
        "Te_low_K": te_low,
        "ne_cm-3": ne,
        "intensity_Hbeta_100": intensity,
    }

    o_double_plus = float(
        O3.getIonAbundance(
            int_ratio=intensity[4959] + intensity[5007],
            tem=te,
            den=ne,
            to_eval="L(4959)+L(5007)",
            Hbeta=100.0,
        )
    )
    result["O_double_plus_H"] = o_double_plus

    if 3727 in intensity and intensity[3727] > 0:
        o_plus = float(
            O2.getIonAbundance(
                int_ratio=intensity[3727],
                tem=te_low,
                den=ne,
                to_eval="L(3726)+L(3729)",
                Hbeta=100.0,
            )
        )
        oxygen = o_plus + o_double_plus
        result.update(
            {
                "O_plus_H": o_plus,
                "O_H": oxygen,
                "12_log_O_H": 12.0 + math.log10(oxygen),
            }
        )

        if 6548 in intensity and 6584 in intensity:
            n_plus = float(
                N2.getIonAbundance(
                    int_ratio=intensity[6548] + intensity[6584],
                    tem=te_low,
                    den=ne,
                    to_eval="L(6548)+L(6584)",
                    Hbeta=100.0,
                )
            )
            if n_plus > 0 and o_plus > 0:
                result.update(
                    {
                        "N_plus_H": n_plus,
                        "log_N_O": math.log10(n_plus / o_plus),
                    }
                )

    return result


def _corrected_limit_intensity(
    wave: Wave,
    upper_limit: float,
    result: Mapping[str, object],
    hbeta_flux: float,
    *,
    extinction_law: str,
    r_v: float,
) -> float:
    rc = pn.RedCorr(
        E_BV=float(result["E_BV_used"]),
        law=extinction_law,
        R_V=r_v,
    )
    return float(
        upper_limit
        * rc.getCorr(wave, rel_wave=4861.0)
        / hbeta_flux
        * 100.0
    )


def _add_one_sided_bounds(
    result: dict[str, object],
    measurements: Mapping[Wave, Measurement],
    hbeta_flux: float,
    *,
    extinction_law: str,
    r_v: float,
) -> None:
    """Append abundance bounds implied by nondetected [O II]/[N II]."""
    if "Te_OIII_K" not in result or "Te_low_K" not in result:
        return

    te = float(result["Te_OIII_K"])
    te_low = float(result["Te_low_K"])
    ne = float(result["ne_cm-3"])

    # [O II] upper limit -> O+ upper limit and total O/H interval.
    m_o2 = measurements[3727]
    if m_o2.status == "upper_limit" and m_o2.upper_limit is not None:
        i_o2_upper = _corrected_limit_intensity(
            3727,
            m_o2.upper_limit,
            result,
            hbeta_flux,
            extinction_law=extinction_law,
            r_v=r_v,
        )
        o_plus_upper = float(
            O2.getIonAbundance(
                int_ratio=i_o2_upper,
                tem=te_low,
                den=ne,
                to_eval="L(3726)+L(3729)",
                Hbeta=100.0,
            )
        )
        result["O_plus_H_upper"] = o_plus_upper
        if "O_double_plus_H" in result:
            o2p = float(result["O_double_plus_H"])
            result["O_H_lower"] = o2p
            result["O_H_upper"] = o2p + o_plus_upper
            result["12_log_O_H_lower"] = 12.0 + math.log10(o2p)
            result["12_log_O_H_upper"] = 12.0 + math.log10(o2p + o_plus_upper)

    # [N II] upper limits -> N+/H+ and N/O upper limit.
    if "O_plus_H" not in result:
        return

    weak = measurements[6548]
    strong = measurements[6584]
    amplitude_limits: list[float] = []

    if weak.status == "upper_limit" and weak.upper_limit is not None:
        amplitude_limits.append(weak.upper_limit)
    if strong.status == "upper_limit" and strong.upper_limit is not None:
        amplitude_limits.append(strong.upper_limit / N2_6584_6548)

    # Only construct an N II bound when neither component was used as a detection.
    if amplitude_limits and weak.status != "detection" and strong.status != "detection":
        amplitude_upper = min(amplitude_limits)
        total_flux_upper = (1.0 + N2_6584_6548) * amplitude_upper

        # Use the 6584 correction for the tied doublet; the difference from
        # using 6548 is negligible because the wavelengths are adjacent.
        i_n2_upper = _corrected_limit_intensity(
            6584,
            total_flux_upper,
            result,
            hbeta_flux,
            extinction_law=extinction_law,
            r_v=r_v,
        )
        n_plus_upper = float(
            N2.getIonAbundance(
                int_ratio=i_n2_upper,
                tem=te_low,
                den=ne,
                to_eval="L(6548)+L(6584)",
                Hbeta=100.0,
            )
        )
        result["N_plus_H_upper"] = n_plus_upper
        o_plus = float(result["O_plus_H"])
        if o_plus > 0 and n_plus_upper > 0:
            result["log_N_O_upper"] = math.log10(n_plus_upper / o_plus)


def _solve_temperature_upper_limit(
    flux: Mapping[Wave, float],
    upper_4363: float,
    *,
    extinction_law: str,
    r_v: float,
    clip_negative_ebv: bool,
    fallback_density: float,
    max_iterations: int,
) -> dict[str, object]:
    """
    Convert an [O III] 4363 upper limit into Te upper and O/H lower bounds.

    Since 4363/(4959+5007) increases with temperature, an upper limit on
    4363 gives an upper limit on Te. For fixed nebular-line fluxes, this in
    turn gives a lower bound on the ionic abundance.
    """
    required = (4861, 4959, 5007, 6563)
    if any(w not in flux or flux[w] <= 0 for w in required):
        raise ValueError(
            "A 4363-limit solution requires positive Halpha, Hbeta, 4959 and 5007."
        )

    te_upper = 20_000.0
    ne = float(fallback_density)
    density_method = f"assumed ne={fallback_density:g} cm^-3"

    for _ in range(max_iterations):
        # Include the upper-limit value only to apply the reddening correction.
        temp_flux = dict(flux)
        temp_flux[4363] = upper_4363
        intensity, intrinsic, ebv_raw, ebv_used, c_hbeta = _deredden(
            temp_flux,
            te_upper,
            ne,
            extinction_law=extinction_law,
            r_v=r_v,
            clip_negative_ebv=clip_negative_ebv,
        )

        if 6716 in intensity and 6731 in intensity:
            ne_try = S2.getTemDen(
                int_ratio=intensity[6731] / intensity[6716],
                tem=te_upper,
                wave1=6731,
                wave2=6716,
                start_x=0.0,
                end_x=8.0,
            )
            if np.isfinite(ne_try) and ne_try > 0:
                ne_new = float(ne_try)
                density_method = "[SII] 6731/6716"
            else:
                ne_new = float(fallback_density)
        else:
            ne_new = float(fallback_density)

        ratio_upper = intensity[4363] / (intensity[4959] + intensity[5007])
        te_new = O3.getTemDen(
            int_ratio=ratio_upper,
            den=ne_new,
            to_eval="L(4363)/(L(4959)+L(5007))",
            start_x=np.log10(3000.0),
            end_x=np.log10(50_000.0),
        )
        if not np.isfinite(te_new) or te_new <= 0:
            raise RuntimeError(
                "The [O III] 4363 upper limit does not map into 3000-50000 K."
            )

        if abs(float(te_new) - te_upper) < 1.0 and abs(ne_new - ne) < 0.1:
            te_upper = float(te_new)
            ne = ne_new
            break
        te_upper = float(te_new)
        ne = ne_new

    te_low_upper = 0.7 * te_upper + 3000.0

    result: dict[str, object] = {
        "temperature_method": "[OIII] 4363 upper limit",
        "density_method": density_method,
        "intrinsic_Ha_Hb": intrinsic,
        "E_BV_raw": ebv_raw,
        "E_BV_used": ebv_used,
        "cHbeta_used": c_hbeta,
        "Te_OIII_upper_K": te_upper,
        "Te_low_upper_K": te_low_upper,
        "ne_cm-3": ne,
        "intensity_Hbeta_100": intensity,
    }

    o_double_plus_lower = float(
        O3.getIonAbundance(
            int_ratio=intensity[4959] + intensity[5007],
            tem=te_upper,
            den=ne,
            to_eval="L(4959)+L(5007)",
            Hbeta=100.0,
        )
    )
    result["O_double_plus_H_lower"] = o_double_plus_lower

    if 3727 in intensity and intensity[3727] > 0:
        o_plus_lower = float(
            O2.getIonAbundance(
                int_ratio=intensity[3727],
                tem=te_low_upper,
                den=ne,
                to_eval="L(3726)+L(3729)",
                Hbeta=100.0,
            )
        )
        oxygen_lower = o_plus_lower + o_double_plus_lower
        result.update(
            {
                "O_plus_H_lower": o_plus_lower,
                "O_H_lower": oxygen_lower,
                "12_log_O_H_lower": 12.0 + math.log10(oxygen_lower),
            }
        )

    return result


def _summarize_samples(samples: list[dict[str, object]]) -> dict[str, dict[str, float]]:
    keys: set[str] = set()
    for sample in samples:
        for key, value in sample.items():
            if isinstance(value, (int, float, np.floating)) and np.isfinite(value):
                keys.add(key)

    summary: dict[str, dict[str, float]] = {}
    for key in sorted(keys):
        values = np.asarray(
            [
                float(s[key])
                for s in samples
                if key in s
                and isinstance(s[key], (int, float, np.floating))
                and np.isfinite(s[key])
            ],
            dtype=float,
        )
        if values.size == 0:
            continue
        q16, q50, q84 = np.percentile(values, [16.0, 50.0, 84.0])
        summary[key] = {
            "median": float(q50),
            "p16": float(q16),
            "p84": float(q84),
            "minus": float(q50 - q16),
            "plus": float(q84 - q50),
            "n": int(values.size),
        }
    return summary


def analyze_pyneb_with_uncertainties(
    flux: Mapping[Wave, float],
    flux_error: Mapping[Wave, float],
    *,
    upper_limits: Mapping[Wave, float] | None = None,
    status_override: Mapping[Wave, LineStatus] | None = None,
    snr_threshold: float = 3.0,
    upper_limit_sigma: float = 3.0,
    n_mc: int = 10_000,
    seed: int = 12345,
    extinction_law: str = "CCM89",
    r_v: float = 3.1,
    clip_negative_ebv: bool = True,
    fallback_density: float = 100.0,
    max_iterations: int = 8,
    tie_fixed_doublets: bool = True,
    te_low_relation_scatter_k: float = 0.0,
    keep_samples: bool = False,
) -> dict[str, object]:
    """
    Run central and Monte-Carlo PyNeb analyses.

    Upper-limit behavior
    --------------------
    - Low-S/N lines are classified as upper limits.
    - Upper-limit lines are excluded from Gaussian Monte-Carlo draws.
    - [N II] limits produce N/O upper bounds when O+ is measured.
    - [O II] limits produce O+ and total-O/H upper bounds.
    - A 4363 limit produces Te upper and O/H lower bounds, not a point Te.

    The Monte-Carlo model assumes independent Gaussian errors for detected
    lines. If your spectral fitter provides correlated posterior samples,
    those should be propagated directly instead of using this approximation.
    """
    if n_mc < 0:
        raise ValueError("n_mc cannot be negative.")
    if fallback_density <= 0:
        raise ValueError("fallback_density must be positive.")
    if te_low_relation_scatter_k < 0:
        raise ValueError("te_low_relation_scatter_k cannot be negative.")

    measurements = classify_lines(
        flux,
        flux_error,
        upper_limits=upper_limits,
        status_override=status_override,
        snr_threshold=snr_threshold,
        upper_limit_sigma=upper_limit_sigma,
    )
    central_flux, doublet_models = _prepare_central_fluxes(
        measurements,
        tie_fixed_doublets=tie_fixed_doublets,
    )

    warnings: list[str] = []
    for wave, m in measurements.items():
        if m.status == "upper_limit":
            warnings.append(
                f"{wave} is treated as an upper limit"
                f" ({m.upper_limit:.4g}; S/N={m.snr:.3g} if available)."
            )

    # Hbeta is essential for all abundance normalizations.
    if measurements[4861].status != "detection":
        raise RuntimeError("Hbeta must be a detection for this analysis.")
    if measurements[6563].status != "detection":
        raise RuntimeError(
            "Halpha must be a detection for the current Balmer-decrement correction."
        )
    if 4959 not in central_flux or 5007 not in central_flux:
        raise RuntimeError("At least one nebular [O III] line must be detected.")

    # Central estimate or one-sided temperature solution.
    if measurements[4363].status == "detection":
        central = _solve_direct_realization(
            central_flux,
            extinction_law=extinction_law,
            r_v=r_v,
            clip_negative_ebv=clip_negative_ebv,
            fallback_density=fallback_density,
            max_iterations=max_iterations,
        )
        _add_one_sided_bounds(
            central,
            measurements,
            hbeta_flux=central_flux[4861],
            extinction_law=extinction_law,
            r_v=r_v,
        )
        mode = "semi-direct"
    else:
        upper_4363 = measurements[4363].upper_limit
        if upper_4363 is None:
            raise RuntimeError("4363 is not detected and has no upper limit.")
        central = _solve_temperature_upper_limit(
            central_flux,
            upper_4363,
            extinction_law=extinction_law,
            r_v=r_v,
            clip_negative_ebv=clip_negative_ebv,
            fallback_density=fallback_density,
            max_iterations=max_iterations,
        )
        mode = "temperature upper limit"

    central["PyNeb_version"] = pn.__version__
    central["analysis_mode"] = mode
    central["extinction_law"] = extinction_law
    central["R_V"] = r_v
    central["low_temperature_method"] = "Te(low)=0.7 Te(OIII)+3000 K"

    rng = np.random.default_rng(seed)
    samples: list[dict[str, object]] = []
    failures: Counter[str] = Counter()

    for _ in range(n_mc):
        draw = _draw_detected_fluxes(
            measurements,
            doublet_models,
            rng,
            tie_fixed_doublets=tie_fixed_doublets,
        )

        # Never redraw until positive: failed draws are counted. Rejection
        # preserves the original Gaussian measurement model.
        required_positive = [4861, 6563, 4959, 5007]
        if measurements[4363].status == "detection":
            required_positive.append(4363)
        if any(w not in draw or draw[w] <= 0 or not np.isfinite(draw[w]) for w in required_positive):
            failures["nonpositive_required_flux"] += 1
            continue

        relation_offset = (
            float(rng.normal(0.0, te_low_relation_scatter_k))
            if te_low_relation_scatter_k > 0
            else 0.0
        )

        try:
            if measurements[4363].status == "detection":
                sample = _solve_direct_realization(
                    draw,
                    extinction_law=extinction_law,
                    r_v=r_v,
                    clip_negative_ebv=clip_negative_ebv,
                    fallback_density=fallback_density,
                    max_iterations=max_iterations,
                    te_low_scatter_offset=relation_offset,
                )
                _add_one_sided_bounds(
                    sample,
                    measurements,
                    hbeta_flux=draw[4861],
                    extinction_law=extinction_law,
                    r_v=r_v,
                )
            else:
                upper_4363 = measurements[4363].upper_limit
                assert upper_4363 is not None
                sample = _solve_temperature_upper_limit(
                    draw,
                    upper_4363,
                    extinction_law=extinction_law,
                    r_v=r_v,
                    clip_negative_ebv=clip_negative_ebv,
                    fallback_density=fallback_density,
                    max_iterations=max_iterations,
                )

            samples.append(sample)
        except RuntimeError:
            failures["diagnostic_inversion"] += 1
        except ValueError:
            failures["invalid_realization"] += 1
        except (FloatingPointError, OverflowError, ZeroDivisionError):
            failures["numerical_failure"] += 1

    mc = {
        "n_requested": int(n_mc),
        "n_valid": int(len(samples)),
        "valid_fraction": float(len(samples) / n_mc) if n_mc else math.nan,
        "failures": dict(failures),
        "summary": _summarize_samples(samples),
    }
    if keep_samples:
        mc["samples"] = samples

    line_report = {
        wave: {
            "flux": m.flux,
            "error": m.error,
            "snr": m.snr,
            "status": m.status,
            "upper_limit": m.upper_limit,
        }
        for wave, m in measurements.items()
    }

    return {
        "central": central,
        "line_report": line_report,
        "monte_carlo": mc,
        "warnings": tuple(warnings),
    }


def _format_summary(summary: Mapping[str, Mapping[str, float]], key: str, fmt: str) -> str | None:
    if key not in summary:
        return None
    row = summary[key]
    return (
        f"{format(row['median'], fmt)} "
        f"(-{format(row['minus'], fmt)}, +{format(row['plus'], fmt)})"
    )


def print_pyneb_uncertainty_result(result: Mapping[str, object]) -> None:
    central = result["central"]
    mc = result["monte_carlo"]
    summary = mc["summary"]

    print(f"PyNeb version           = {central['PyNeb_version']}")
    print(f"Analysis mode           = {central['analysis_mode']}")
    print(
        f"Extinction law          = {central['extinction_law']} "
        f"(R_V={central['R_V']})"
    )
    print(f"Temperature method      = {central['temperature_method']}")
    print(f"Density method          = {central['density_method']}")

    for key, label, fmt in (
        ("Te_OIII_K", "Te([O III]) [K]", ".0f"),
        ("Te_OIII_upper_K", "Te([O III]) upper [K]", ".0f"),
        ("Te_low_K", "Te(low) [K]", ".0f"),
        ("ne_cm-3", "ne [cm^-3]", ".2g"),
        ("12_log_O_H", "12+log(O/H)", ".3f"),
        ("12_log_O_H_lower", "12+log(O/H) lower", ".3f"),
        ("12_log_O_H_upper", "12+log(O/H) upper", ".3f"),
        ("log_N_O", "log(N/O)", ".3f"),
        ("log_N_O_upper", "log(N/O) upper", ".3f"),
    ):
        text = _format_summary(summary, key, fmt)
        if text is not None:
            print(f"{label:24s}= {text}")
        elif key in central:
            print(f"{label:24s}= {format(float(central[key]), fmt)}")

    print(
        f"Valid MC realizations   = {mc['n_valid']}/{mc['n_requested']} "
        f"({100.0 * mc['valid_fraction']:.1f}%)"
    )
    if mc["failures"]:
        print(f"MC failures             = {mc['failures']}")

    line_report = result["line_report"]
    limits = [
        (wave, row)
        for wave, row in line_report.items()
        if row["status"] == "upper_limit"
    ]
    if limits:
        print("Upper-limit lines:")
        for wave, row in limits:
            print(
                f"  {wave}: flux<{row['upper_limit']:.4g} "
                f"(measured S/N={row['snr']:.3g})"
            )

    warnings = result["warnings"]
    if warnings:
        print("Warnings:")
        for warning in warnings:
            print(f"  - {warning}")


if __name__ == "__main__":
    # Demonstration values only. Replace errors and limits with those from
    # your spectral fit / variance spectrum.
    flux = {
        3727: 3,
        4363: 2.8,
        4861: 19.4,
        4959: 36.3,
        5007: 98.2,
        6548: 1.4,
        6563: 52.5,
        6584: 1.,
        6716: 1.,
        6731: 1.,
    }
    
    flux_error = {
        3727: 0.6,
        4363: 0.5,
        4861: 0.5,
        4959: 0.6,
        5007: 1.5,
        6548: 1.4,
        6563: 0.9,
        6584: 1.4,
        6716: 1.6,
        6731: 1.7,
    }

    result = analyze_pyneb_with_uncertainties(
        flux,
        flux_error,
        n_mc=100,
        extinction_law="G03 LMC",
        r_v=2.77,
        # Manual overrides are useful when visual inspection shows that a
        # nominally positive fitted feature is not trustworthy.
        status_override={6548: "upper_limit", 6584: "upper_limit", 6731: "upper_limit"},
    )
    print_pyneb_uncertainty_result(result)



#==================================================================
#========================  END of Functions  ======================
#==================================================================
def flux_from_df(index, df, source='MAST', printWaves=False):

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




