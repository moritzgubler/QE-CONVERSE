#!/usr/bin/env python3
"""
parse_nmr.py — Parse QE-CONVERSE output files, assemble NMR shielding tensors,
compute quadrupolar coupling constants, and plot predicted NMR spectra.

Usage:
    parse_nmr.py run_x.out run_y.out run_z.out --json params.json --B0 9.4 --output spectrum.pdf

The three output files must correspond to converse runs with m_0(1)=1, m_0(2)=1, m_0(3)=1
(one per Cartesian direction of the magnetic dipole).

JSON format (element-specific parameters):
    {
        "Si": { "spin": 0.5,  "ref": 79494000.0,  "lw": 200.0  },
        "O":  { "spin": 2.5,  "ref": 54246000.0,  "lw": 500.0  }
    }

    spin : nuclear spin quantum number I
    ref  : absolute resonance frequency (Hz) of the reference compound at field B0
           (e.g. the TMS peak frequency measured on your spectrometer for 29Si)
    lw   : optional Gaussian linewidth (Hz); default 200 Hz
"""

import argparse
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


DEFAULT_LW = 200.0   # Hz

# Cq (MHz) = CQ_FACTOR * Vzz (Ha/bohr²) * Q (barn)
# Derived from: Cq = e * Q * Vzz_SI / h
#   1 a.u. EFG = 9.7173624e21 V/m², 1 barn = 1e-28 m²
CQ_FACTOR = (1.602176634e-19 * 9.7173624e21 * 1e-28
             / 6.62607015e-34 / 1e6)   # ≈ 234.96 MHz/(Ha/bohr² · barn)


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def _parse_direction(path):
    """
    Parse one QE-CONVERSE output file.

    Returns a dict with:
        dipole_dir : int   1, 2 or 3 (direction of m_0)
        shifts     : dict  {atom_index: {'element': str, 'shift': np.ndarray(3),
                                         'core': float}}
        efg        : dict  {atom_index: {'element': str, 'Vzz': float,
                                         'eta': float, 'Cq': float or None}}
    """
    text = Path(path).read_text()

    # --- dipole direction from output header ---
    # Format: "m_0 = (  1.00  0.00  0.00)  mu_b"
    dipole_dir = None
    m = re.search(r'm_0\s*=\s*\(\s*([\d.eE+\-]+)\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)\s*\)', text)
    if m:
        vec = [float(m.group(i)) for i in range(1, 4)]
        for i, v in enumerate(vec):
            if abs(v) > 1e-6:
                dipole_dir = i + 1
                break
    if dipole_dir is None:
        # fallback: lambda_so for EPR runs
        m = re.search(r'lambda_so\s*=\s*\(\s*([\d.eE+\-]+)\s+([\d.eE+\-]+)\s+([\d.eE+\-]+)\s*\)', text)
        if m:
            vec = [float(m.group(i)) for i in range(1, 4)]
            for i, v in enumerate(vec):
                if abs(v) > 1e-6:
                    dipole_dir = i + 1
                    break
    if dipole_dir is None:
        raise ValueError(f"{path}: cannot determine dipole direction (m_0 not set?)")

    # --- chemical shifts ---
    # Output format (one block per target atom per run):
    #   NUCLEAR DIPOLE ON ATOM   1 (Si ):
    #   m_0                =       1.414214      0.000000      0.000000
    #   Chemical shift (ppm):     -409.5718        0.7873      -10.6118
    #   Core shift     (ppm):      836.9099

    shifts = {}
    core_shifts = {}

    num = r'([\d.\-eE+]+)'
    block_re = re.compile(
        r'NUCLEAR DIPOLE ON ATOM\s+(\d+)\s+\((\w+)\s*\).*?'
        r'Chemical shift \(ppm\):\s*' + num + r'\s+' + num + r'\s+' + num,
        re.DOTALL
    )
    core_re = re.compile(r'Core shift\s+\(ppm\):\s*' + num)

    for m in block_re.finditer(text):
        aidx = int(m.group(1))
        elem = m.group(2).strip()
        vec  = np.array([float(m.group(3)), float(m.group(4)), float(m.group(5))])
        shifts[aidx] = {'element': elem, 'shift': vec}

    for m in core_re.finditer(text):
        # core is scalar (isotropic); store under the last found atom or all atoms
        pass   # filled below together with the block

    # Re-parse blocks with core included
    full_block_re = re.compile(
        r'NUCLEAR DIPOLE ON ATOM\s+(\d+)\s+\((\w+)\s*\).*?'
        r'Chemical shift \(ppm\):\s*' + num + r'\s+' + num + r'\s+' + num +
        r'.*?Core shift\s+\(ppm\):\s*' + num,
        re.DOTALL
    )
    for m in full_block_re.finditer(text):
        aidx = int(m.group(1))
        core_shifts[aidx] = float(m.group(6))

    # --- EFG ---
    # Format: "  Si      1  Vzz=    0.0871    eta= 0.67021"
    # Optional: "  Cq=  X.XX MHz" appended to same line
    efg = {}
    efg_re = re.compile(
        r'^\s+(\w+)\s+(\d+)\s+Vzz=\s*' + num + r'\s+eta=\s*' + num +
        r'(?:\s+Cq=\s*' + num + r'\s*MHz)?',
        re.MULTILINE
    )
    for m in efg_re.finditer(text):
        elem = m.group(1)
        aidx = int(m.group(2))
        Vzz  = float(m.group(3))
        eta  = float(m.group(4))
        Cq   = float(m.group(5)) if m.group(5) else None
        efg[aidx] = {'element': elem, 'Vzz': Vzz, 'eta': eta, 'Cq': Cq}

    return {
        'dipole_dir': dipole_dir,
        'shifts': shifts,
        'core_shifts': core_shifts,
        'efg': efg,
        'path': str(path),
    }


def parse_files(paths):
    """
    Parse a list of output files and return per-atom assembled data.

    Returns:
        atoms : dict  {atom_index: {
                    'element': str,
                    'sigma': np.ndarray(3,3) or None,   # full shielding tensor
                    'core':  float or None,
                    'efg':   dict or None,
               }}
        missing_dirs : set of atom indices that lack all 3 directions
    """
    # runs[(atom_idx, dipole_dir)] = parsed result
    runs = {}
    efg_store = {}

    for path in paths:
        r = _parse_direction(path)
        d = r['dipole_dir']
        # The target atom is the one whose shift was printed in this run
        if not r['shifts']:
            print(f"WARNING: {path}: no chemical shift found, skipping.")
            continue
        for aidx in r['shifts']:
            key = (aidx, d)
            if key in runs:
                print(f"WARNING: duplicate (atom={aidx}, dir={d}), using {path}")
            runs[key] = r
        efg_store.update(r['efg'])

    # Collect all atom indices seen
    all_atoms = set(aidx for (aidx, _) in runs)

    atoms = {}
    missing_dirs = set()

    for aidx in sorted(all_atoms):
        elem = None
        for (a, d), r in runs.items():
            if a == aidx and aidx in r['shifts']:
                elem = r['shifts'][aidx]['element']
                break

        # Assemble 3x3 tensor: column beta = shift vector from beta-run
        has_all = all((aidx, d) in runs for d in (1, 2, 3))
        if has_all:
            sigma = np.zeros((3, 3))
            core = None
            for d in (1, 2, 3):
                r = runs[(aidx, d)]
                sigma[:, d - 1] = r['shifts'][aidx]['shift']
                if core is None and aidx in r.get('core_shifts', {}):
                    core = r['core_shifts'][aidx]
        else:
            sigma = None
            missing_dirs.add(aidx)
            core = None

        atoms[aidx] = {
            'element': elem,
            'sigma': sigma,
            'core': core,
            'efg': efg_store.get(aidx),
        }

    return atoms, missing_dirs


# ---------------------------------------------------------------------------
# Tensor analysis
# ---------------------------------------------------------------------------

def principal_values_haeberlen(sigma):
    """
    Symmetrize sigma, diagonalize, return principal values in Haeberlen convention.

    Returns:
        sigma_iso   : isotropic shift (ppm)
        delta_aniso : anisotropy = sigma_zz - sigma_iso  (Haeberlen delta)
        eta_cs      : asymmetry = (sigma_yy - sigma_xx) / delta_aniso
        eigenvalues : sorted eigenvalues [sigma_xx, sigma_yy, sigma_zz] (Haeberlen)
    """
    sym = 0.5 * (sigma + sigma.T)
    evals = np.linalg.eigvalsh(sym)
    sigma_iso = np.mean(evals)

    # Sort by |sigma_ii - sigma_iso|, largest last (= sigma_zz in Haeberlen)
    order = np.argsort(np.abs(evals - sigma_iso))
    evals = evals[order]   # [sigma_xx, sigma_yy, sigma_zz]

    delta_aniso = evals[2] - sigma_iso
    if abs(delta_aniso) > 1e-10:
        eta_cs = (evals[1] - evals[0]) / delta_aniso
    else:
        eta_cs = 0.0

    return sigma_iso, delta_aniso, eta_cs, evals


# ---------------------------------------------------------------------------
# Spectrum simulation
# ---------------------------------------------------------------------------

def _compute_cq(efg, elem_params):
    """Return Cq in MHz from parsed EFG and element params, or None."""
    if efg is None or efg.get('Vzz') is None:
        return None
    Q = elem_params.get('Q')
    if Q is None:
        return None
    return CQ_FACTOR * efg['Vzz'] * Q


def _larmor(elem_params, B0):
    """Larmor frequency in Hz for element at field B0 (Tesla)."""
    gamma = elem_params.get('gamma')
    if gamma is None:
        raise KeyError("No 'gamma' for this element in the JSON.")
    return abs(gamma) * B0 / (2.0 * np.pi)


def _sigma_to_freq(sigma_ppm, nu_L, nu_ref):
    """
    Convert absolute shielding (ppm) to frequency (Hz) relative to reference.

        nu = nu_L * (1 - sigma*1e-6)   [absolute]
        nu - nu_ref                    [relative to reference compound]
    """
    return nu_L * (1.0 - sigma_ppm * 1e-6) - nu_ref


def _csa_powder_static(sigma_iso, delta_aniso, eta_cs, nu_L, nu_ref,
                       freqs, lw):
    """
    Static CSA powder pattern via (theta, phi) grid.
    Returns intensity array on `freqs` grid.
    """
    ntheta, nphi = 200, 400
    theta = np.linspace(0, np.pi, ntheta)
    phi   = np.linspace(0, 2 * np.pi, nphi)
    TH, PH = np.meshgrid(theta, phi, indexing='ij')
    sin_th = np.sin(TH)
    cos_th = np.cos(TH)

    # Orientation-dependent shielding (Mehring convention)
    sigma_orient = (sigma_iso
                    + delta_aniso * (3 * cos_th**2 - 1 - eta_cs * sin_th**2 * np.cos(2 * PH)) / 2.0)

    nu_orient = _sigma_to_freq(sigma_orient, nu_L, nu_ref)   # Hz, shape (ntheta, nphi)
    weights   = sin_th   # solid-angle weight

    intensity = np.zeros_like(freqs)
    sigma_g   = lw / (2.0 * np.sqrt(2.0 * np.log(2.0)))   # FWHM -> sigma
    for i in range(ntheta):
        for j in range(nphi):
            w = weights[i, j]
            intensity += w * np.exp(-0.5 * ((freqs - nu_orient[i, j]) / sigma_g)**2)

    intensity /= intensity.max() if intensity.max() > 0 else 1.0
    return intensity


def _quadrupolar_satellites_static(Cq_hz, eta_q, spin, nu_L, nu_ref,
                                    sigma_iso, freqs, lw):
    """
    First-order quadrupolar satellites for a spin-I nucleus (static powder).

    For each pair of satellite transitions (m <-> m-1, m != 1/2 <-> -1/2),
    add a powder-broadened line.

    Returns intensity array.
    """
    if spin <= 0.5 or Cq_hz == 0:
        return np.zeros_like(freqs)

    I = spin
    nu_Q = 3.0 * Cq_hz / (2.0 * I * (2.0 * I - 1.0))   # quadrupolar frequency

    intensity = np.zeros_like(freqs)
    sigma_g   = lw / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    # Satellite transitions: m <-> m-1 for m = -I+1 .. I, excluding central (1/2 <-> -1/2)
    m_vals = np.arange(-I + 1, I + 0.5, 1.0)
    for m in m_vals:
        if abs(m - 0.5) < 1e-9:
            continue   # skip central transition

        ntheta, nphi = 150, 300
        theta = np.linspace(0, np.pi, ntheta)
        phi   = np.linspace(0, 2 * np.pi, nphi)
        TH, PH = np.meshgrid(theta, phi, indexing='ij')
        sin_th = np.sin(TH)
        cos_th = np.cos(TH)

        # First-order quadrupolar shift for transition m <-> m-1
        nu_sat = (_sigma_to_freq(sigma_iso, nu_L, nu_ref)
                  + nu_Q * (m - 0.5) * (3 * cos_th**2 - 1 + eta_q * sin_th**2 * np.cos(2 * PH)))

        weights = sin_th
        for i in range(ntheta):
            for j in range(nphi):
                intensity += weights[i, j] * np.exp(
                    -0.5 * ((freqs - nu_sat[i, j]) / sigma_g)**2)

    if intensity.max() > 0:
        intensity /= intensity.max()
    return intensity


def simulate_atom_spectrum(atom_data, elem_params, B0, freqs):
    """
    Simulate the NMR spectrum for one atom.

    Returns (central, satellites) intensity arrays on `freqs`.
    """
    sigma   = atom_data['sigma']
    efg     = atom_data['efg']
    spin    = elem_params['spin']
    nu_ref  = elem_params['ref']
    lw      = elem_params.get('lw', DEFAULT_LW)
    nu_L    = _larmor(elem_params, B0)

    sigma_iso, delta_aniso, eta_cs, _ = principal_values_haeberlen(sigma)

    central = _csa_powder_static(sigma_iso, delta_aniso, eta_cs,
                                  nu_L, nu_ref, freqs, lw)

    satellites = np.zeros_like(freqs)
    Cq_mhz = _compute_cq(efg, elem_params)
    if Cq_mhz is not None and spin > 0.5:
        Cq_hz = Cq_mhz * 1e6   # MHz -> Hz
        eta_q = efg['eta'] if efg is not None and efg['eta'] is not None else 0.0
        satellites = _quadrupolar_satellites_static(
            Cq_hz, eta_q, spin, nu_L, nu_ref, sigma_iso, freqs, lw)

    return central, satellites


# ---------------------------------------------------------------------------
# Frequency axis helpers
# ---------------------------------------------------------------------------

def _auto_range(peak_freqs, atoms, params, margin_factor=3.0):
    """
    Auto-determine frequency range covering all peaks and their quadrupolar satellites.
    """
    if not peak_freqs:
        return -1e4, 1e4

    extras = []
    for _, adat in atoms.items():
        elem = adat['element']
        if elem not in params or adat['sigma'] is None:
            continue
        efg = adat['efg']
        if efg is not None and efg.get('Cq') is not None:
            extras.append(efg['Cq'] * 1e6)   # rough satellite width in Hz

    half_width = max(max(extras) * 2.0 if extras else 0,
                     abs(max(peak_freqs) - min(peak_freqs)) * 0.5) + 2000.0

    centre = 0.5 * (max(peak_freqs) + min(peak_freqs))
    margin = half_width * margin_factor
    return centre - margin, centre + margin


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def plot_spectrum(atoms, params, B0, output_pdf):
    # Group atoms by element
    by_elem = defaultdict(list)
    for aidx, adat in atoms.items():
        elem = adat['element']
        if elem in params and adat['sigma'] is not None:
            by_elem[elem].append((aidx, adat))

    if not by_elem:
        print("ERROR: no atoms with complete tensor and matching JSON params.")
        sys.exit(1)

    n_panels = len(by_elem)
    fig, axes = plt.subplots(n_panels, 1, figsize=(10, 4 * n_panels), squeeze=False)

    for panel_idx, (elem, atom_list) in enumerate(sorted(by_elem.items())):
        ax = axes[panel_idx, 0]
        ep = params[elem]
        nu_L   = _larmor(ep, B0)
        nu_ref = ep['ref']

        # Build frequency axis
        peak_freqs = [_sigma_to_freq(
                          principal_values_haeberlen(adat['sigma'])[0],
                          nu_L, nu_ref)
                      for _, adat in atom_list]
        fmin, fmax = _auto_range(peak_freqs, dict(atom_list), params)
        lw = ep.get('lw', DEFAULT_LW)
        npts = max(4096, int((fmax - fmin) / (lw / 10.0)))
        freqs = np.linspace(fmin, fmax, npts)

        total = np.zeros_like(freqs)
        for aidx, adat in atom_list:
            central, sats = simulate_atom_spectrum(adat, ep, B0, freqs)
            total += central + 0.5 * sats   # satellites at half weight

        total /= total.max() if total.max() > 0 else 1.0

        ax.plot(freqs / 1e3, total, color='navy', lw=1.2)
        ax.set_xlabel("Frequency offset from reference (kHz)", fontsize=11)
        ax.set_ylabel("Intensity (arb.)", fontsize=11)
        ax.set_title(f"Element: {elem}  |  B₀ = {B0:.2f} T  |  "
                     f"ν_L = {nu_L/1e6:.3f} MHz", fontsize=11)

        # Annotate each atom with its isotropic frequency
        for aidx, adat in atom_list:
            sigma_iso, _, _, _ = principal_values_haeberlen(adat['sigma'])
            nu_iso = _sigma_to_freq(sigma_iso, nu_L, nu_ref)
            efg = adat['efg']
            label = f"atom {aidx}\nσ_iso={sigma_iso:.1f} ppm"
            if efg and efg.get('Cq') is not None:
                label += f"\nCq={efg['Cq']:.2f} MHz\nη={efg['eta']:.2f}"
            ymax = total.max()
            ax.annotate(label,
                        xy=(nu_iso / 1e3, ymax * 0.95),
                        xytext=(nu_iso / 1e3, ymax * 1.05),
                        fontsize=7, ha='center',
                        arrowprops=dict(arrowstyle='->', lw=0.8))

        ax.set_xlim(fmin / 1e3, fmax / 1e3)
        ax.set_ylim(bottom=0)
        ax.tick_params(labelsize=10)

    fig.tight_layout()
    fig.savefig(output_pdf, format="pdf", bbox_inches="tight")
    print(f"Spectrum written to {output_pdf}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Parse QE-CONVERSE output, assemble NMR tensors, plot spectrum.")
    parser.add_argument("files", nargs="+",
                        help="QE-CONVERSE output files (one per dipole direction x/y/z)")
    parser.add_argument("--json", required=True, metavar="PARAMS",
                        help=(
                            "JSON file with element-specific parameters. "
                            "Keys starting with '_' are ignored (use for comments). "
                            "Required per element: "
                            "gamma (gyromagnetic ratio, rad/s/T), "
                            "spin (nuclear spin I), "
                            "ref (resonance frequency in Hz of the reference compound "
                            "at the given B0, measured on your spectrometer). "
                            "Optional: "
                            "Q (nuclear quadrupole moment in barn, Stone 2016; "
                            "omit for spin-1/2 nuclei), "
                            "lw (Gaussian linewidth in Hz, default 200). "
                            "Example: "
                            '{\"O\": {\"gamma\": -36.281e6, \"spin\": 2.5, '
                            '\"Q\": -0.02558, \"ref\": 54234000.0, \"lw\": 500.0}}'
                        ))
    parser.add_argument("--B0", type=float, required=True, metavar="TESLA",
                        help="Static magnetic field strength in Tesla")
    parser.add_argument("--output", required=True, metavar="FILE.pdf",
                        help="Output PDF filename")
    args = parser.parse_args()

    with open(args.json) as fh:
        params = {k: v for k, v in json.load(fh).items()
                  if not k.startswith('_')}

    atoms, missing = parse_files(args.files)

    if missing:
        print(f"WARNING: atoms {sorted(missing)} are missing one or more dipole "
              f"directions — shielding tensor incomplete, these atoms will be skipped.")

    # Print summary table
    print(f"\n{'Atom':>6}  {'Elem':>4}  {'σ_iso (ppm)':>12}  "
          f"{'δ_aniso':>10}  {'η_cs':>6}  {'Cq (MHz)':>10}  {'η_q':>6}")
    print("-" * 70)
    for aidx, adat in sorted(atoms.items()):
        if adat['sigma'] is None:
            continue
        sig_iso, d_aniso, eta_cs, _ = principal_values_haeberlen(adat['sigma'])
        efg = adat['efg']
        Cq_val  = _compute_cq(efg, params.get(adat['element'], {}))
        Cq_str  = f"{Cq_val:10.3f}" if Cq_val is not None else f"{'—':>10}"
        eta_str = f"{efg['eta']:6.3f}" if efg and efg.get('eta') is not None else f"{'—':>6}"
        print(f"{aidx:>6}  {adat['element']:>4}  {sig_iso:12.2f}  "
              f"{d_aniso:10.2f}  {eta_cs:6.3f}  {Cq_str}  {eta_str}")
    print()

    plot_spectrum(atoms, params, args.B0, args.output)


if __name__ == "__main__":
    main()
