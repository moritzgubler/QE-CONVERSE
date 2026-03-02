"""
Parse QE-CONVERSE text output.

All public functions operate on a string (stdout/file contents) and return
plain Python dicts with float values or lists of floats.
"""
import re


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _floats(text):
    """Return all floats found in *text*, handling -0.0000 and D-exponents."""
    # Fortran sometimes writes 1.23D-04 — normalise to E
    text = text.replace('D', 'E').replace('d', 'e')
    return [float(x) for x in re.findall(r'[-+]?\d+\.\d+(?:[eE][+-]?\d+)?', text)]


def _first_float(pattern, text):
    """Return the first float on the line matched by *pattern*, or None."""
    m = re.search(pattern, text)
    if m is None:
        return None
    vals = _floats(m.group(1))
    return vals[0] if vals else None


def _first_vec3(pattern, text):
    """Return the first 3-vector matched by *pattern*, or None."""
    m = re.search(pattern, text)
    if m is None:
        return None
    vals = _floats(m.group(1))
    return vals[:3] if len(vals) >= 3 else None


def _last_vec3(pattern, text):
    """Return the last 3-vector matched by *pattern* (useful for duplicate keys)."""
    matches = re.findall(pattern, text)
    if not matches:
        return None
    vals = _floats(matches[-1])
    return vals[:3] if len(vals) >= 3 else None


# ---------------------------------------------------------------------------
# EPR g-tensor quantities
# ---------------------------------------------------------------------------

def parse_gtensor(text):
    """
    Parse EPR g-tensor results from QE-CONVERSE output text.

    Returns a dict with the following keys (all values are floats or
    lists of 3 floats):

        delta_g_rmc         – scalar RMC correction (ppm)
        delta_g_rmc_gipaw   – scalar GIPAW correction to RMC (ppm)
        delta_g_so          – [x, y, z] SO contribution (ppm)
        delta_g_total       – [x, y, z] total Δg (ppm)
        m_lc                – [x, y, z] LC magnetization (with Berry curvature)
        m_ic                – [x, y, z] IC magnetization (with Berry curvature)
        m_tot               – [x, y, z] total magnetization (with Berry curvature)
        lambda_so           – [x, y, z] spin-orbit direction
    """
    result = {}

    # --- Scalar g-tensor corrections ---
    # "delta_g RMC\s+=" requires whitespace between RMC and =, which
    # distinguishes it from "delta_g RMC(GIPAW) =" (no whitespace after RMC).
    result['delta_g_rmc'] = _first_float(
        r'delta_g RMC\s+=\s+(.*)', text)
    result['delta_g_rmc_gipaw'] = _first_float(
        r'delta_g RMC\(GIPAW\)\s+=\s+(.*)', text)

    # --- Vector g-tensor components ---
    result['delta_g_so'] = _first_vec3(
        r'delta_g SO\s+=\s+(.*)', text)
    result['delta_g_total'] = _first_vec3(
        r'delta_g total\s+=\s+(.*)', text)

    # --- Magnetization (with Berry curvature block, last occurrence) ---
    # The output prints M_LC etc. twice: without and with Berry curvature.
    # We take the last match, which corresponds to the "with Berry curvature" block.
    result['m_lc'] = _last_vec3(r'M_LC\s+=\s+(.*)', text)
    result['m_ic'] = _last_vec3(r'M_IC\s+=\s+(.*)', text)
    result['m_tot'] = _last_vec3(r'M_tot\s+=\s+(.*)', text)

    # --- Run settings ---
    result['lambda_so'] = _first_vec3(r'lambda_so\s+=\s+(.*)', text)

    return result


# ---------------------------------------------------------------------------
# NMR chemical shift quantities
# ---------------------------------------------------------------------------

def parse_chemical_shift(text, atom_index=None):
    """
    Parse NMR chemical shift tensors from QE-CONVERSE output text.

    If *atom_index* is given (1-based), return only that atom's tensor.
    Otherwise return a list of 3×3 tensors (one per atom per direction).

    Currently returns:
        shifts  – list of dicts, each with keys:
                    'sigma'  – 3×3 tensor as list-of-lists (ppm)
                    'isotropic' – Tr(sigma)/3 (ppm)
    """
    # Locate all "Chemical shift (ppm)" blocks
    # Format (example):
    #   Chemical shift (ppm):
    #      -107.53  ...
    shifts = []
    blocks = re.findall(
        r'Chemical shift\s*\(ppm\):\s*\n((?:\s+.*\n){3})', text)
    for block in blocks:
        rows = []
        for line in block.strip().splitlines():
            vals = _floats(line)
            if vals:
                rows.append(vals[:3])
        if len(rows) == 3:
            iso = sum(rows[i][i] for i in range(3)) / 3.0
            shifts.append({'sigma': rows, 'isotropic': iso})

    if atom_index is not None:
        idx = atom_index - 1
        return shifts[idx] if idx < len(shifts) else None
    return shifts


# ---------------------------------------------------------------------------
# NMR chemical shift (single-direction run)
# ---------------------------------------------------------------------------

def parse_nmr_shift(text):
    """
    Parse NMR chemical shift from a single-direction QE-CONVERSE output.

    Each NMR run sets m_0(dir) = 1.0 for one atom and one Cartesian direction.
    The output contains a single "Chemical shift" line with three components.

    Returns a dict with:
        shift   – [s1, s2, s3] chemical shift components (ppm), or None
        core    – core shift (ppm), or None
        atom    – atom index (1-based int), or None
        element – element symbol, or None
    """
    import re as _re
    result = {}

    # "Chemical shift (ppm):    -404.2546    0.5417    -9.9074"
    m = _re.search(r'Chemical shift\s*\(ppm\):\s*(.*)', text)
    if m:
        vals = _floats(m.group(1))
        result['shift'] = vals[:3] if len(vals) >= 3 else None
    else:
        result['shift'] = None

    # "Core shift     (ppm):    837.9127"
    m = _re.search(r'Core shift\s*\(ppm\):\s*(.*)', text)
    if m:
        vals = _floats(m.group(1))
        result['core'] = vals[0] if vals else None
    else:
        result['core'] = None

    # "NUCLEAR DIPOLE ON ATOM   1 (Si ):"
    m = _re.search(r'NUCLEAR DIPOLE ON ATOM\s+(\d+)\s+\((\w+)', text)
    if m:
        result['atom'] = int(m.group(1))
        result['element'] = m.group(2).strip()
    else:
        result['atom'] = None
        result['element'] = None

    return result


# ---------------------------------------------------------------------------
# Generic run metadata
# ---------------------------------------------------------------------------

def parse_metadata(text):
    """
    Return a dict with basic run metadata:
        converged  – True if "JOB DONE" appears
        n_scf_iter – number of SCF iterations
        cpu_time   – total CPU time string
        wall_time  – total wall time string
    """
    result = {}
    result['converged'] = 'JOB DONE' in text

    m = re.search(r'convergence has been achieved in\s+(\d+)\s+iterations', text)
    result['n_scf_iter'] = int(m.group(1)) if m else None

    m = re.search(r'QE-CONVERSE\s*:\s*([\d.]+\w+)\s+CPU\s+([\d.]+\w+)\s+WALL', text)
    if m:
        result['cpu_time'] = m.group(1)
        result['wall_time'] = m.group(2)

    return result
