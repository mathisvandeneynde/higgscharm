"""Event weight for the data-driven jet -> l (fake-lepton) estimate in e-mu.

For every event in a fake category (e.g. emu_qcd = the AR selection), the weight is

    w = s * F_T(dR, N_jet) * C_e(pT_e) * C_mu(pT_mu) * b

with s = +1 for data and s = -1 for MC (the MC events also keep their usual weights),
so that data - MC (genuine leptons) is extrapolated from same sign to opposite sign.

F_T, C_e, C_mu and b come from the correctionlib file made by scripts/qcd/make_fake_cset.py.
dR, N_jet and the lepton pTs are computed from exactly the same collections as the fake_*
histograms of the workflow (via pair_key / jets_key), so the factors are applied in the cells
in which they were measured. Values outside a table are clamped, like the histogram filling.
"""
import numpy as np
import awkward as ak
import correctionlib
from analysis.corrections.correctionlib_files import correction_files
from analysis.selections.emu_regions import pair_key, jets_key


def add_fake_factor_weights(events, year, category, bias_setting, weights_container, shift, is_mc):
    """Add the fake-factor weight (and its variations) to the weights container.

    bias_setting: 'measured' -> b from the correction file; a number -> fixed b (e.g. 1.0)
    """
    cset = correctionlib.CorrectionSet.from_file(correction_files["fake_factors"][year])

    # the selected pair (exactly one per event in these categories) and the region's jets
    pair = ak.firsts(events[f"selected_{pair_key(category)}"])
    jets = events[f"selected_{jets_key(category)}"]
    l1, l2 = pair.l1, pair.l2
    dr = ak.to_numpy(l1.delta_r(l2))
    njet = np.minimum(ak.to_numpy(ak.num(jets, axis=1)), 2)
    pt_e = ak.to_numpy(ak.where(np.abs(l1.pdgId) == 11, l1.pt, l2.pt))
    pt_mu = ak.to_numpy(ak.where(np.abs(l1.pdgId) == 13, l1.pt, l2.pt))

    ft = cset["FT"].evaluate("nominal", dr, njet)
    c_e = cset["C_e"].evaluate("nominal", pt_e)
    c_mu = cset["C_mu"].evaluate("nominal", pt_mu)
    use_measured_bias = str(bias_setting).lower() == "measured"
    if use_measured_bias:
        b = cset["bias"].evaluate("nominal", njet)
    else:
        b = np.full(len(dr), float(bias_setting))

    sign = -1.0 if is_mc else 1.0
    weights_container.add("fake_factor", sign * ft * c_e * c_mu * b)

    # variations: only for the nominal shift (like all other weight variations)
    if shift is None:
        ones = np.ones(len(dr))
        for j in range(3):
            up = cset["FT"].evaluate(f"up_nj{j}", dr, njet) / ft
            down = cset["FT"].evaluate(f"down_nj{j}", dr, njet) / ft
            weights_container.add(f"CMS_fake_FT_nj{j}_{year}", ones, weightUp=up, weightDown=down)
        for name, pt, label in (("C_e", pt_e, "e"), ("C_mu", pt_mu, "m")):
            nom = cset[name].evaluate("nominal", pt)
            weights_container.add(
                f"CMS_fake_closure_{label}_{year}",
                ones,
                weightUp=cset[name].evaluate("up", pt) / nom,
                weightDown=cset[name].evaluate("down", pt) / nom,
            )
        if use_measured_bias:
            weights_container.add(
                f"CMS_fake_bias_{year}",
                ones,
                weightUp=cset["bias"].evaluate("up", njet) / b,
                weightDown=cset["bias"].evaluate("down", njet) / b,
            )
