"""Mapping from category to the lepton pair and the jet collection that belong to it (xtoemu).

Used in the histogram expressions of xtoemu.yaml (via pair_key / jets_key) and later by the
fake-factor correction, so that both use exactly the same objects per region.
An unknown category raises an error instead of silently falling back to 'dileptons'.
"""

PAIRS = {
    # isolated e and mu (SR and control regions)
    "emu": "dileptons",
    "emu_nodzeta": "dileptons",
    "emu_nobveto": "dileptons",
    "emu_ttbar_cr": "dileptons",
    "emu_diboson_cr": "dileptons",
    # application region (AR): isolated, same sign
    "emu_ss": "dileptons_ss",
    # determination region (DR): anti-isolated muon, isolated electron
    "emu_dr_os": "dileptons_ai",
    "emu_dr_ss": "dileptons_ai_ss",
    # bias region B_iso: anti-isolated electron, isolated muon
    "emu_bi_os": "dileptons_bi",
    "emu_bi_ss": "dileptons_bi_ss",
    # bias region B_anti: anti-isolated electron and anti-isolated muon
    "emu_ba_os": "dileptons_ba",
    "emu_ba_ss": "dileptons_ba_ss",
}

JETS = {
    "emu": "jets",
    "emu_nodzeta": "jets",
    "emu_nobveto": "jets",
    "emu_ttbar_cr": "jets",
    "emu_diboson_cr": "jets",
    "emu_ss": "jets",
    "emu_dr_os": "jets_dr",
    "emu_dr_ss": "jets_dr",
    "emu_bi_os": "jets_bi",
    "emu_bi_ss": "jets_bi",
    "emu_ba_os": "jets_ba",
    "emu_ba_ss": "jets_ba",
}


def pair_key(category):
    """Name of the dilepton collection of this category."""
    try:
        return PAIRS[category]
    except KeyError:
        raise KeyError(
            f"category '{category}' has no entry in PAIRS (analysis/selections/emu_regions.py)"
        ) from None


def jets_key(category):
    """Name of the jet collection (cleaned against this region's leptons) of this category."""
    try:
        return JETS[category]
    except KeyError:
        raise KeyError(
            f"category '{category}' has no entry in JETS (analysis/selections/emu_regions.py)"
        ) from None
