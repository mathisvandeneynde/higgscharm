import copy
import correctionlib
import numpy as np
import awkward as ak
from pathlib import Path


def get_lepton_os_fr(lepton, year, channel):
    barrel_cut = {"electron": 1.479, "muon": 1.2}
    pt_cut = {"electron": 7.0, "muon": 5.0}
    pdgid_cut = {"electron": 11, "muon": 13}

    json_path = Path.cwd() / "analysis" / "data" / f"{year}_os_fr_weight.json.gz"
    cset = correctionlib.CorrectionSet.from_file(str(json_path))

    lepton_pt = lepton.pt
    lepton_is_barrel = np.abs(lepton.eta) < barrel_cut[channel]

    lepton_in_binning = (
        (lepton_pt >= pt_cut[channel])
        & (lepton_pt <= 80.0)
        & (np.abs(lepton.pdgId) == pdgid_cut[channel])
    )
    lepton_in_pt = ak.fill_none(lepton_pt.mask[lepton_in_binning], 10.0)
    lepton_in_is_barrel = ak.fill_none(lepton_is_barrel.mask[lepton_in_binning], True)

    lepton_fr = cset["os_fake_rate"].evaluate(
        lepton_in_pt, ak.values_astype(lepton_in_is_barrel, float), channel
    )
    return ak.fill_none(
        ak.where(lepton_in_binning, lepton_fr, ak.zeros_like(lepton_fr)), 0
    )


def get_os_fr(lepton, year):
    muon_fr = get_lepton_os_fr(lepton=lepton, year=year, channel="muon")
    electron_fr = get_lepton_os_fr(lepton=lepton, year=year, channel="electron")
    return ak.where(np.abs(lepton.pdgId) == 11, electron_fr, muon_fr)


def add_os_2p2f_weights(events, year, weights_container):
    """Add wweights needed to compute the 2P2F component of the Z+X background estimation"""
    is_2fcr = ak.firsts(events.selected_best_zllcandidate_2fcr.z2.is_2fcr)
    lep1 = ak.firsts(events.selected_best_zllcandidate_2fcr.z2.l1)
    lep2 = ak.firsts(events.selected_best_zllcandidate_2fcr.z2.l2)

    fr_lep1 = get_os_fr(lep1, year)
    fr_lep2 = get_os_fr(lep2, year)

    weight_lep1 = fr_lep1 / (1 - fr_lep1)
    weight_lep2 = fr_lep2 / (1 - fr_lep2)

    weights_container.add(
        name="os_2p2f_w1plusw2",
        weight=ak.where(is_2fcr, weight_lep1 + weight_lep2, 1),
    )
    weights_container.add(
        name="os_2p2f_w1timesw2",
        weight=ak.where(is_2fcr, weight_lep1 * weight_lep2, 1),
    )


def add_os_3p1f_weights(events, year, weights_container):
    """Add wweights needed to compute the 3P1F component of the Z+X background estimation"""
    is_1fcr = ak.firsts(events.selected_best_zllcandidate_1fcr.z2.is_1fcr)
    lep1 = ak.firsts(events.selected_best_zllcandidate_1fcr.z2.l1)
    lep2 = ak.firsts(events.selected_best_zllcandidate_1fcr.z2.l2)

    lep = ak.where(lep1.is_relaxed, lep1, lep2)
    fr = get_os_fr(lep, year)
    weights_container.add(
        name="os_3p1f",
        weight=ak.where(is_1fcr, fr / (1 - fr), 1),
    )