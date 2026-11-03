import yaml
import json
import numpy as np
import awkward as ak
import importlib.resources
from coffea.lumi_tools import LumiMask
from analysis.selections.trigger import trigger_mask, trigger_match_mask, zzto4l_trigger


def get_lumi_mask(events, year):
    year_map = {
        "2016": "analysis/data/lumimask/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt",
        "2017": "analysis/data/lumimask/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt",
        "2018": "analysis/data/lumimask/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt",
        "2022": "analysis/data/lumimask/Cert_Collisions2022_355100_362760_Golden.txt",
        "2023": "analysis/data/lumimask/Cert_Collisions2023_366442_370790_Golden.txt",
        "2024": "analysis/data/lumimask/Cert_Collisions2024_378981_386951_Golden.txt",
    }
    for key in year_map:
        if year.startswith(key):
            goldenjson = year_map[key]
            break
    else:
        raise ValueError(f"Unrecognized year format: '{year}'")

    if hasattr(events, "genWeight"):
        lumi_mask = np.ones(len(events), dtype="bool")
    else:
        lumi_info = LumiMask(goldenjson)
        lumi_mask = lumi_info(events.run, events.luminosityBlock)
    return lumi_mask == 1


def get_zzto4l_trigger_mask(events, hlt_paths, dataset_key, year):
    return zzto4l_trigger(events, hlt_paths, dataset_key, year)


def get_trigger_mask(events, hlt_paths, dataset_key, year):
    return trigger_mask(events, hlt_paths, dataset_key, year)


def get_trigger_match_mask(events, hlt_paths, year, leptons):
    mask = trigger_match_mask(events, hlt_paths, year, leptons)
    return ak.sum(mask, axis=-1) > 0


def get_metfilters_mask(events, year):
    with importlib.resources.path("analysis.data.met", "metfilters.json") as path:
        with open(path, "r") as handle:
            metfilters = json.load(handle)[year]
    metfilters_mask = np.ones(len(events), dtype="bool")
    metfilterkey = "mc" if hasattr(events, "genWeight") else "data"
    for mf in metfilters[metfilterkey]:
        if mf in events.Flag.fields:
            metfilters_mask = metfilters_mask & events.Flag[mf]
    return metfilters_mask


def get_stitching_mask(events, dataset, dataset_key, ht_value):
    stitching_mask = np.ones(len(events), dtype="bool")
    if dataset.startswith(dataset_key):
        stitching_mask = events.LHE.HT < ht_value
    return stitching_mask

def gen_split_mask(events, objects, dataset):
    """Split DY into ZTT and ZLL using the generator origin of the selected e and mu.

    ZTT: both the selected electron and muon come from a tau decay (genPartFlav == 15).
    ZLL: the complement, so ZTT + ZLL is exactly the full DY sample.
    Returns True everywhere for data and for all samples without _ZTT/_ZLL in their name.
    """
    is_ztt_key = "_ZTT" in dataset
    is_zll_key = "_ZLL" in dataset
    if not (is_ztt_key or is_zll_key):
        return np.ones(len(events), dtype=bool)
    e_flav = ak.firsts(objects["electrons"].genPartFlav)
    mu_flav = ak.firsts(objects["muons"].genPartFlav)
    ztt = ak.to_numpy(ak.fill_none((e_flav == 15) & (mu_flav == 15), False))
    return ztt if is_ztt_key else ~ztt

def no_extra_loose_lepton(objects, pair_key, dr_match=0.01):
    """Veto on extra loose leptons, relative to the selected pair.

    Counts the loose leptons (electrons_loose, muons_loose) that are NOT one of the two
    leptons of the selected pair in `pair_key`, and returns True if there are none.
    A loose lepton belongs to the pair if it has the same flavour and lies within
    dr_match of the pair lepton. Because both come from the same corrected collection,
    a selected lepton lies at dR = 0 from itself.

    The criterion is the same in every region, independent of the isolation of the
    selected leptons (SR/AR: isolated, DR and bias regions: non-isolated).
    Events without a pair return True; the pair requirement removes them anyway.
    """
    pair = ak.firsts(objects[pair_key])
    n_extra = 0
    for key, pdg in (("electrons_loose", 11), ("muons_loose", 13)):
        loose = objects[key]
        is_l1 = (loose.delta_r(pair.l1) < dr_match) & (np.abs(pair.l1.pdgId) == pdg)
        is_l2 = (loose.delta_r(pair.l2) < dr_match) & (np.abs(pair.l2.pdgId) == pdg)
        n_extra = n_extra + ak.sum(~(is_l1 | is_l2), axis=1)
    return ak.to_numpy(ak.fill_none(n_extra == 0, True))