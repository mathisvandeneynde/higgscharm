"""Build the 'QCD' (fake-lepton) process from the fake categories.

In a fake category (e.g. emu_qcd) data have weight +w and MC events -w (times their usual
weight). The fake estimate in the target region is therefore Data + sum(MC) in that category:

    QCD[target, var] = Data[fake_cat, var_data] + sum_MC MC[fake_cat, var]

Data have no systematic shifts of their own: for every variation that is not a fake-factor
variation (CMS_fake_*), the nominal data are used (the data entries under object shifts come
from uncorrected objects and must not be used). For an MC process that lacks a variation,
its nominal is used. The same holds for a variation that exists on the axis but was not
filled in this fake category (e.g. CMS_fake_bias in a category with a fixed bias).
"""
import numpy as np
from analysis.selections.emu_regions import FAKE_TARGETS

FAKE_PREFIX = "CMS_fake_"


def _slicer(hist, category, variation):
    names = list(hist.axes.name)
    sl = [slice(None)] * len(names)
    sl[names.index("category")] = hist.axes["category"].index(category)
    sl[names.index("variation")] = hist.axes["variation"].index(variation)
    return tuple(sl)


def _filled(hist, category, variation):
    """True if this (category, variation) slice of the histogram has any entries."""
    if variation not in set(hist.axes["variation"]):
        return False
    view = hist.view(flow=True)[_slicer(hist, category, variation)]
    return bool(np.any(view.value != 0) or np.any(view.variance != 0))


def add_fake_process(processed_histograms, name="QCD", targets=None):
    """Add processed_histograms[name] built from the fake categories; returns the dict."""
    targets = FAKE_TARGETS if targets is None else targets
    data = processed_histograms["Data"]
    mc_processes = [p for p in processed_histograms if p not in ("Data", name)]
    out = {}
    for hname, hdata in data.items():
        if "category" not in hdata.axes.name or "variation" not in hdata.axes.name:
            continue
        cats = list(hdata.axes["category"])
        if not any(src in cats and tgt in cats for src, tgt in targets.items()):
            continue
        mc_hists = [processed_histograms[p][hname] for p in mc_processes if hname in processed_histograms[p]]
        # template with the largest set of variations
        base = max(mc_hists, key=lambda h: len(h.axes["variation"])) if mc_hists else hdata
        qcd = base.copy()
        qcd.reset()
        qview = qcd.view(flow=True)
        for src, tgt in targets.items():
            if src not in cats or tgt not in cats:
                continue
            for var in qcd.axes["variation"]:
                qsl = _slicer(qcd, tgt, var)
                use = var.startswith(FAKE_PREFIX) and _filled(hdata, src, var)
                dvar = var if use else "nominal"
                dview = hdata.view(flow=True)[_slicer(hdata, src, dvar)]
                qview.value[qsl] += dview.value
                qview.variance[qsl] += dview.variance
                for h in mc_hists:
                    mvar = var if _filled(h, src, var) else "nominal"
                    mview = h.view(flow=True)[_slicer(h, src, mvar)]
                    qview.value[qsl] += mview.value
                    qview.variance[qsl] += mview.variance
        out[hname] = qcd
    if out:
        processed_histograms[name] = out
    return processed_histograms
