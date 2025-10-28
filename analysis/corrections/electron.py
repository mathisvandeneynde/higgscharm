import correctionlib
import numpy as np
import awkward as ak
from typing import Type
from coffea.analysis_tools import Weights
from analysis.selections.trigger import trigger_match_mask
from analysis.selections.event_selections import get_trigger_mask
from analysis.corrections.met import update_met
from analysis.corrections.utils import (
    get_pog_json,
    get_egamma_json,
    unflat_sf,
    get_electron_hlt_json,
)


class ElectronWeights:
    """
    Electron ID, Reco and HLT weights class

    Parameters:
    -----------
        events:
            events collection
        weights:
            Weights container
        year:
            Year of the dataset {2016preVFP, 2016postVFP, 2017, 2018, 2022preEE, 2022postEE, 2023preBPix, 2023postBPix}
        variation:
            syst variation
        id_wp:
            ID working point {wp80iso, wp90iso}
        nano_version:
            NanoAOD version. '9' for Run2 datasets, and '12' for Run3 datasets

    more info: https://twiki.cern.ch/twiki/bin/view/CMS/EgammSFandSSRun3#Scale_factors_and_correction_AN1
    """

    def __init__(
        self,
        events: ak.Array,
        weights: Type[Weights],
        year: str,
        nano_version: str,
        variation: str,
    ) -> None:
        self.events = events
        self.electrons = events.selected_electrons
        self.weights = weights
        self.year = year
        self.variation = variation
        self.nano_version = nano_version

        self.flat_electrons = ak.flatten(self.electrons)
        self.electrons_counts = ak.num(self.electrons)

        # set id working points map
        self.id_map = {
            "wp80iso": "wp80iso",
            "wp90iso": "wp90iso",
        }
        self.year_map = {
            "2016preVFP": "2016preVFP",
            "2016postVFP": "2016postVFP",
            "2017": "2017",
            "2018": "2018",
            "2022postEE": "2022Re-recoE+PromptFG",
            "2022preEE": "2022Re-recoBCD",
            "2023preBPix": "2023PromptC",
            "2023postBPix": "2023PromptD",
        }
        self.year_key = year[:4]

    def add_id_weights(self, id_wp):
        """
        add electron ID weights to weights container
        """
        if self.nano_version == "12":
            nominal_weights = self.get_id_weights_run3(variation="sf", id_wp=id_wp)
            up_weights = self.get_id_weights_run3(variation="sfup", id_wp=id_wp)
            down_weights = self.get_id_weights_run3(variation="sfdown", id_wp=id_wp)
        else:
            nominal_weights = self.get_id_weights_run2(variation="sf", id_wp=id_wp)
            up_weights = self.get_id_weights_run2(variation="sfup", id_wp=id_wp)
            down_weights = self.get_id_weights_run2(variation="sfdown", id_wp=id_wp)
        if self.variation == "nominal":
            # add scale factors to weights container
            self.weights.add(
                name=f"CMS_eff_e_id_{self.year_key}",
                weight=nominal_weights,
                weightUp=up_weights,
                weightDown=down_weights,
            )
        else:
            self.weights.add(
                name=f"CMS_eff_e_id_{self.year_key}",
                weight=nominal_weights,
            )

    def add_reco_weights(self, reco_range):
        """
        add electron Reco weights to weights container
        """
        var_naming_map = {
            "RecoBelow20": f"CMS_eff_e_reco_below20_{self.year_key}",
            "RecoAbove20": f"CMS_eff_e_reco_above20_{self.year_key}",
            "Reco20to75": f"CMS_eff_e_reco_20to75_{self.year_key}",
            "RecoAbove75": f"CMS_eff_e_reco_above75_{self.year_key}",
        }
        if self.nano_version == "12":
            nominal_weights = self.get_reco_weights_run3(
                variation="sf", reco_range=reco_range
            )
            up_weights = self.get_reco_weights_run3(
                variation="sfup", reco_range=reco_range
            )
            down_weights = self.get_reco_weights_run3(
                variation="sfdown", reco_range=reco_range
            )
        else:
            nominal_weights = self.get_reco_weights_run2(
                variation="sf", reco_range=reco_range
            )
            up_weights = self.get_reco_weights_run2(
                variation="sfup", reco_range=reco_range
            )
            down_weights = self.get_reco_weights_run2(
                variation="sfdown", reco_range=reco_range
            )
        if self.variation == "nominal":
            # add scale factors to weights container
            self.weights.add(
                name=var_naming_map[reco_range],
                weight=nominal_weights,
                weightUp=up_weights,
                weightDown=down_weights,
            )
        else:
            self.weights.add(
                name=var_naming_map[reco_range],
                weight=nominal_weights,
            )

    def add_hlt_weights(self, id_wp):
        """
        add electron HLT weights to weights container
        """
        if self.nano_version == "12":
            nominal_weights = self.get_hlt_weights_run3(variation="nom", id_wp=id_wp)
        else:
            nominal_weights = self.get_hlt_weights_run2(variation="nom", id_wp=id_wp)
        if self.variation == "nominal":
            # add scale factors to weights container
            self.weights.add(
                name=f"CMS_eff_e_trigger_{self.year_key}",
                weight=nominal_weights,
                # weightUp=up_weights,
                # weightDown=down_weights,
            )
        else:
            self.weights.add(
                name=f"CMS_eff_e_trigger_{self.year_key}",
                weight=nominal_weights,
            )

    def get_id_weights_run3(self, variation, id_wp):
        """Compute electron ID weights for Run3 datasets"""
        # get electron correction set
        cset = correctionlib.CorrectionSet.from_file(
            get_pog_json(json_name="electron_id", year=self.year)
        )
        # get electrons that pass the id wp, and within SF binning
        electron_pt_mask = self.flat_electrons.pt > 10.0
        in_electron_mask = electron_pt_mask
        in_electrons = self.flat_electrons.mask[in_electron_mask]

        # get electrons pT and abseta (replace None values with some 'in-limit' value)
        electron_pt = ak.fill_none(in_electrons.pt, 15.0)
        electron_eta = ak.fill_none(in_electrons.eta, 0)
        electron_phi = ak.fill_none(in_electrons.phi, 0)

        cset_args = [
            self.year_map[self.year],
            variation,
            self.id_map[id_wp],
            electron_eta,
            electron_pt,
        ]
        if self.year.startswith("2023"):
            cset_args += [electron_phi]
        weights = unflat_sf(
            cset["Electron-ID-SF"].evaluate(*cset_args),
            in_electron_mask,
            self.electrons_counts,
        )
        return weights

    def get_id_weights_run2(self, variation, id_wp):
        """Compute electron ID weights for Run2 datasets"""
        # get electron correction set
        cset = correctionlib.CorrectionSet.from_file(
            get_pog_json(json_name="electron_id", year=self.year)
        )
        # get electrons that pass the id wp, and within SF binning
        electron_pt_mask = (self.flat_electrons.pt > 10.0) & (
            self.flat_electrons.pt < 499.999
        )  # potential problems with pt > 500 GeV
        in_electron_mask = electron_pt_mask
        in_electrons = self.flat_electrons.mask[in_electron_mask]

        # get electrons pT and abseta (replace None values with some 'in-limit' value)
        electron_pt = ak.fill_none(in_electrons.pt, 15.0)
        electron_eta = ak.fill_none(in_electrons.eta + in_electrons.deltaEtaSC, 0)
        electron_phi = ak.fill_none(in_electrons.phi, 0)

        cset_args = [
            self.year_map[self.year],
            variation,
            self.id_map[id_wp],
            electron_eta,
            electron_pt,
        ]
        if self.year.startswith("2023"):
            cset_args += [electron_phi]
        weights = unflat_sf(
            cset["UL-Electron-ID-SF"].evaluate(*cset_args),
            in_electron_mask,
            self.electrons_counts,
        )
        return weights

    def get_reco_weights_run3(self, variation, reco_range):
        """Compute electron Reco weights for Run3 datasets"""
        # get electron correction set
        cset = correctionlib.CorrectionSet.from_file(
            get_pog_json(json_name="electron_id", year=self.year)
        )
        # get electrons that pass the id wp, and within SF binning
        electron_pt_mask = {
            "RecoBelow20": (self.flat_electrons.pt > 10.0)
            & (self.flat_electrons.pt < 20.0),
            "Reco20to75": (self.flat_electrons.pt > 20.0)
            & (self.flat_electrons.pt < 75.0),
            "RecoAbove75": self.flat_electrons.pt > 75,
        }
        in_electrons_mask = electron_pt_mask[reco_range]
        in_electrons = self.flat_electrons.mask[in_electrons_mask]

        # get electrons pT and abseta (replace None values with some 'in-limit' value)
        electron_pt_limits = {
            "RecoBelow20": 15,
            "Reco20to75": 30,
            "RecoAbove75": 80,
        }
        electron_pt = ak.fill_none(in_electrons.pt, electron_pt_limits[reco_range])
        electron_eta = ak.fill_none(in_electrons.eta, 0)
        electron_phi = ak.fill_none(in_electrons.phi, 0)

        cset_args = [
            self.year_map[self.year],
            variation,
            reco_range,
            electron_eta,
            electron_pt,
        ]
        if self.year.startswith("2023"):
            cset_args += [electron_phi]
        weights = unflat_sf(
            cset["Electron-ID-SF"].evaluate(*cset_args),
            in_electrons_mask,
            self.electrons_counts,
        )
        return weights

    def get_reco_weights_run2(self, variation, reco_range):
        """Compute electron Reco weights for Run2 datasets"""
        # get electron correction set
        cset = correctionlib.CorrectionSet.from_file(
            get_pog_json(json_name="electron_id", year=self.year)
        )
        # get electrons that pass the id wp, and within SF binning
        electron_pt_mask = {
            "RecoAbove20": (self.flat_electrons.pt > 20)
            & (self.flat_electrons.pt < 499.999),
            "RecoBelow20": (self.flat_electrons.pt > 10)
            & (self.flat_electrons.pt < 20),
        }
        in_electrons_mask = electron_pt_mask[reco_range]
        in_electrons = self.flat_electrons.mask[in_electrons_mask]

        # get electrons pT and abseta (replace None values with some 'in-limit' value)
        electrons_pt_limits = {
            "RecoAbove20": 21,
            "RecoBelow20": 15,
            "Reco20to75": 30,
            "RecoAbove75": 80,
        }
        electron_pt = ak.fill_none(in_electrons.pt, electrons_pt_limits[reco_range])
        electron_eta = ak.fill_none(in_electrons.eta + in_electrons.deltaEtaSC, 0)
        electron_phi = ak.fill_none(in_electrons.phi, 0)

        cset_args = [
            self.year_map[self.year],
            variation,
            reco_range,
            electron_eta,
            electron_pt,
        ]
        weights = unflat_sf(
            cset["UL-Electron-ID-SF"].evaluate(*cset_args),
            in_electrons_mask,
            self.electrons_counts,
        )
        return weights

    def get_hlt_weights_run3(self, variation, id_wp):
        """Compute electron HLT weights for Run3 datasets"""
        # get electron correction set
        cset = correctionlib.CorrectionSet.from_file(
            get_pog_json(json_name="electron_hlt", year=self.year)
        )
        # get electrons that pass the id wp, and within SF binning
        electron_pt_mask = self.flat_electrons.pt > 25.0

        in_electrons_mask = electron_pt_mask
        in_electrons = self.flat_electrons.mask[in_electrons_mask]

        # get electrons pT and abseta (replace None values with some 'in-limit' value)
        electron_pt = ak.fill_none(in_electrons.pt, 25)
        electron_eta = ak.fill_none(in_electrons.eta, 0)

        hlt_path_id_map = {
            "wp80iso": "HLT_SF_Ele30_MVAiso80ID",
            "wp90iso": "HLT_SF_Ele30_MVAiso90ID",
        }
        sf_variations_map = {"nom": "sf", "up": "sfup", "down": "sfdown"}

        kind = "single" if ak.all(ak.num(self.electrons) == 1) else "double"
        if kind == "single":
            # for single electorn events, compute SF from POG SF
            sf = cset["Electron-HLT-SF"].evaluate(
                self.year_map[self.year],
                sf_variations_map[variation],
                hlt_path_id_map[id_wp],
                electron_eta,
                electron_pt,
            )
            sf = ak.where(in_electrons_mask, sf, ak.ones_like(sf))
            sf = ak.fill_none(ak.unflatten(sf, self.electrons_counts), value=1)
            nominal_sf = ak.firsts(nominal_sf)

        elif kind == "double":
            # for double electron events, compute SF from electrons' efficiencies
            data_eff = cset["Electron-HLT-DataEff"].evaluate(
                self.year_map[self.year],
                variation,
                hlt_path_id_map[id_wp],
                electron_eta,
                electron_pt,
            )
            data_eff = ak.where(in_electrons_mask, data_eff, ak.ones_like(data_eff))
            data_eff = ak.unflatten(data_eff, self.electrons_counts)
            data_eff_leading = ak.firsts(data_eff)
            data_eff_subleading = ak.pad_none(data_eff, target=2)[:, 1]
            full_data_eff = (
                data_eff_leading
                + data_eff_subleading
                - data_eff_leading * data_eff_subleading
            )
            full_data_eff = ak.fill_none(full_data_eff, 1)

            mc_eff = cset["Electron-HLT-McEff"].evaluate(
                self.year_map[self.year],
                variation,
                hlt_path_id_map[id_wp],
                electron_eta,
                electron_pt,
            )
            mc_eff = ak.where(in_electrons_mask, mc_eff, ak.ones_like(mc_eff))
            mc_eff = ak.unflatten(mc_eff, self.electrons_counts)
            mc_eff_leading = ak.firsts(mc_eff)
            mc_eff_subleading = ak.pad_none(mc_eff, target=2)[:, 1]
            full_mc_eff = (
                mc_eff_leading + mc_eff_subleading - mc_eff_leading * mc_eff_subleading
            )
            full_mc_eff = ak.fill_none(full_mc_eff, 1)

            # compute SF from efficiencies
            nominal_sf = full_data_eff / full_mc_eff

        return nominal_sf

    def get_hlt_weights_run2(self, variation, id_wp):
        """Compute electron HLT weights for Run2 datasets"""
        hlt_paths = {
            "2016preVFP": "Ele27_WPTight_Gsf_OR_Photon175",
            "2016postVFP": "Ele27_WPTight_Gsf_OR_Photon175",
            "2017": "Ele35_WPTight_Gsf_OR_Photon200",
            "2018": "Ele32_WPTight_Gsf",
        }
        hlt_path_id_map = {
            "wp80iso": "HLT_SF_Ele30_MVAiso80ID",
            "wp90iso": "HLT_SF_Ele30_MVAiso90ID",
        }
        pt_limits = {
            "2016preVFP": 27.0,
            "2016postVFP": 27.0,
            "2017": 35.0,
            "2018": 32.0,
        }
        electron_pt_mask = (self.flat_electrons.pt > pt_limits[self.year]) & (
            self.flat_electrons.pt < 499.9
        )
        eta_mask = (
            np.abs(self.flat_electrons.eta + self.flat_electrons.deltaEtaSC) < 2.5
        )
        in_electrons_mask = electron_pt_mask & eta_mask
        in_electrons = self.flat_electrons.mask[in_electrons_mask]

        # get electrons pT and abseta (replace None values with some 'in-limit' value)
        electron_pt = ak.fill_none(in_electrons.pt, 40)
        electron_eta = ak.fill_none(
            self.flat_electrons.eta + self.flat_electrons.deltaEtaSC, 0
        )
        # check whether there are single or double electron events
        kind = "single" if ak.all(ak.num(self.electrons) == 1) else "double"

        if kind == "single":
            cset = correctionlib.CorrectionSet.from_file(
                get_electron_hlt_json("SF", self.year)
            )
            sf = cset[f"HLT_SF_{hlt_paths[self.year]}_MVAiso80ID"].evaluate(
                electron_eta, electron_pt
            )
            sf = ak.where(in_electrons_mask, sf, ak.ones_like(sf))
            sf = ak.fill_none(ak.unflatten(sf, self.electrons_counts), value=1)
            nominal_sf = ak.firsts(sf)

        elif kind == "double":
            cset = correctionlib.CorrectionSet.from_file(
                get_electron_hlt_json("DataEff", self.year)
            )
            data_eff = cset[f"HLT_DataEff_{hlt_paths[self.year]}_MVAiso80ID"].evaluate(
                electron_eta, electron_pt
            )
            data_eff = ak.where(in_electrons_mask, data_eff, ak.ones_like(data_eff))
            data_eff = ak.unflatten(data_eff, self.electrons_counts)
            data_eff_leading = ak.firsts(data_eff)
            data_eff_subleading = ak.pad_none(data_eff, target=2)[:, 1]
            full_data_eff = (
                data_eff_leading
                + data_eff_subleading
                - data_eff_leading * data_eff_subleading
            )
            full_data_eff = ak.fill_none(full_data_eff, 1)

            cset = correctionlib.CorrectionSet.from_file(
                get_electron_hlt_json("MCEff", self.year)
            )
            mc_eff = cset[f"HLT_MCEff_{hlt_paths[self.year]}_MVAiso80ID"].evaluate(
                electron_eta, electron_pt
            )
            mc_eff = ak.where(in_electrons_mask, mc_eff, ak.ones_like(mc_eff))
            mc_eff = ak.unflatten(mc_eff, self.electrons_counts)
            mc_eff_leading = ak.firsts(mc_eff)
            mc_eff_subleading = ak.pad_none(mc_eff, target=2)[:, 1]
            full_mc_eff = (
                mc_eff_leading + mc_eff_subleading - mc_eff_leading * mc_eff_subleading
            )
            full_mc_eff = ak.fill_none(full_mc_eff, 1)

            nominal_sf = full_data_eff / full_mc_eff

        return nominal_sf
