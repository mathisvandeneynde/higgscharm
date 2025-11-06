import correctionlib
import awkward as ak
from analysis.corrections.utils import correction_files


def get_ctag_mask(jets: ak.Array, year: str, wp: str, tagger: str):
    """

    Parameters:
    -----------
      jets: Jet collection
      year: {2016preVFP, 2016postVFP, 2017, 2018, 2022preEE, 2022postEE, 2023preBPix, 2023postBPix, 2024}
      wp: {loose, medium, tight}
      tagger: {deepjet, pnet, upart}
    """
    tagger_map = {
        "deepjet": "deepJet_wp_values",
        "pnet": "particleNet_wp_values",
        "upart": "UParTAK4_wp_values"
    }
    wp_map = {"loose": "L", "medium": "M", "tight": "T"}
    
    cset = correctionlib.CorrectionSet.from_file(correction_files["ctagging"][year])
    ctag_wps_evaluator = cset[tagger_map[tagger]]
    cvsb_wp = ctag_wps_evaluator.evaluate(wp_map[wp], "CvB")
    cvsl_wp = ctag_wps_evaluator.evaluate(wp_map[wp], "CvL")
    
    if tagger == "deepjet": 
        pass_ctag_wp = (jets.btagDeepFlavCvB > cvsb_wp) & (jets.btagDeepFlavCvL > cvsl_wp)
    elif tagger == "pnet": 
        pass_ctag_wp = (jets.btagPNetCvB > cvsb_wp) & (jets.btagPNetCvL > cvsl_wp)
    elif tagger == "upart": 
        pass_ctag_wp = (jets.btagUParTAK4CvB > cvsb_wp) & (jets.btagUParTAK4CvL > cvsl_wp)
    
    return pass_ctag_wp