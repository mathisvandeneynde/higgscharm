import argparse
import subprocess
from pathlib import Path
from analysis.filesets.utils import check_nano_version


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-y",
        "--year",
        dest="year",
        type=str,
        choices=[
            "2016preVFP",
            "2016postVFP",
            "2017",
            "2018",
            "2022preEE",
            "2022postEE",
            "2023preBPix",
            "2023postBPix",
        ],
    )
    parser.add_argument(
        "--image",
        dest="image",
        type=str,
        default="/cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-dask-almalinux9:2025.10.1-py3.10",        
    )
    parser.add_argument(
        "--samples",
        nargs="*",
        type=str,
        help="(Optional) List of samples to use. If omitted, all available samples will be used",
    )
    parser.add_argument(
        "--nanov",
        dest="nanov",
        type=str,
        choices=["9", "12", "15"],
        default="12",
        help="NanoAOD version",
    )
    args = parser.parse_args()

    try:
        subprocess.run("voms-proxy-info -exists", shell=True, check=True)
    except subprocess.CalledProcessError:
        raise Exception(
            "VOMS proxy expired or non-existing: please run 'voms-proxy-init --voms cms'"
        )

    check_nano_version(args.year, args.nanov)

    sites_file = Path.cwd() / "analysis" / "filesets" / f"{args.year}_sites.yaml"
    if not sites_file.exists():
        cmd = f"python3 analysis/filesets/build_sites.py --year {args.year}"
        subprocess.run(cmd, shell=True)


    samples_str = " ".join(args.samples) if args.samples else ""
    cmd = f"singularity exec -B /afs -B /cvmfs {args.image} python3 analysis/filesets/make_filesets.py --year {args.year} --nanov {args.nanov} --samples {samples_str}"
    subprocess.run(cmd, shell=True)
