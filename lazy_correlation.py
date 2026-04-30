"""
This script assumes that the Harvester step has already been run, and hence
begins from the T2W step.
"""

import subprocess
import argparse
import os
import glob

parser = argparse.ArgumentParser()
parser.add_argument("--datacards", "-d", required=True, help="Path (relative) to the input folder containing the datacards.")
args = parser.parse_args()

subprocess.run([
        "combineTool.py", 
        "-m", "125",
        "-M", "T2W",
        "-P", "CombineHarvester.Combine_HtautauCP.CPMixtureDecays:CPMixtureDecays",
        "-i", args.datacards,
        "-o", "ws.root",
        "--parallel", "8"
    ], check=True)

cwd = os.path.abspath(args.datacards)
ws_path = os.path.join(cwd, "ws.root")
if not os.path.isfile(ws_path):
    raise FileNotFoundError(f"Workspace file not found at {ws_path}. Please check the path and ensure that the T2W step has been completed successfully.")

subprocess.run([
            "combineTool.py",
            "-m", "125",
            "-M", "MultiDimFit",
            "--setParameters", "muV=1,alpha=0,muggH=1,mutautau=1",
            "--redefineSignalPOIs", "muV,muggH",
            "--setParameterRanges", "muV=-5,4:muggH=-2,4",
            "--algo", "none",
            "--saveFitResult",
            "--robustHesse=1",
            "-d", ws_path,
            "-t", "-1",
            "--there",
            "-n", ".muVsmu",
            "--alignEdges", "1",
            "--cminDefaultMinimizerStrategy=0",
            "--cminDefaultMinimizerTolerance=0.1",
            "--cminFallbackAlgo", "Minuit2,Migrad,0:1",
            "--cminFallbackAlgo", "Minuit2,Migrad,0:2",
            "--cminFallbackAlgo", "Minuit2,Migrad,0:4",
            "--cminFallbackAlgo", "Minuit2,Migrad,0:10",
            "--job-mode", "condor",
            "--task-name", "condor-run3-muVsmu",
            "--sub-opts=+MaxRuntime=10799",
            "--split-points", "5"
        ], cwd=cwd, check=True)

# Wait for the condor jobs to finish before proceeding
print("Waiting for condor jobs to finish...")
log_files = sorted(glob.glob(os.path.join(cwd, "*.log")))
for log_file in log_files:
    subprocess.run(["condor_wait", log_file], check=True)

subprocess.run([
    "python3", "correlation_np.py",
    "--root_file", os.path.join(cwd, "multidimfit.muVsmu.root"),
    "--output_dir", cwd
], check=True)
