import argparse
import os
import subprocess
import sys

# gets changed automatically by this script
HARVEST_CFG = "configs/harvestDatacards.yml"
# path to the base directory where the datacards are located (TIDAL to avoid copying stuff)
DATACARD_BASE = "/vols/cms/lcr119/offline/HiggsCP/TIDAL/Draw/Plots/TauIDStudyDatacards"


def run(cmd):
    print(f"\n>>> {cmd}")
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        print(f"\nERROR: command failed with exit code {result.returncode}")
        sys.exit(result.returncode)


def update_harvest_config(output, channel="tt"):
    cfg = os.path.join(os.path.dirname(__file__), "..", HARVEST_CFG)
    cfg = os.path.normpath(cfg)
    content = (
        f'input_folder: "{DATACARD_BASE}/{output}"\n'
        f'output_folder: "output_26June_OptTauID_{channel}/{output}"\n'
        f'channels: "{channel}"\n'
        f'merge_mode: 1\n'
    )
    with open(cfg, "w") as f:
        f.write(content)
    print(f"Updated {HARVEST_CFG} for output: {output}")


def main():
    parser = argparse.ArgumentParser(description="Run all HiggsCP analysis steps")
    parser.add_argument("output", help="Output folder name, e.g. PNet_7")
    parser.add_argument("--channel", default="tt", help="Channel to analyze (default: tt)")
    args = parser.parse_args()
    output = args.output
    channel = args.channel

    update_harvest_config(output, channel)

    run(f"python3 scripts/convertCards.py -f {DATACARD_BASE}/{output}/added_histo_{channel}.root")

    run("python3 scripts/harvestDatacards.py -c configs/harvestDatacards.yml")

    run(
        f"combineTool.py -m 125 -M T2W"
        f" -P CombineHarvester.Combine_HtautauCP.CPMixtureDecays:CPMixtureDecays"
        f" -i output_26June_OptTauID_{channel}/{output}/{channel}"
        f" -o ws.root --parallel 8"
    )

    run(
        f"combineTool.py -m 125 -M MultiDimFit"
        f" --setParameters muV=1,alpha=0,muggH=1,mutautau=1"
        f" --setParameterRanges alpha=-90,90"
        f" --points 100 --redefineSignalPOIs alpha"
        f" -d output_26June_OptTauID_{channel}/{output}/{channel}/ws.root"
        f" --algo grid -t -1 --there -n .alphaEXPECTED --alignEdges 1"
        f" --cminDefaultMinimizerStrategy=0 --cminDefaultMinimizerTolerance=0.1"
        f" --cminFallbackAlgo Minuit2,Migrad,0:1"
        f" --cminFallbackAlgo Minuit2,Migrad,0:2"
        f" --cminFallbackAlgo Minuit2,Migrad,0:4"
        f" --cminFallbackAlgo Minuit2,Migrad,0:10"
    )

    run(
        f"python3 scripts/plot1DScan.py"
        f" --main=output_26June_OptTauID_{channel}/{output}/{channel}/higgsCombine.alphaEXPECTED.MultiDimFit.mH125.root"
        f" --POI=alpha"
        f" --output=output_26June_OptTauID_{channel}/{output}/{channel}/alpha_{channel}_EXPECTED"
        f" --no-numbers --no-box --x-min=-90 --x-max=90 --y-max=9"
    )

    print("\nAll steps completed successfully.")


if __name__ == "__main__":
    main()
