import subprocess
import argparse
import yaml
import os
import glob
import shutil


def update_config(input_folder, output_folder, merge_mode):
    with open('configs/harvestDatacards.yml', 'r') as file:
        config = yaml.safe_load(file)

    config['input_folder'] = input_folder
    config['output_folder'] = output_folder
    config['channels'] = "tt" # Can be changed to "all" if needed
    config['merge_mode'] = merge_mode # 1 for symmetrised templates, 0 for non-symmetrised templates 

    with open('configs/harvestDatacards.yml', 'w') as file:
        yaml.dump(config, file, sort_keys=False)


def main():
    parser = argparse.ArgumentParser(description="The lazy way to go from the TIDAL output to likelihood scans.")
    parser.add_argument("--input", "-i", required=True, help="Path (absolute) to the input folder containing added_histo.root")
    parser.add_argument("--output", "-o", required=True, help="Path (relative) to the output folder where the results will be stored.")
    parser.add_argument("--alpha", action='store_true', help="Run alpha scan")
    parser.add_argument("--muvsmu", action='store_true', help="Run mu ggH vs mu V scan")
    parser.add_argument("--qqH", action='store_true', help="Run 2 signal category datacard harvesting")
    parser.add_argument("--lumi", default=None, help="Luminosity to use in plot labels")
    parser.add_argument("--channel", required=False, help="Channel to determine label if not in output name")
    parser.add_argument("--symmetrise", action='store_true', help="Produces symmetrised templates")
    args = parser.parse_args()

    if args.symmetrise:
        if not args.qqH:
            subprocess.run([
                "python3", "scripts/convertCards.py",
                "-f", f"{args.input}/added_histo.root"
            ], check=True)
        else:
            subprocess.run([
                "python3", "scripts/convertCards_qqH.py",
                "-f", f"{args.input}/added_histo.root"
            ], check=True)
        update_config(args.input, args.output, merge_mode=1)
    else:
        update_config(args.input, args.output, merge_mode=0)

    # Produce text datacards
    if args.qqH:
        subprocess.run([
            "python3", "scripts/harvestDatacards_qqH.py",
            "-c", "configs/harvestDatacards.yml"
        ], check=True)
    else:
        subprocess.run([
            "python3", "scripts/harvestDatacards.py",
            "-c", "configs/harvestDatacards.yml"
        ], check=True)

    # Creating workspaces
    sub_dir = "cmb" if not args.channel else args.channel
    subprocess.run([
        "combineTool.py", 
        "-m", "125",
        "-M", "T2W",
        "-P", "CombineHarvester.Combine_HtautauCP.CPMixtureDecays:CPMixtureDecays",
        "-i", f"{args.output}/{sub_dir}",
        "-o", "ws.root",
        "--parallel", "8"
    ], check=True)

    if args.alpha:

        # Run maximum likelihood fits
        subprocess.run([
            "combineTool.py",
            "-m", "125",
            "-M", "MultiDimFit",
            "--setParameters", "muV=1,alpha=0,muggH=1,mutautau=1",
            "--setParameterRanges", "alpha=-90,90",
            "--points", "41",
            "--redefineSignalPOIs", "alpha",
            "-d", f"{args.output}/{sub_dir}/ws.root",
            "--algo", "grid",
            "-t", "-1",
            "--there",
            "-n", ".alpha",
            "--alignEdges", "1"
        ], check=True)

        # Make plot of alpha scan
        command = [
            "python3", "scripts/plot1DScan.py",
            "--main", f"{args.output}/{sub_dir}/higgsCombine.alpha.MultiDimFit.mH125.root",
            "--POI", "alpha",
            "--output", f"{args.output}/{sub_dir}/alpha_{sub_dir}",
            "--no-numbers",
            "--no-box",
            "--x-min=-90",
            "--x-max=90",
            "--y-max=20"
        ]

        if args.lumi:
            command.append("--luminosity")
            command.append(f"{args.lumi}")

        if args.channel:
            command.append("--channel")
            command.append(args.channel)

        subprocess.run(command, check=True)
        
    if args.muvsmu:

        out_dir = os.path.abspath(os.path.join(args.output, sub_dir))
        ws_path = os.path.join(out_dir, "ws.root")
        if not os.path.exists(ws_path):
            raise FileNotFoundError(f"Workspace not found at {ws_path}")

        # Run maximum likelihood fits
        subprocess.run([
            "combineTool.py",
            "-m", "125",
            "-M", "MultiDimFit",
            "--setParameters", "muV=1,alpha=0,muggH=1,mutautau=1",
            "--redefineSignalPOIs", "muV,muggH",
            "--setParameterRanges", "muV=-5,4:muggH=-2,4",
            "--points", "2000",
            "-d", ws_path,
            "--algo", "grid",
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
        ], cwd=out_dir, check=True)

        # Wait for condor jobs to finish (this is a simple approach, can be improved by checking job status)
        print("Waiting for condor jobs to finish...")
        log_files = sorted(glob.glob(os.path.join(out_dir, "*.log")))
        for log_file in log_files:
            subprocess.run(["condor_wait", log_file], check=True)

        # Combine the outputs from condor jobs
        point_files = sorted(glob.glob(os.path.join(out_dir, "higgsCombine.muVsmu.POINTS*.MultiDimFit.mH125.root")))
        if point_files:
            subprocess.run([
                "hadd",
                "-v", "1",
                "-f", "higgsCombine.muVsmu.MultiDimFit.mH125",
                *point_files
            ], cwd=out_dir, check=True)
        else:
            print("No POINTS files found for hadd.")

        # Clean up condor output
        condor_dir = os.path.join(out_dir, ".condor")
        condor_log_dir = os.path.join(out_dir, ".condor", "log")
        condor_out_dir = os.path.join(out_dir, ".condor", "out")
        condor_err_dir = os.path.join(out_dir, ".condor", "err")
        os.makedirs(condor_dir, exist_ok=True)
        os.makedirs(condor_log_dir, exist_ok=True)
        os.makedirs(condor_out_dir, exist_ok=True)
        os.makedirs(condor_err_dir, exist_ok=True)

        for path in glob.glob(os.path.join(out_dir, "*.sh")):
            shutil.move(path, os.path.join(condor_dir, os.path.basename(path)))
        for path in glob.glob(os.path.join(out_dir, "*.sub")):
            shutil.move(path, os.path.join(condor_dir, os.path.basename(path)))
        for path in glob.glob(os.path.join(out_dir, "*.log")):
            shutil.move(path, os.path.join(condor_log_dir, os.path.basename(path)))
        for path in glob.glob(os.path.join(out_dir, "*.out")):
            shutil.move(path, os.path.join(condor_out_dir, os.path.basename(path)))
        for path in glob.glob(os.path.join(out_dir, "*.err")):
            shutil.move(path, os.path.join(condor_err_dir, os.path.basename(path)))

        # Clean up individual point outputs
        multidimfit_dir = os.path.join(out_dir, ".MultiDimFit")
        os.makedirs(multidimfit_dir, exist_ok=True)
        for path in glob.glob(os.path.join(out_dir, "higgsCombine.muVsmu.POINTS.*")):
            shutil.move(path, os.path.join(multidimfit_dir, os.path.basename(path)))

        # Make plot of mu ggH vs mu V scan
        subprocess.run([
            "python3", "scripts/plot_2D_scans.py",
            "--file", f"{args.output}/{sub_dir}/higgsCombine.muVsmu.MultiDimFit.mH125",
            "--muvsmu"
        ], check=True)


if __name__ == "__main__":
    main()