from prettytable import PrettyTable
import subprocess
import argparse
import yaml

BLUE = "\033[94m"
END = "\033[0m"


def get_args():
    parser = argparse.ArgumentParser(description="Run cut optimisation")
    parser.add_argument("--input", "-i", type=str, default="/vols/cms/dmw25/TIDAL/Draw/plots/cuts_for_LateRun3")
    parser.add_argument("--output", "-o", type=str, default="output/datacards/cuts_for_LateRun3/full")
    parser.add_argument("--step", type=str, required=True, help="'fit', 'plot', and then 'summarise'")
    parser.add_argument("--channels", type=str, required=True, help="e.g. et,mt,tt")
    args = parser.parse_args()
    args.channels = args.channels.split(",")
    return args


def to_string(value: float) -> str:
    return str(value).replace(".", "p")


def update_config(input: str, output: str, channels: list[str]) -> None:
    config = {
        "input_folder": input,
        "output_folder": output,
        "channels": channels,
        "merge_mode": 1
    }
    with open("configs/harvestDatacards.yml", "w") as f:
        yaml.dump(config, f, sort_keys=False)


def notify(message: str) -> None:
    print("—" * 80)
    print(f"{BLUE}{message}{END}")
    print("—" * 80)


def submit_fit_job(IPsig, Esplit, output_dir, ch, args):
    input_dir = f"{args.input}/IPsig_{IPsig}_Esplit_{Esplit}/Combined"

    # —————————————————————————————————————————————————————————————————————————
    # Step 1: Symmetrise bins
    # —————————————————————————————————————————————————————————————————————————
    notify(f"Symmetrising bins for IPsig={IPsig}, Esplit={Esplit}...")
    for channel in args.channels:
        symmetrise_command = [
        "python3", "scripts/convertCards.py",
        "-f", f"{input_dir}/added_histo_{channel}.root", 
        ]
        subprocess.run(symmetrise_command, check=True)
    
    # —————————————————————————————————————————————————————————————————————————
    # Step 2: Update config
    # —————————————————————————————————————————————————————————————————————————
    update_config(input_dir, output_dir, args.channels)

    # —————————————————————————————————————————————————————————————————————————
    # Step 3: Run harvestDatacards.py
    # —————————————————————————————————————————————————————————————————————————
    notify("Running harvestDatacards.py...")
    harvest_command = [
        "python3", "scripts/harvestDatacards.py",
        "-c", "configs/harvestDatacards.yml"
    ]
    subprocess.run(harvest_command, check=True)

    # —————————————————————————————————————————————————————————————————————————
    # Step 4: Create workspaces
    # —————————————————————————————————————————————————————————————————————————
    notify("Creating workspaces...")
    t2w_command = [
        "combineTool.py",
        "-m", "125",
        "-M", "T2W",
        "-P", "CombineHarvester.Combine_HtautauCP.CPMixtureDecays:CPMixtureDecays",
        "-i", f"{output_dir}/{ch}",
        "-o", "ws.root",
        "--parallel", "8"
    ]
    subprocess.run(t2w_command, check=True)

    # —————————————————————————————————————————————————————————————————————————
    # Step 5: Submitting fit job to condor
    # —————————————————————————————————————————————————————————————————————————
    notify("Submitting fit job to condor...")
    job_dir = f"{output_dir}/{ch}"
    fit_command = [
        "combineTool.py",
        "-m", "125",
        "-M", "MultiDimFit",
        "--setParameters", "muV=1,alpha=0,muggH=1,mutautau=1",
        "--setParameterRanges", "alpha=-90,90",
        "--points", "51",
        "--redefineSignalPOIs", "alpha",
        "-d", f"{job_dir}/ws.root",
        "--algo", "grid",
        "-t", "-1",
        "--there",
        "-n", ".alpha",
        "--alignEdges", "1",
        "--job-mode", "condor",
        "--task-name", f"IPsig_{IPsig}_Esplit_{Esplit}",
        "--sub-opts", "+MaxRuntime=3600",
    ]
    subprocess.run(fit_command, check=True, cwd=job_dir)


def plot_alpha_scan(output_dir, ch):
    notify("Making plot of alpha scan...")
    plot_command = [
        "python3", "scripts/plot1DScan.py",
        f"--main={output_dir}/{ch}/higgsCombine.alpha.MultiDimFit.mH125.root",
        "--POI=alpha",
        f"--output={output_dir}/alpha_{ch}",
        "--no-numbers",
        "--no-box",
        "--x-min=-90",
        "--x-max=90",
        "--y-max=30"
    ]
    plot_output = subprocess.run(plot_command, capture_output=True, text=True)
    with open(f"{output_dir}/RESULT.txt", "w") as f:
        f.write(plot_output.stdout)
    print(f"Output written to {output_dir}/RESULT.txt")


def extract_sensitivity(result_file):
    notify(f"Extracting sensitivity from {result_file}...")
    with open(result_file, "r") as f:
        lines = f.readlines()

    sig_lines = []
    for line in lines:
        if "max sigma" in line:
            cp_odd_excl = float(line.split()[-1])
        if "valid_lo" in line:
            sig_lines.append(line)
    one_sig_line = sig_lines[0]
    one_sig_lo = -1 * float(one_sig_line.split()[1][:-1])
    one_sig_hi = float(one_sig_line.split()[3][:-1])
    alpha_err = (one_sig_hi + one_sig_lo) / 2
    
    return alpha_err, cp_odd_excl
    

def main(args):
    IPsig_values = [1.25, 1.35, 1.45, 1.55, 1.65]
    Esplit_values = [0.1, 0.125, 0.15, 0.175, 0.2]

    if args.step == "summarise":
        table = PrettyTable()
        table.field_names = ["cut_IPsig", "cut_Esplit", "alpha_err", "cp_odd_excl"]

    for IPsig_ in IPsig_values:
        for Esplit_ in Esplit_values:
            IPsig, Esplit = to_string(IPsig_), to_string(Esplit_)
            output_dir = f"{args.output}/IPsig_{IPsig}_Esplit_{Esplit}"
            ch = "cmb" if len(args.channels) > 1 else args.channels[0]

            if args.step == "fit":
                submit_fit_job(IPsig, Esplit, output_dir, ch, args)

            elif args.step == "plot":
                plot_alpha_scan(output_dir, ch)

            elif args.step == "summarise":
                result_file = f"{output_dir}/RESULT.txt"
                alpha_err, cp_odd_excl = extract_sensitivity(result_file)
                table.add_row([
                    f"{IPsig_:.2f}", f"{Esplit_:.2f}", f"{alpha_err:.4f}", f"{cp_odd_excl:.4f}"
                ])

    if args.step == "summarise":
        print(table)


if __name__ == "__main__":
    args = get_args()
    main(args)

