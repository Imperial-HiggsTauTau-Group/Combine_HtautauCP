import ROOT

# -----------------------------
# user inputs
# -----------------------------
toy=1
workspace_file = "outputs/Feb02_Unblinding/cmb/ws.root"
workspace_name = "w"              # often "w"
toy_file = "outputs/Feb02_Unblinding/cmb/higgsCombine.500Toys.GenerateOnly.mH125.123456.root"
toy_dataset_name = f"toys/toy_{toy}"   # directory + dataset name inside file
output_file = workspace_file.replace(".root", f"_toy{toy}_data.root")  # e.g. ws_with_toy.root

# -----------------------------
# open workspace file
# -----------------------------
f_ws = ROOT.TFile.Open(workspace_file, "READ")
w = f_ws.Get(workspace_name)

if not w:
    raise RuntimeError(f"Workspace '{workspace_name}' not found in {workspace_file}")

print("Loaded workspace:", w.GetName())

d_old = w.data("data_obs")
if not d_old:
    raise RuntimeError("data_obs not found")


print("Old name:", d_old.GetName())

d_old.SetName("data_obs_orig")

print("Renamed old dataset to:", d_old.GetName())


# -----------------------------
# open toy file and get dataset
# -----------------------------
f_toy = ROOT.TFile.Open(toy_file, "READ")
toy_data = f_toy.Get(toy_dataset_name)


if not toy_data:
    raise RuntimeError(f"Toy dataset '{toy_dataset_name}' not found in {toy_file}")

print("Loaded toy dataset:", toy_data.GetName())

# -----------------------------
# clone toy dataset and rename it to data_obs
# -----------------------------
toy_clone = toy_data.Clone("data_obs")

# -----------------------------
# import into workspace
# overwrite existing data_obs if present
# -----------------------------
getattr(w, "import")(toy_clone, ROOT.RooFit.Rename("data_obs"), ROOT.RooFit.RecycleConflictNodes())

print("Imported toy dataset as data_obs into workspace.")

# -----------------------------
# write new workspace file
# -----------------------------
f_out = ROOT.TFile.Open(output_file, "RECREATE")
w.Write()
f_out.Close()

print("Wrote modified workspace to:", output_file)

# cleanup
f_ws.Close()
f_toy.Close()
