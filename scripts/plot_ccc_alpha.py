import argparse
import os
import json

import ROOT
ROOT.gROOT.SetBatch()

import CombineHarvester.CombineTools.plotting as plot


label_dict = {
    "htt_et": "#tau_{e}#tau_{h}",
    "htt_mt": "#tau_{#mu}#tau_{h}",
    "htt_tt": "#tau_{h}#tau_{h}",
}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", type=str,
                        help="Input root file.")
    parser.add_argument("-p", "--poi",
                        type=str,
                        default="alpha",
                        help="POI that is fit.")
    parser.add_argument("-o", "--output",
                        default="ChannelCompatibilityCheck",
                        help="Name of the plot file")
    parser.add_argument("-r", "--frame-range",
                        type=lambda ranges: [float(edge.replace("m", "-")) for edge in ranges.split(",")],
                        default=[-90, 120],
                        help="Range of the fit values")
    parser.add_argument("--separate-nominal-fit",
                        type=str,
                        default=None,
                        help="Fallback solution, take nominal fit result from different fit")
    parser.add_argument("--toy-json",
                        type=str,
                        default=None,
                        help="Path to the json file containing toy results and p-value")
    parser.add_argument("--legend-position",
                        type=int,
                        default=0,
                        choices=[0,1,2,3],
                        help="Switch between legend positioning in upper or lower half.")
    return parser.parse_args()

def main(args):
    plot.ModTDRStyle(l=0.14, b=0.12)
    # Read in input root file.
    if not os.path.exists(args.input):
        print("Input file {} does not exist.".format(args.input))
        raise ValueError
    infile = ROOT.TFile(args.input, "read")

    # Load the fit results
    if args.separate_nominal_fit is None:
        fit_nominal = infile.Get("fit_nominal")
    else:
        if not os.path.exists(args.separate_nominal_fit):
            print("Input file {} for separate fit does not exist.".format(args.input))
            raise ValueError
        infile_separate = ROOT.TFile(args.separate_nominal_fit, "read")
        fit_nominal = infile_separate.Get("fit_nominal")
    fit_alternate = infile.Get("fit_alternate")

    # Get POI result RooFit object
    res_poi = fit_nominal.floatParsFinal().find(args.poi)

    # Get number of channels from fit result
    prefix = "_ChannelCompatibilityCheck_{}_".format(args.poi)
    parameters = fit_alternate.floatParsFinal()
    num_pars = parameters.getSize()
    channels = [parameters.at(x).GetName() for x in range(num_pars) if parameters.at(x).GetName().startswith(prefix)]
    num_chans = len(channels)
    print(channels)    

    canv = ROOT.TCanvas("canv")
    # Create the plot object.
    axis_label = "#alpha^{H#tau#tau} (degrees)"
    frame = ROOT.TH2F("frame",
                      ";{};".format(axis_label),
                      1,
                      # res_poi.getVal() + 5*res_poi.getAsymErrorLo(),
                      args.frame_range[0],
                      args.frame_range[1],
                      num_chans,
                      0,
                      num_chans
                      )

    print("Nominal fit: {} {}/+{}".format(res_poi.getVal(), res_poi.getAsymErrorLo(), res_poi.getAsymErrorHi()))
    # Fill TGraphAsymmErrors with fit values.
    points = ROOT.TGraphAsymmErrors()

    chan_it = sorted(channels)

    for i, ch_ in enumerate(chan_it):
        ri = parameters.find(ch_)
        p = points.GetN()
        points.SetPoint(p, ri.getVal(), i+0.5)
        points.SetPointError(p, -ri.getAsymErrorLo(), ri.getAsymErrorHi(), 0, 0)
        print("Alternate fit: ", ch_, ri.getVal())

        if ri.getVal()<0:
            points.SetPoint(p+1, ri.getVal()+180, i+0.5)
            points.SetPointError(p+1, -ri.getAsymErrorLo(), ri.getAsymErrorHi(), 0, 0)
        elif ri.getVal()>0:
            points.SetPoint(p+1, ri.getVal()-180, i+0.5)
            points.SetPointError(p+1, -ri.getAsymErrorLo(), ri.getAsymErrorHi(), 0, 0)
        frame.GetYaxis().SetBinLabel(i+1, label_dict[ch_.replace(prefix, "")])

    points.SetLineColor(ROOT.kBlack)
    points.SetLineWidth(3)
    points.SetMarkerStyle(20)
    points.SetMarkerSize(1.3)
    # frame.GetXaxis().SetNdivisions(505)
    frame.GetXaxis().SetTitleSize(0.05)
    frame.GetXaxis().SetLabelSize(0.04)
    if '_pt_' in channels[1]:
        frame.GetYaxis().SetLabelSize(0.04)
    else: frame.GetYaxis().SetLabelSize(0.06)
    frame.Draw()

    # ROOT.gStyle.SetOptStat(0)
    globalFitBand = ROOT.TBox(res_poi.getVal()+res_poi.getAsymErrorLo(), 0, res_poi.getVal()+res_poi.getAsymErrorHi(), num_chans)
    globalFitBand.SetFillStyle(1001)
    # globalFitBand.SetFillColor(17)
    # globalFitBand.SetFillColor(ROOT.TColor.GetColor("#034732"))
    globalFitBand.SetFillColor(ROOT.TColor.GetColor("#73a3d4"))
    globalFitBand.SetLineStyle(0)
    globalFitBand.SetLineColor(ROOT.TColor.GetColor("#73a3d4"))
    globalFitBand.Draw("same")

    globalFitLine = ROOT.TLine(res_poi.getVal(), 0, res_poi.getVal(), num_chans)
    globalFitLine.SetLineWidth(3)
    # globalFitLine.SetLineColor(12)
    globalFitLine.SetLineColor(ROOT.TColor.GetColor("#1d3d5e"))
    globalFitLine.Draw("same")
    points.Draw("0PZ SAME")

    ROOT.gPad.SetTickx()
    ROOT.gPad.SetTicky()
    ROOT.gPad.RedrawAxis()

    globalFitLine_gr = ROOT.TGraph()
    globalFitLine_gr.SetMarkerStyle(20)
    globalFitLine_gr.SetLineWidth(3)

    if args.legend_position == 0:
        #legend = ROOT.TLegend(0.18, 0.77, 0.48, 0.92)
        # legend above plot with 2 columns
        #legend = ROOT.TLegend(0.18, 0.95, 0.9, 1.0)
        #legend.SetNColumns(2)
        legend = ROOT.TLegend(0.17, 0.82, 0.35, 0.92)
        legend.SetFillStyle(0)
        legend.SetNColumns(1)
    elif args.legend_position == 1:
        legend = ROOT.TLegend(0.65, 0.77, 0.95, 0.92)
    elif args.legend_position == 2:
        legend = ROOT.TLegend(0.65, 0.13, 0.95, 0.28)
    else:
        legend = ROOT.TLegend(0.18, 0.13, 0.48, 0.28)

    legend.AddEntry(globalFitLine_gr, "Global Best Fit", "lp")
    legend.AddEntry(globalFitBand, "Global Best Fit #pm 1 #sigma", "f")
    legend.SetTextSize(0.032)
    legend.Draw()

    extra_ts = 0.030
    limit_tree = infile.Get("limit")
    limit_tree.GetEntry(0)
    latex2 = ROOT.TLatex()
    latex2.SetNDC()
    latex2.SetTextAngle(0)
    latex2.SetTextColor(ROOT.kBlack)
    latex2.SetTextSize(extra_ts)
    begin_left = 0.175 #0.200
    begin_right = 0.68
    up_pos = 0.82 #0.690
    spacing = 0.04
    low_pos = 0.470

    if args.toy_json is not None:
        with open(args.toy_json, "r") as fi:
            res = json.load(fi)["125.0"]
        label = "#font[42]{{p-value = {:.2f}}}".format(res["p"])
        if '_pt_' in channels[1]:
            latex2.DrawLatex(0.26, up_pos-1*spacing, label)
        elif args.legend_position == 0:
            latex2.DrawLatex(begin_left, up_pos-1*spacing, label)
        elif args.legend_position == 1:
            latex2.DrawLatex(begin_right, up_pos-1*spacing, label)
        elif args.legend_position == 2:
            latex2.DrawLatex(begin_right, low_pos-1*spacing, label)
        elif args.legend_position == 3:
            latex2.DrawLatex(begin_left, low_pos-1*spacing, label)

    # Draw CMS logo in upper left corner
    cms_outside = True
    if cms_outside:
        #plot.DrawCMSLogo(canv, 'CMS', 'Supplementary', 0,
        #                   0.095, 0.05, 1.0, '', 0.6)
        plot.DrawCMSLogo(canv, 'CMS', 'Supplementary', 0,
                           0.14, 0.05, 2.0, '', 0.85)
    else:
        plot.DrawCMSLogo(canv, 'CMS', 'Supplementary', 11,
                           0.045, 0.05, 1.0, '', 0.6)
    # Draw luminosity label
    #styles.DrawTitle(canv, "62.4 fb^{-1} (13.6 TeV)", 3, 0.6)

    canv.Print(args.output + ".pdf", "pdf")
    return

if __name__ == "__main__":
    args = parse_args()
    main(args)
