import CombineHarvester.CombineTools.ch as ch
import yaml
# Note CMS common systematics should be named following: https://gitlab.cern.ch/cms-analysis/general/systematics/-/blob/master/systematics_master.yml?ref_type=heads, analysis specific ones should eventually have the CADI number in the names

def AddSMRun3Systematics(cb, YR18_uncertainties=False):

    if YR18_uncertainties:
        # for YR18 scheme for HL-scaling most experimental uncertainties are scaled by the 1/sqrt(L), except for a few that are scaled down to floor values
        # theory uncertainties are scale by 50%
        #YR18_exp_scale = (62.4/3000)**.5
        #YR18_theory_scale = 0.5
        YR18_exp_scale = 1.0
        YR18_theory_scale = 1.0  

    else: 
        YR18_exp_scale = 1.0
        YR18_theory_scale = 1.0

    cats_tt = {
        1:'tt_mva_tau',
        2:'tt_mva_fake',
        3:'tt_higgs_rhorho',
        4:'tt_higgs_rhoa11pr',
        5:'tt_higgs_rhoa1',
        6:'tt_higgs_a1a1',
        7:'tt_higgs_pirho',
        8:'tt_higgs_pipi',
        9:'tt_higgs_pia1',
        10: 'tt_higgs_pia11pr',
        11: 'tt_higgs_a11pra1',
        }

    eras = ['2022', '2022EE', '2023', '2023BPix']
    
    # define processes lists
    sig_procs = {}
    sig_procs['ggH'] = ['ggH_sm_prod_sm_htt','ggH_ps_prod_sm_htt','ggH_mm_prod_sm_htt', 'ggH_sm_prod_ps_htt','ggH_ps_prod_ps_htt','ggH_mm_prod_ps_htt', 'ggH_sm_prod_mm_htt','ggH_ps_prod_mm_htt','ggH_mm_prod_mm_htt']
    sig_procs['VBF'] = ['qqH_sm_htt','qqH_ps_htt','qqH_mm_htt']
    sig_procs['ZH'] = ['ZH_sm_htt','ZH_ps_htt','ZH_mm_htt']
    sig_procs['WH'] = ['WH_sm_htt','WH_ps_htt','WH_mm_htt']
    # sig_procs['qqH'] = ['qqH_sm_htt','qqH_ps_htt','qqH_mm_htt','WH_sm_htt','WH_ps_htt','WH_mm_htt','ZH_sm_htt','ZH_ps_htt','ZH_mm_htt']
    
    dy_procs = ['ZTT', 'ZL']
    ttbar_procs = ['TTT']
    vv_procs = ['VVT']
    bkg_mc_procs = dy_procs + ttbar_procs + vv_procs #+ ['JetFakesSublead']
    
    mc_procs = bkg_mc_procs
    for p in sig_procs.values(): mc_procs+=p

    recoil_procs = ['ZTT','ZL']
    for p in sig_procs.values(): recoil_procs+=p
   
    ###############################################
    # Luminosity
    ###############################################

    
    # lumi uncertainty from here: https://cms-talk.web.cern.ch/t/luminosity-uncertainty-correlations-between-run-2-and-2022-and-2023/132007
    if not YR18_uncertainties: cb.cp().process(mc_procs).AddSyst(cb, 'lumi_13p6TeV_2223', 'lnN', ch.SystMap()(1.0102))
    else: cb.cp().process(mc_procs).AddSyst(cb, 'lumi', 'lnN', ch.SystMap()(1.01))

    ###############################################
    # Pileup
    ###############################################
    cb.cp().process(mc_procs).AddSyst(cb, 'CMS_pileup', "shape", ch.SystMap()(YR18_exp_scale))

    ###############################################
    # Cross sections and BRs
    ###############################################

    # Cross-sections uncertainties - we keep naming consistent with Run-2 analysis for now

    # DY XS uncertainties from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MATRIXCrossSectionsat13p6TeV
    # Quadrature sum of scale, PDF, and difference between NLO additive vs multiplicative
    cb.cp().process(dy_procs).AddSyst(cb, 'cross_section_Z', 'lnN', ch.SystMap()(((1-0.984)*YR18_theory_scale + 1,(1-1.013)*YR18_theory_scale + 1)))
    
    # ttbar cross-section uncertainties from here: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
    # Quadrature sum of scale, PDF, and mass uncerts
    cb.cp().process(ttbar_procs).AddSyst(cb, 'cross_section_ttbar', 'lnN', ch.SystMap()(((1-0.949)*YR18_theory_scale + 1,(1-1.044)*YR18_theory_scale + 1)))

    #For VV in principle can take NNLO numbers from here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/MATRIXCrossSectionsat13p6TeV
    #But since we mix together all VV + rare procs into this template, we take a conservative 5% uncertainty (same as for Run-2)
    cb.cp().process(vv_procs).AddSyst(cb, 'cross_section_VV', 'lnN', ch.SystMap()((1-1.05)*YR18_theory_scale + 1))

    # Higgs cross-section uncertainties from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHWG136TeVxsec_extrap

    # QCD scale uncertainties
    cb.cp().process(sig_procs['ggH']).AddSyst(cb, 'QCDscale_ggH', 'lnN', ch.SystMap()(((1-1.039)*YR18_theory_scale + 1)))
    
    cb.cp().process(sig_procs['WH']).AddSyst(cb, 'QCDscale_VH', 'lnN', ch.SystMap()(((1-0.993)*YR18_theory_scale + 1,(1-1.004)*YR18_theory_scale + 1)))
    
    cb.cp().process(sig_procs['ZH']).AddSyst(cb, 'QCDscale_VH', 'lnN', ch.SystMap()(((1-0.968)*YR18_theory_scale + 1,(1-1.038)*YR18_theory_scale + 1)))
    
    cb.cp().process(sig_procs['VBF']).AddSyst(cb, 'QCDscale_VH', 'lnN', ch.SystMap()(((1-0.997)*YR18_theory_scale + 1,(1-1.005)*YR18_theory_scale + 1)))

    # PDF uncertainties
    cb.cp().process(sig_procs['ggH']).AddSyst(cb, 'pdf_Higgs_gg', 'lnN', ch.SystMap()((1-1.032)*YR18_theory_scale + 1))
    
    cb.cp().process(sig_procs['VBF']).AddSyst(cb, 'pdf_Higgs_qqbar', 'lnN', ch.SystMap()((1-1.032)*YR18_theory_scale + 1))
    
    cb.cp().process(sig_procs['WH']).AddSyst(cb, 'pdf_Higgs_qqbar', 'lnN', ch.SystMap()((1-1.016)*YR18_theory_scale + 1))
    
    cb.cp().process(sig_procs['ZH']).AddSyst(cb, 'pdf_Higgs_qqbar', 'lnN', ch.SystMap()((1-1.013)*YR18_theory_scale + 1))
    
    # H->tautau BR uncertainties from: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR#SM_Higgs_Branching_Ratios_and_To
    cb.cp().process(sig_procs['ggH']+sig_procs['VBF']+sig_procs['WH']+sig_procs['ZH']).AddSyst(cb, 'BR_htt_THU', 'lnN', ch.SystMap()(((1-0.984)*YR18_theory_scale + 1,(1-1.017)*YR18_theory_scale + 1)))
    cb.cp().process(sig_procs['ggH']+sig_procs['VBF']+sig_procs['WH']+sig_procs['ZH']).AddSyst(cb, 'BR_htt_PU_mq', 'lnN', ch.SystMap()((1-1.010)*YR18_theory_scale + 1))
    cb.cp().process(sig_procs['ggH']+sig_procs['VBF']+sig_procs['WH']+sig_procs['ZH']).AddSyst(cb, 'BR_htt_PU_alphas', 'lnN', ch.SystMap()((1-1.006)*YR18_theory_scale + 1))
    
    
    ###############################################
    # Shape/acceptance theory uncertainties
    ###############################################

    # signal theory uncertainties

    # QCD scale variations
    cb.cp().process(sig_procs['ggH']+sig_procs['VBF']).AddSyst(cb, "QCDscale_ren_signal", "shape", ch.SystMap()(1.0))
    cb.cp().process(sig_procs['ggH']+sig_procs['VBF']).AddSyst(cb, "QCDscale_fac_signal", "shape", ch.SystMap()(1.0)) # note these 2 will be scaled down in harvesting code instead

    # parton shower variations
    cb.cp().process(sig_procs['ggH']+sig_procs['VBF']).AddSyst(cb, "ps_isr_signal", "shape", ch.SystMap()(YR18_theory_scale))
    cb.cp().process(sig_procs['ggH']+sig_procs['VBF']).AddSyst(cb, "ps_fsr_signal", "shape", ch.SystMap()(YR18_theory_scale))

    # DY shape uncertainty (e.g from pT/mass reweighting)
    cb.cp().process(dy_procs).AddSyst(cb, "CMS_HIG25012_Z_pt_reweighting", "shape", ch.SystMap()(YR18_exp_scale))

    # ttbar shape uncertainty (e.g from pT reweighting)
    cb.cp().process(ttbar_procs).AddSyst(cb, "top_pt_reweighting", "shape", ch.SystMap()(1./3))
    
    ###############################################
    # Offline object identification
    ###############################################

    # muon ID and isolation
    if YR18_uncertainties:
        cb.cp().process(mc_procs).channel(['mt']).AddSyst(cb, "CMS_eff_m", "lnN", ch.SystMap()(1.005))
    else:
        cb.cp().process(mc_procs).channel(['mt']).AddSyst(cb, "CMS_eff_m_id", "shape", ch.SystMap()(1.0))
        cb.cp().process(mc_procs).channel(['mt']).AddSyst(cb, "CMS_eff_m_iso", "shape", ch.SystMap()(1.0))

    # electron ID and reco
    if YR18_uncertainties:
        cb.cp().process(mc_procs).channel(['et']).AddSyst(cb, "CMS_eff_e", "lnN", ch.SystMap()(1.005))
    else: 
        cb.cp().process(mc_procs).channel(['et']).AddSyst(cb, "CMS_eff_e_reco", "shape", ch.SystMap()(1.0))
        cb.cp().process(mc_procs).channel(['et']).AddSyst(cb, "CMS_eff_e_id", "shape", ch.SystMap()(1.0))

    # electron scale and smearing
    cb.cp().process(mc_procs).channel(['et']).AddSyst(cb, "CMS_scale_e", "shape", ch.SystMap()(YR18_exp_scale))

    # mu->tau fakes (ZL only in mt)
    for era in eras:
        for eta in ['0p0','0p4','0p8','1p2','1p7']:
            cb.cp().process(mc_procs).process(['ZL']).channel(['mt']).AddSyst(cb, f'CMS_fake_t_DeepTau2018v2p5_VSmu_{era}_eta_{eta}', 'shape', ch.SystMap()(YR18_exp_scale))
            # decided not to increase or decouple

    # e->tau fakes (ZL only in et)
    for era in eras:
        for eta in ['0p0', '1p5']:
            for dm in ['0', '1', '2', '10']:
                cb.cp().process(mc_procs).process(['ZL']).bin_id([1,2]).channel(['et']).AddSyst(cb, f'CMS_fake_t_DeepTau2018v2p5_VSe_{era}_eta_{eta}_DM{dm}PNet', 'shape', ch.SystMap()(YR18_exp_scale))

            cb.cp().process(mc_procs).process(['ZL']).bin_id([4]).channel(['et']).AddSyst(cb, f'CMS_fake_t_DeepTau2018v2p5_VSe_{era}_eta_{eta}_DM0PNet', 'shape', ch.SystMap()(YR18_exp_scale))
            cb.cp().process(mc_procs).process(['ZL']).bin_id([3]).channel(['et']).AddSyst(cb, f'CMS_fake_t_DeepTau2018v2p5_VSe_{era}_eta_{eta}_DM1PNet', 'shape', ch.SystMap()(YR18_exp_scale))
            cb.cp().process(mc_procs).process(['ZL']).bin_id([6]).channel(['et']).AddSyst(cb, f'CMS_fake_t_DeepTau2018v2p5_VSe_{era}_eta_{eta}_DM2PNet', 'shape', ch.SystMap()(YR18_exp_scale))
            cb.cp().process(mc_procs).process(['ZL']).bin_id([5]).channel(['et']).AddSyst(cb, f'CMS_fake_t_DeepTau2018v2p5_VSe_{era}_eta_{eta}_DM10PNet', 'shape', ch.SystMap()(YR18_exp_scale))

    # Genuine Tau ID

    for era in eras:
        # statistical uncertainties from fitted function parameters
        for u in ['stat1','stat2']:
            for dm in ['0', '1', '2', '10']:
                cb.cp().process(mc_procs).process(['ZL'], False).bin_id([1,2]).AddSyst(cb, f'CMS_HIG25012_eff_t_DeepTau2018v2p5_VSjet_{u}_DM{dm}PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
                # uncomment the line below if want ZL to have Tau ID unct. as well
                # cb.cp().process(['ZL']).channel(['tt']).bin_id([1,2]).AddSyst(cb, f'CMS_HIG25012_eff_t_DeepTau2018v2p5_VSjet_{u}_DM{dm}PNet_{era}', 'shape', ch.SystMap()(1.0))

            # Add DM specific uncertainties for other bins (tt channels)
            cb.cp().process(mc_procs).process(['ZL'], False).channel(['tt']).bin_id([7,8,9,10]).AddSyst(cb, f'CMS_HIG25012_eff_t_DeepTau2018v2p5_VSjet_{u}_DM0PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
            cb.cp().process(mc_procs).process(['ZL'], False).channel(['tt']).bin_id([3,4,5,7]).AddSyst(cb, f'CMS_HIG25012_eff_t_DeepTau2018v2p5_VSjet_{u}_DM1PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
            cb.cp().process(mc_procs).process(['ZL'], False).channel(['tt']).bin_id([4,10,11]).AddSyst(cb, f'CMS_HIG25012_eff_t_DeepTau2018v2p5_VSjet_{u}_DM2PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
            cb.cp().process(mc_procs).process(['ZL'], False).channel(['tt']).bin_id([5,6,9,11]).AddSyst(cb, f'CMS_HIG25012_eff_t_DeepTau2018v2p5_VSjet_{u}_DM10PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
            # Add DM specific uncertainties for other bins (lt channels)
            cb.cp().process(mc_procs).process(['ZL'], False).channel(['tt'], False).bin_id([4]).AddSyst(cb, f'CMS_HIG25012_eff_t_DeepTau2018v2p5_VSjet_{u}_DM0PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
            cb.cp().process(mc_procs).process(['ZL'], False).channel(['tt'], False).bin_id([3]).AddSyst(cb, f'CMS_HIG25012_eff_t_DeepTau2018v2p5_VSjet_{u}_DM1PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
            cb.cp().process(mc_procs).process(['ZL'], False).channel(['tt'], False).bin_id([6]).AddSyst(cb, f'CMS_HIG25012_eff_t_DeepTau2018v2p5_VSjet_{u}_DM2PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
            cb.cp().process(mc_procs).process(['ZL'], False).channel(['tt'], False).bin_id([5]).AddSyst(cb, f'CMS_HIG25012_eff_t_DeepTau2018v2p5_VSjet_{u}_DM10PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))


        # systematic that is correlated across decay modes but decorrelated across eras
        cb.cp().process(mc_procs).process(['ZL'], False).AddSyst(cb, f'CMS_HIG25012_eff_t_DeepTau2018v2p5_VSjet_syst_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
    # systematic that is correlated across eras and decay modes
    cb.cp().process(mc_procs).process(['ZL'], False).AddSyst(cb, f'CMS_HIG25012_eff_t_DeepTau2018v2p5_VSjet_syst_alleras', 'shape', ch.SystMap()(YR18_exp_scale))

    if YR18_uncertainties:
        # add floor values of 2.5% per tau
        cb.cp().process(mc_procs).channel(['tt'],False).process(['ZL'], False).AddSyst(cb, f'eff_t_syst_floor', 'lnN', ch.SystMap()(1.025))
        cb.cp().process(mc_procs).channel(['tt']).process(['ZL'], False).AddSyst(cb, f'eff_t_syst_floor', 'lnN', ch.SystMap()(1.051))
   
    ###############################################
    # Trigger
    ###############################################

    # muon trigger
    if YR18_uncertainties:
        cb.cp().process(mc_procs).channel(['mt']).AddSyst(cb, f'CMS_HIG25012_eff_m_trigger', 'lnN', ch.SystMap()(1.005))
    else:
        cb.cp().process(mc_procs).channel(['mt']).AddSyst(cb, f'CMS_HIG25012_eff_m_trigger', 'shape', ch.SystMap()(1.0))
    
    # electron trigger
    if YR18_uncertainties:
        cb.cp().process(mc_procs).channel(['et']).AddSyst(cb, f'CMS_HIG25012_eff_e_trigger', 'lnN', ch.SystMap()(1.005))
    else:
        cb.cp().process(mc_procs).channel(['et']).AddSyst(cb, f'CMS_HIG25012_eff_e_trigger', 'shape', ch.SystMap()(1.0))

    # tau trigger
    # statistical uncertainties
    for era in eras:
        # tau leg uncertainties
        for trig in ['ditau','ditaujet']:
            for dm in ['0', '1', '2', '10']:
                cb.cp().process(mc_procs).channel(['tt']).bin_id([1,2]).AddSyst(cb, f'CMS_HIG25012_trig_t_{trig}_VTight_DM{dm}PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
            cb.cp().process(mc_procs).channel(['tt']).bin_id([7,8,9,10]).AddSyst(cb, f'CMS_HIG25012_trig_t_{trig}_VTight_DM0PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
            cb.cp().process(mc_procs).channel(['tt']).bin_id([3,4,5,7]).AddSyst(cb, f'CMS_HIG25012_trig_t_{trig}_VTight_DM1PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
            cb.cp().process(mc_procs).channel(['tt']).bin_id([4,10,11]).AddSyst(cb, f'CMS_HIG25012_trig_t_{trig}_VTight_DM2PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
            cb.cp().process(mc_procs).channel(['tt']).bin_id([5,6,9,11]).AddSyst(cb, f'CMS_HIG25012_trig_t_{trig}_VTight_DM10PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))

        # jet leg uncertainties
        cb.cp().process(mc_procs).channel(['tt']).AddSyst(cb, f'CMS_HIG25012_trig_j_ditaujet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))

    # We also add a 3% systematic uncertainty due to the background modelling in the SF extraction based on the studies in https://indico.cern.ch/event/1263107/contributions/5306043/attachments/2606862/4503028/tautriggerSF_checks.pdf
    cb.cp().process(mc_procs).channel(['tt']).AddSyst(cb, "CMS_HIG25012_trig_t_ditau_syst", "lnN", ch.SystMap()((1-1.03)*YR18_exp_scale + 1))
    
    ###############################################
    # Lepton/Tau energy scales
    ###############################################



    # e->tau fake energy scale
    for era in eras:
        for dm in ['0', '1', '2', '10']:
            cb.cp().process(mc_procs).process(['ZL']).channel(['et']).bin_id([1,2]).AddSyst(cb, f'CMS_HIG25012_scale_t_eFake_DM{dm}PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))

        cb.cp().process(mc_procs).process(['ZL']).channel(['et']).bin_id([4]).AddSyst(cb, f'CMS_HIG25012_scale_t_eFake_DM0PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
        cb.cp().process(mc_procs).process(['ZL']).channel(['et']).bin_id([3]).AddSyst(cb, f'CMS_HIG25012_scale_t_eFake_DM1PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
        cb.cp().process(mc_procs).process(['ZL']).channel(['et']).bin_id([6]).AddSyst(cb, f'CMS_HIG25012_scale_t_eFake_DM2PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
        cb.cp().process(mc_procs).process(['ZL']).channel(['et']).bin_id([5]).AddSyst(cb, f'CMS_HIG25012_scale_t_eFake_DM10PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))


    # mu->tau fake energy scale
    for era in eras:
        for dm in ['0', '1', '2', '10']:
            cb.cp().process(mc_procs).process(['ZL']).channel(['mt']).bin_id([1,2]).AddSyst(cb, f'CMS_HIG25012_scale_t_muFake_DM{dm}PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))

        cb.cp().process(mc_procs).process(['ZL']).channel(['mt']).bin_id([4]).AddSyst(cb, f'CMS_HIG25012_scale_t_muFake_DM0PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
        cb.cp().process(mc_procs).process(['ZL']).channel(['mt']).bin_id([3]).AddSyst(cb, f'CMS_HIG25012_scale_t_muFake_DM1PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
        cb.cp().process(mc_procs).process(['ZL']).channel(['mt']).bin_id([6]).AddSyst(cb, f'CMS_HIG25012_scale_t_muFake_DM2PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
        cb.cp().process(mc_procs).process(['ZL']).channel(['mt']).bin_id([5]).AddSyst(cb, f'CMS_HIG25012_scale_t_muFake_DM10PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))

    # genuine tau energy scale
    for era in eras:
        for dm in ['0', '1', '2', '10']:
            cb.cp().process(mc_procs).process(['ZL'], False).bin_id([1,2]).AddSyst(cb, f'CMS_HIG25012_scale_t_DM{dm}PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))

        cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt']).bin_id([7,8,9,10]).AddSyst(cb, f'CMS_HIG25012_scale_t_DM0PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
        cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt']).bin_id([3,4,5,7]).AddSyst(cb, f'CMS_HIG25012_scale_t_DM1PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
        cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt']).bin_id([4,10,11]).AddSyst(cb, f'CMS_HIG25012_scale_t_DM2PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
        cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt']).bin_id([5,6,9,11]).AddSyst(cb, f'CMS_HIG25012_scale_t_DM10PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))

        cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt'], False).bin_id([4]).AddSyst(cb, f'CMS_HIG25012_scale_t_DM0PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
        cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt'], False).bin_id([3]).AddSyst(cb, f'CMS_HIG25012_scale_t_DM1PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
        cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt'], False).bin_id([6]).AddSyst(cb, f'CMS_HIG25012_scale_t_DM2PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))
        cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt'], False).bin_id([5]).AddSyst(cb, f'CMS_HIG25012_scale_t_DM10PNet_{era}', 'shape', ch.SystMap()(YR18_exp_scale))

    ###############################################
    # IP significance cut
    ###############################################
    for tautype in ['prompt', 'tauDecay']:
        for eta in ['Lt1p0', '1p0to1p6', 'Gt1p6']:
            cb.cp().process(mc_procs).channel(['tt'], False).AddSyst(cb, f'CMS_HIG25012_eff_IPSigCut_{tautype}_eta_{eta}', 'shape', ch.SystMap()(YR18_exp_scale))


    ###############################################
    # Jet/MET scale/resolutions
    ###############################################

    # Combined
    # cb.cp().process(mc_procs).AddSyst(cb, 'CMS_scale_j_13p6TeV', 'shape', ch.SystMap()(1.0))
    # cb.cp().process(mc_procs).AddSyst(cb, 'CMS_res_j_13p6TeV', 'shape', ch.SystMap()(1.0))

    # Regrouped JEC
    jec_variations_correlated = ['Absolute', 'BBEC1', 'EC2', 'FlavorQCD', 'HF', 'RelativeBal']
    for var in jec_variations_correlated:
        cb.cp().process(mc_procs).AddSyst(cb, f'CMS_scale_j_{var}', 'shape', ch.SystMap()(0.5 if YR18_uncertainties else 1.0))

    for era in eras:
        # uncorrelated JEC uncertainties
        jec_variations_uncorrelated = [f'Absolute_{era}', f'BBEC1_{era}', f'EC2_{era}', f'HF_{era}', f'RelativeSample_{era}']
        for var in jec_variations_uncorrelated:
            cb.cp().process(mc_procs).AddSyst(cb, f'CMS_scale_j_{var}', 'shape', ch.SystMap()(0.5 if YR18_uncertainties else 1.0))
        # uncorrelated JER
        cb.cp().process(mc_procs).AddSyst(cb, f'CMS_res_j_{era}', 'shape', ch.SystMap()(0.5 if YR18_uncertainties else 1.0))

    # JER uncorrelated between eras
    

    # TODO: MET uncl
    
    # MET recoil uncertainties

    cb.cp().process(recoil_procs).AddSyst(cb,'CMS_HIG25012_scale_met', 'shape', ch.SystMap()(0.5 if YR18_uncertainties else 1.0))
    cb.cp().process(recoil_procs).AddSyst(cb,'CMS_HIG25012_res_met', 'shape', ch.SystMap()(0.5 if YR18_uncertainties else 1.0))
    
    ###############################################
    # jet->tau fake-factors
    ###############################################

    # tt channel statistical uncertainties
    for njets in [0,1,2]:
        cb.cp().process(['JetFakes']).channel(['tt']).bin_id([1,2,7,8,9,10]).AddSyst(cb, f"CMS_HIG25012_fake_t_stat_dm0_{njets}j", "shape", ch.SystMap()(YR18_exp_scale))
        cb.cp().process(['JetFakes']).channel(['tt']).bin_id([1,2,3,4,5,7]).AddSyst(cb, f"CMS_HIG25012_fake_t_stat_dm1_{njets}j", "shape", ch.SystMap()(YR18_exp_scale))
        cb.cp().process(['JetFakes']).channel(['tt']).bin_id([1,2,4,10,11]).AddSyst(cb, f"CMS_HIG25012_fake_t_stat_dm2_{njets}j", "shape", ch.SystMap()(YR18_exp_scale))
        cb.cp().process(['JetFakes']).channel(['tt']).bin_id([1,2,5,6,9,11]).AddSyst(cb, f"CMS_HIG25012_fake_t_stat_dm10_{njets}j", "shape", ch.SystMap()(YR18_exp_scale))

    # add lnN FF uncertainty from yml file located in configs/ff_lnN_uncertainties.yml
    # open yml file and read uncertainties
    with open('configs/ff_lnN_uncertainties.yml', 'r') as f:
        ff_lnN_uncertainties = yaml.safe_load(f)
        cb.cp().process(['JetFakes']).channel(['tt']).AddSyst(cb, "CMS_HIG25012_fake_t_syst_lnN", "lnN", ch.SystMap()((1-ff_lnN_uncertainties['tt']['correlated'])*YR18_exp_scale + 1))
        for i in range(1, 12):
            # add lnN uncertainty for each decay mode
            cb.cp().process(['JetFakes']).channel(['tt']).bin_id([i]).AddSyst(cb, "CMS_HIG25012_fake_t_syst_lnN_$BIN", "lnN", ch.SystMap()((1-ff_lnN_uncertainties['tt'][cats_tt[i]])*YR18_exp_scale + 1))

    # add shape uncertainties for BDT score, this is decorrelated between tau, fake, and higgs categories
    cb.cp().process(['JetFakes']).channel(['tt']).bin_id([1]).AddSyst(cb, "CMS_HIG25012_fake_t_syst_BDTshape_tau", "shape", ch.SystMap()(YR18_exp_scale))
    cb.cp().process(['JetFakes']).channel(['tt']).bin_id([2]).AddSyst(cb, "CMS_HIG25012_fake_t_syst_BDTshape_fake", "shape", ch.SystMap()(YR18_exp_scale))
    cb.cp().process(['JetFakes']).channel(['tt']).bin_id([1,2], False).AddSyst(cb, "CMS_HIG25012_fake_t_syst_BDTshape_higgs", "shape", ch.SystMap()(YR18_exp_scale))

    # add shape uncertainties for aco angle, this is decorrelated between categories
    cb.cp().process(['JetFakes']).channel(['tt']).bin_id([1,2], False).AddSyst(cb, "CMS_HIG25012_fake_t_syst_acoshape_$BIN", "shape", ch.SystMap()(YR18_exp_scale))
    # uncertainty due to subtracted real taus in the tt channel (uncorrelated in the harvester)
    cb.cp().process(['JetFakes']).channel(['tt']).AddSyst(cb, "CMS_HIG25012_fake_t_sub_syst", "shape", ch.SystMap()(YR18_exp_scale))
    # lnN uncertainty for the JetFakesSublead as it is estimated from MC
    cb.cp().process(['JetFakesSublead']).channel(['tt']).AddSyst(cb, "CMS_HIG25012_fake_t_mc", "lnN", ch.SystMap()((1-1.3)*YR18_exp_scale + 1))

    # mt and et channel FF uncertainties (uncorrelated in the harvester)
    for ff_type in ['wj','qcd','mc_top']:
        # statistics
        for njets in [0,1,2]:
            cb.cp().process(['JetFakes']).channel(['mt','et']).bin_id([1,2,4]).AddSyst(cb, f"CMS_HIG25012_fake_t_{ff_type}_stat_dm0_{njets}j", "shape", ch.SystMap()(YR18_exp_scale))
            cb.cp().process(['JetFakes']).channel(['mt','et']).bin_id([1,2,3]).AddSyst(cb, f"CMS_HIG25012_fake_t_{ff_type}_stat_dm1_{njets}j", "shape", ch.SystMap()(YR18_exp_scale))
            cb.cp().process(['JetFakes']).channel(['mt','et']).bin_id([1,2,6]).AddSyst(cb, f"CMS_HIG25012_fake_t_{ff_type}_stat_dm2_{njets}j", "shape", ch.SystMap()(YR18_exp_scale))
            cb.cp().process(['JetFakes']).channel(['mt','et']).bin_id([1,2,5]).AddSyst(cb, f"CMS_HIG25012_fake_t_{ff_type}_stat_dm10_{njets}j", "shape", ch.SystMap()(YR18_exp_scale))

        # systematic # uncorrelate et and mt for qcd only
        cb.cp().process(['JetFakes']).channel(['mt','et']).AddSyst(cb, f"CMS_HIG25012_fake_t_{ff_type}_syst", "shape", ch.SystMap()(YR18_exp_scale))

    # uncertainty due to subtracted real taus in the mt and et channels (uncorrelated in the harvester)
    cb.cp().process(['JetFakes']).channel(['mt','et']).AddSyst(cb, "CMS_HIG25012_fake_t_sub_syst", "shape", ch.SystMap()(YR18_exp_scale))

    ###############################################
    # 4-vectors for CP angle reconstruction
    ###############################################

    # IP direction/scale
    cb.cp().process(mc_procs).channel(['tt']).bin_id([7,8,9,10]).AddSyst(cb, "CMS_HIG25012_res_IP", "shape", ch.SystMap()(YR18_exp_scale))
    cb.cp().process(mc_procs).channel(['tt'], False).bin_id([3,4,5,6]).AddSyst(cb, "CMS_HIG25012_res_IP", "shape", ch.SystMap()(YR18_exp_scale))
    
    # TODO: pi0 direction/scale (not included for Run-2 but could add)
    
    # TODO: pi direction/scale (not included for Run-2 but could add)
    
    # SV vertex resolution uncertainty
    cb.cp().process(mc_procs).channel(['tt']).bin_id([5,6,9,11]).AddSyst(cb, "CMS_HIG25012_res_SV", "shape", ch.SystMap()(YR18_exp_scale))
    cb.cp().process(mc_procs).channel(['tt'], False).bin_id([5]).AddSyst(cb, "CMS_HIG25012_res_SV", "shape", ch.SystMap()(YR18_exp_scale))


    ###############################################
    # DM-migration uncertainties
    ###############################################

    cb.cp().process(mc_procs).process(['ZL'],False).channel(['mt','et']).bin_id([4]).AddSyst(cb, "CMS_HIG25012_DM_migrations_GenDM0", "shape", ch.SystMap()(YR18_exp_scale))
    cb.cp().process(mc_procs).process(['ZL'],False).channel(['mt','et']).bin_id([3]).AddSyst(cb, "CMS_HIG25012_DM_migrations_GenDM1", "shape", ch.SystMap()(YR18_exp_scale))
    cb.cp().process(mc_procs).process(['ZL'],False).channel(['mt','et']).bin_id([6]).AddSyst(cb, "CMS_HIG25012_DM_migrations_GenDM2", "shape", ch.SystMap()(YR18_exp_scale))
    cb.cp().process(mc_procs).process(['ZL'],False).channel(['mt','et']).bin_id([5]).AddSyst(cb, "CMS_HIG25012_DM_migrations_GenDM10", "shape", ch.SystMap()(YR18_exp_scale))

    cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt']).bin_id([7,8,9,10]).AddSyst(cb, "CMS_HIG25012_DM_migrations_GenDM0", "shape", ch.SystMap()(YR18_exp_scale))
    cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt']).bin_id([3,4,5,7]).AddSyst(cb, "CMS_HIG25012_DM_migrations_GenDM1", "shape", ch.SystMap()(YR18_exp_scale))
    cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt']).bin_id([4,10,11]).AddSyst(cb, "CMS_HIG25012_DM_migrations_GenDM2", "shape", ch.SystMap()(YR18_exp_scale))
    cb.cp().process(mc_procs).process(['ZL'],False).channel(['tt']).bin_id([5,6,9,11]).AddSyst(cb, "CMS_HIG25012_DM_migrations_GenDM10", "shape", ch.SystMap()(YR18_exp_scale))

    return cb
