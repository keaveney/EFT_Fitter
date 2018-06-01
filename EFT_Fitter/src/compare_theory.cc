#include "compare_theory.h"

int main(){
    
    make_plot("p_{T}^{t}","GeV","data/abs/DiffXS_HypToppT_source.root","abs", preds_top_pt, pred_names_top_pt,divide_options_top_pt, plot_options_top_pt);
    
    make_plot("p_{T}^{t}","GeV","data/norm/DiffXS_HypToppT_source.root","norm", preds_top_pt, pred_names_top_pt,divide_options_top_pt, plot_options_top_pt_norm);
    
    make_plot("p_{T}^{#bar{t}}","GeV","data/abs/DiffXS_HypAntiToppT_source.root","abs", preds_antitop_pt, pred_names_antitop_pt,divide_options_antitop_pt, plot_options_top_pt);
    
    make_plot("p_{T}^{#bar{t}}","GeV","data/norm/DiffXS_HypAntiToppT_source.root","norm", preds_antitop_pt, pred_names_antitop_pt,divide_options_antitop_pt, plot_options_top_pt_norm);
    
    make_plot("y_{t}","","data/abs/DiffXS_HypTopRapidity_source.root","abs", preds_top_y, pred_names_top_y,divide_options_top_y, plot_options_top_y);
    make_plot("y_{t}","","data/norm/DiffXS_HypTopRapidity_source.root","norm", preds_top_y, pred_names_top_y,divide_options_top_y, plot_options_top_y_norm);
    
    make_plot("y_{#bar{t}}","","data/abs/DiffXS_HypAntiTopRapidity_source.root","abs", preds_antitop_y, pred_names_antitop_y,divide_options_antitop_y, plot_options_antitop_y);
    
    make_plot("y_{#bar{t}}","","data/norm/DiffXS_HypAntiTopRapidity_source.root","norm", preds_antitop_y, pred_names_antitop_y,divide_options_antitop_y, plot_options_antitop_y_norm);
    
    make_plot("m_{tt}","GeV","data/abs/DiffXS_HypTTBarMass_source.root","abs", preds_mtt, pred_names_mtt,divide_options_mtt, plot_options_mtt);
    
    make_plot("m_{tt}","GeV","data/norm/DiffXS_HypTTBarMass_source.root","norm", preds_mtt, pred_names_mtt,divide_options_mtt, plot_options_mtt_norm);

    make_plot("p^{tt}_{T}","GeV","data/abs/DiffXS_HypTTBarpT_source.root","abs", preds_pttt, pred_names_pttt,divide_options_pttt, plot_options_pttt);
    
    make_plot("p^{tt}_{T}","GeV","data/norm/DiffXS_HypTTBarpT_source.root","norm", preds_pttt, pred_names_pttt,divide_options_pttt, plot_options_pttt_norm);

    make_plot("y_{tt}","","data/abs/DiffXS_HypTTBarRapidity_source.root","abs", preds_ytt, pred_names_ytt,divide_options_ytt, plot_options_ytt);
    make_plot("y_{tt}","","data/norm/DiffXS_HypTTBarRapidity_source.root","norm", preds_ytt, pred_names_ytt,divide_options_ytt, plot_options_ytt_norm);
    
    make_plot("#Delta|y|(t, #bar{t})","","data/abs/DiffXS_HypTTBarDeltaRapidity_source.root ","abs", preds_dytt, pred_names_dytt,divide_options_dytt, plot_options_dytt);
    make_plot("#Delta|y|(t, #bar{t})","","data/norm/DiffXS_HypTTBarDeltaRapidity_source.root ","norm", preds_dytt, pred_names_dytt,divide_options_dytt, plot_options_dytt_norm);
    
    return 1;
}

