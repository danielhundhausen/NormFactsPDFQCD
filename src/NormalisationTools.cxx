#include "UHH2/NormFactsQCDPDF/include/NormalisationTools.h"
#include "UHH2/core/include/Event.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>

using namespace std;
using namespace uhh2;

NormalisationHists::NormalisationHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
  init(ctx); 
}

void NormalisationHists::init(Context & ctx){
  cout << "... Initialize PreselectionHists ..." << endl;

  take_ntupleweights = !(m_oname.Contains("QCD"));
  TString m_pdfname = "NNPDF31_nnlo_hessian_pdfas";

  is_dy = ctx.get("dataset_version").find("DYJets") == 0;
  is_wjets = ctx.get("dataset_version").find("WJets") == 0;
  is_azh = (ctx.get("dataset_version").find("AZH") == 0) ||
           (ctx.get("dataset_version").find("AToZHToLLTTbar") == 0);
  is_mc = ctx.get("dataset_type") == "MC";

  if(is_mc && !take_ntupleweights) m_pdfweights.reset(new PDFWeights(m_pdfname));

  pdf_idx = 9;
  if(is_dy || is_wjets || is_azh) pdf_idx = 45;

  sum_event_weights = book<TH1F>("sum_event_weights", "counting experiment", 1, 0.5, 1.5);

  h_murmuf_upup = book<TH1F>("h_murmuf_upup", "h_murmuf_upup", 1,0.5,1.5);
  h_murmuf_upnone = book<TH1F>("h_murmuf_upnone", "h_murmuf_upnone", 1,0.5,1.5);
  h_murmuf_noneup = book<TH1F>("h_murmuf_noneup", "h_murmuf_noneup", 1,0.5,1.5);
  h_murmuf_downdown = book<TH1F>("h_murmuf_downdown", "h_murmuf_downdown", 1,0.5,1.5);
  h_murmuf_downnone = book<TH1F>("h_murmuf_downnone", "h_murmuf_downnone", 1,0.5,1.5);
  h_murmuf_nonedown = book<TH1F>("h_murmuf_nonedown", "h_murmuf_nonedown", 1,0.5,1.5);

  if (is_dy || is_wjets || is_azh)
  {
    id_upup = 20;
    id_upnone = 5;
    id_noneup = 15;
    id_downdown = 40;
    id_downnone = 10;
    id_nonedown = 30;  
  } else {
    id_upup = 4;
    id_upnone = 1;
    id_noneup = 3;
    id_downdown = 8;
    id_downnone = 2;
    id_nonedown = 6;
  }

  for (int i = 0; i < 100; ++i)
  {
    h_pdf_hists.emplace_back(book<TH1F>(Form("h_pdf_%d", i), Form("h_pdf_%d", i), 1,0.5,1.5));
  }
}

void NormalisationHists::fill(const Event & event){
  double weight = event.weight;

  // Get muR, muF variations
  // event.genInfo->originalXWGTUP() = event.genInfo->systweights().at(0)
  double nominal = 1.0;

  double pdf_weight = 1.0;
  double murmuf_upup = 1.0;
  double murmuf_upnone = 1.0;
  double murmuf_noneup = 1.0;
  double murmuf_downdown = 1.0;
  double murmuf_downnone = 1.0;
  double murmuf_nonedown = 1.0;

  if(event.genInfo->systweights().size() >= 9){
      nominal = event.genInfo->originalXWGTUP();
      murmuf_upup = event.genInfo->systweights().at(id_upup);
      murmuf_upnone = event.genInfo->systweights().at(id_upnone);
      murmuf_noneup = event.genInfo->systweights().at(id_noneup);
      murmuf_downdown = event.genInfo->systweights().at(id_downdown);
      murmuf_downnone = event.genInfo->systweights().at(id_downnone);
      murmuf_nonedown = event.genInfo->systweights().at(id_nonedown);
  }

  // Now fill histograms with event.weight and scale variations
  sum_event_weights->Fill(1., weight);
  h_murmuf_upup->Fill(1, weight * murmuf_upup / nominal);
  h_murmuf_upnone->Fill(1, weight * murmuf_upnone / nominal);
  h_murmuf_noneup->Fill(1, weight * murmuf_noneup / nominal);
  h_murmuf_downdown->Fill(1, weight * murmuf_downdown / nominal);
  h_murmuf_downnone->Fill(1, weight * murmuf_downnone / nominal);
  h_murmuf_nonedown->Fill(1, weight * murmuf_nonedown / nominal);

  // if(event.genInfo->systweights().size() < 100 && take_ntupleweights) throw runtime_error("In NormalisationTools.cxx: Systweights in event.genInfo() is too small but ntupleweights shall be taken. Is this correct? In this case add exception to take_ntupleweights.");
  if(event.genInfo->systweights().size() > 110 && (!take_ntupleweights)) throw runtime_error("In NormalisationTools.cxx: Systweights in event.genInfo() is NOT empty but take_ntupleweights is set to 'false'. Is this correct? In this case Thomas says the genInfo weight should be used. Add this sample to take_ntupleweights");

  if(take_ntupleweights){
    for(int i=0; i<100; i++){
      if(event.genInfo->systweights().size()>=9){
        pdf_weight = event.genInfo->systweights().at(i+pdf_idx) / event.genInfo->originalXWGTUP();
      }
      h_pdf_hists[i]->Fill(1, weight * pdf_weight);
    }
  }
  else{
    std::vector<double> pdf_weights = m_pdfweights->GetWeightList(event);
    for(int i=0; i<100; i++){
      h_pdf_hists[i]->Fill(1, weight * pdf_weights[i]);
    }
  }
}

NormalisationHists::~NormalisationHists(){}


// PDF weights handler
PDFWeightHandleProducer::PDFWeightHandleProducer(Context & ctx){
  std::cout << "... Initialize PDF handler ..." << std::endl;

  is_mc = ctx.get("dataset_type") == "MC";
  is_dy = ctx.get("dataset_version").find("DYJets") == 0;
  is_wjets = ctx.get("dataset_version").find("WJets") == 0;
  is_azh = (ctx.get("dataset_version").find("AZH") == 0) ||
           (ctx.get("dataset_version").find("AToZHToLLTTbar") == 0);

  m_oname = ctx.get("dataset_version");
  TString m_pdfname = "NNPDF31_nnlo_hessian_pdfas";
  // TO CLARIFY:
  // Alternatively : PDF4LHC15_nnlo_100

  take_ntupleweights = !(m_oname.Contains("QCD"));

  pdf_idx = 9;
  if(is_dy || is_wjets || is_azh) pdf_idx = 45;

  if(is_mc && !take_ntupleweights) m_pdfweights.reset(new PDFWeights(m_pdfname));
  h_pdfweights = ctx.declare_event_output<vector<float>>("weight_pdf");

}

bool PDFWeightHandleProducer::process(Event & event){
  vector<float> pdf_weights;

  if((!is_mc) || (event.genInfo->systweights().size()<9)){
    for(int i=0; i<100; i++){
      pdf_weights.emplace_back(1.);
    }
    event.set(h_pdfweights, pdf_weights);
    return false;
  }

  // if(event.genInfo->systweights().size() < 100 && take_ntupleweights) throw runtime_error("In NormalisationTools.cxx: Systweights in event.genInfo() is too small but ntupleweights shall be taken. Is this correct? In this case add exception to take_ntupleweights.");
  if(event.genInfo->systweights().size() > 110 && (!take_ntupleweights)) throw runtime_error("In NormalisationTools.cxx: Systweights in event.genInfo() is NOT empty but take_ntupleweights is set to 'false'. Is this correct? In this case Thomas says the genInfo weight should be used. Add this sample to take_ntupleweights");

  double pdf_weight = 1.0;
  if(take_ntupleweights){
    for(int i=0; i<100; i++){
      if(event.genInfo->systweights().size() >= 9){
          pdf_weight = event.genInfo->systweights().at(i+pdf_idx) / event.genInfo->originalXWGTUP();
      }
      pdf_weights.emplace_back(pdf_weight);
    }
  }
  else{
    std::vector<double> weights = m_pdfweights->GetWeightList(event);
    for(int i=0; i<100; i++){
      pdf_weights.emplace_back(weights[i]);
    }
  }
  event.set(h_pdfweights, pdf_weights);
  return true;
}
