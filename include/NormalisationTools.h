#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/common/include/PDFWeights.h"

class NormalisationHists: public uhh2::Hists {
public:
  explicit NormalisationHists(uhh2::Context&, const std::string&);
  virtual void fill(const uhh2::Event&) override;

protected:
  void init(uhh2::Context & ctx);
  bool is_mc, is_dy, is_wjets, is_azh, take_ntupleweights;

  TH1F *sum_event_weights;

  TH1F *h_murmuf_upup;
  TH1F *h_murmuf_upnone;
  TH1F *h_murmuf_noneup;
  TH1F *h_murmuf_downdown;
  TH1F *h_murmuf_downnone;
  TH1F *h_murmuf_nonedown;

  int pdf_idx;
  int id_upup;
  int id_upnone;
  int id_noneup;
  int id_downdown;
  int id_downnone;
  int id_nonedown;

  std::vector<TH1F *> h_pdf_hists;
  std::unique_ptr<PDFWeights> m_pdfweights;

  TString m_oname;

  virtual ~NormalisationHists();
};

// PDF weights handler
class PDFWeightHandleProducer: public uhh2::AnalysisModule{

public:
  explicit PDFWeightHandleProducer(uhh2::Context & ctx);
  virtual bool process(uhh2::Event & event) override;

private:
  uhh2::Event::Handle<std::vector<float>> h_pdfweights;
  bool is_mc, is_dy, is_wjets, is_azh, take_ntupleweights;
  int pdf_idx;
  TString m_oname;
  std::unique_ptr<PDFWeights> m_pdfweights;
};
