#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/NormFactsQCDPDF/include/NormalisationTools.h"


using namespace std;
using namespace uhh2;


// Prefix h_ denotes histogram
// Prefix s_ denotes selection
// Prefix handle_ denotes handle


class NormFacts: public AnalysisModule {
  public:

    explicit NormFacts(Context & ctx);
    virtual bool process(Event & event) override;

  private:
    std::unique_ptr<Hists> h_unc_norm;
    // Histograms for Cutflow
    bool is_mc;
};


NormFacts::NormFacts(Context & ctx){
  is_mc = ctx.get("dataset_type") == "MC";
  h_unc_norm.reset(new NormalisationHists(ctx, "UncNorms"));
  ctx.undeclare_all_event_output();
  }


bool NormFacts::process(Event & event) {
  // Fill histograms before applying any analysis cut
  // Fill histograms with nominal weight (event.weight)
  // and for all the variations (e.g. muR, muF) needed
  // to compute normalisation effects on these systematics
  if (is_mc) h_unc_norm->fill(event);
  return false;
}

UHH2_REGISTER_ANALYSIS_MODULE(NormFacts)

