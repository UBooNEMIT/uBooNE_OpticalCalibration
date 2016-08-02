Double_t MPEfit(Double_t * x, Double_t * par) {

// Defines the multi-photoelectron (MPE) fit used to fit charge and amplitude histograms.
// Fit is comprised of a series of poisson convolved gaussians.

  double ped_mu     = par[0];     // Centroid of the gaussian fit to the pedestal
  double ped_sig    = par[1];     // FWHM of the gaussian fit to the pedestal
  double mu_pe      = par[2];     // Average number of photoelectrons seen, determines the behavior of the poisson terms
  double spe        = par[3];     // Centroid of the 1-PE gaussian peak, used to extract the SPE value
  double pe_sig     = par[4];     // Intrincic FWHM of the photoelectron peaks, actual peaks will have FWHM including the baseline spread (ped_sig)
  double gain       = par[5];     // Peak to peak distance of photoelectron peaks, i.e. distance between the n PE and n+1 PE centroids
  double mpe_norm   = par[6];     // Scaling constant for the sum of poisson convolved gaussians

  double MPE = 0;

  int NPEtoFit = 1;

  for (int npe = 1; npe <= NPEtoFit ; npe++) {

    double subgauss_sigma = TMath::Sqrt(npe*pe_sig*pe_sig + ped_sig*ped_sig);
    MPE = MPE + TMath::Exp(-mu_pe)*TMath::Power(mu_pe,npe)/TMath::Factorial(npe) *TMath::Gaus(x[0],spe + ped_mu + (npe-1)*gain,subgauss_sigma,kTRUE);

  }

  MPE = MPE*mpe_norm;

  return MPE;

    }

//--------------------------------------------------------------------------------------------------------------------------------------------
// Actual fitting code follows                                                                                                             ---
//--------------------------------------------------------------------------------------------------------------------------------------------

std::vector<double> CalibrationFit(TString input_filename, int channel_num) {

  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TString input_tree_name = "specalib/pulsetree";                        // Whatever the tree from the swizzler is called
  string charge           = "charge";                                    // Whatever integrated peaks are called in the swizzled tree
  string amplitude        = "maxamp";                                    // Whatever peak maxima are called in the swizzled tree
  string baselinerms      = "baselinerms";                               // Whatever the baseline rms for the pre charge region is called in swizzled tree
  string baselinerms2     = "baselinerms2";                              // Whatever the baseline rms for the post charge region is called in swizzled tree
  string channel          = "opchannel";                                 // Whatever the channel number is called in the swizzled tree

  char base_ratio[100];
  sprintf(base_ratio,"%s / %s",baselinerms.c_str(),baselinerms2.c_str());


  int amp_min             = 10;                                          // Determines where the pedestal region ends for amplitude
  int amp_max             = 30;                                          // Determines where the region to be fit ends for amplitude
  int baseline_cutoff     = 3;                                           // Determines the maximum variability in baseline allowed
  double ratio_threshold  = 0.05;                                        // Determines the amount by which the post-charge baseline can vary from the pre-charge baseline
  

  TFile* in_file = new TFile(input_filename.Data());                     // Get input file output from swizzler
  TTree* in_tree = (TTree*)in_file->Get(input_tree_name);                // Get tree with PMT waveforms and data from input file


//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// This following segment histograms the charge and amplitudes and obtains an estimate for the position of the SPE peak to be                                  
// used to seed the actual fit later.                                                                                                                          

  TH1D* amp_histo    = new TH1D("amp_histo","Amplitude Histogram",100,0,100);

  char buffer1[300];
  sprintf(buffer1,"%s>%i && %s<%i && %s<%i && %s<(1+%f) && %s>(1-%f) && %s==%i",amplitude.c_str(),amp_min,baselinerms.c_str(),baseline_cutoff,baselinerms2.c_str(),baseline_cutoff,base_ratio,ratio_threshold,base_ratio,ratio_threshold,channel.c_str(),channel_num);

  in_tree->Project("amp_histo",amplitude.c_str(),buffer1,"goff",in_tree->GetEntries());

  double spe_amp_estimate = amp_histo->GetXaxis()->GetBinCenter(amp_histo->GetMaximumBin());

//----------------------------------------------------------------------------------------------------------------------------------------------------------------
// The following couple lines gives an estimate of the mpe normalization by simply finding the size of the 1 PE hump, note the                                 
// factor of 2.5066 is to account for the sqrt(2pi) norm factor in the gaussians, doesn't account for the sigma as we dont know                                
// that until we fit anyway, simply a rough initial starting point.                                                                                            

  double amp_mpe_norm_estimate = 2.5066*amp_histo->GetMaximum();
    

//---------------------------------------------------------------------------------------------------------------------------------------------------------
// Performs the actual fitting, seeding with the estimates made above                                                                                         

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(15000);  // Raises maximum minimization calls, decreases likelikhood of early termination

  TH1F* amp_fit_histo = new TH1F("amp_fit_histo","Amplitude Histogram With Fit",100,0,50);

  amp_fit_histo->SetTitle("Amplitude Distribution With Fit");
  amp_fit_histo->GetYaxis()->SetTitle("Counts");
  amp_fit_histo->GetXaxis()->SetTitle("Amplitude (ADCs)");

  char buffer8[300];
  sprintf(buffer8,"%s < %i && %s < %i && %s < %i && %s<(1+%f) && %s>(1-%f) && %s==%i",amplitude.c_str(),amp_max,baselinerms.c_str(),baseline_cutoff,baselinerms2.c_str(),baseline_cutoff,base_ratio,ratio_threshold,base_ratio,ratio_threshold,channel.c_str(),channel_num);

  in_tree->Project("amp_fit_histo",amplitude.c_str(),buffer8,"goff",in_tree->GetEntries());

  
  TF1* amp_fit  = new TF1("amp_fit",MPEfit,amp_min,amp_max,7);
  amp_fit->SetNpx(1000);
  Double_t amp_par[7] = {0,2,0.2,spe_amp_estimate,5,5,amp_mpe_norm_estimate};
  amp_fit->SetParNames("ped_mu","ped_sig","mu_pe","spe","pe_sig","gain","mpe_norm");
  amp_fit->SetParameters(amp_par);

  amp_fit->SetParLimits(0,-2,2);
  amp_fit->SetParLimits(1,0,4);
  amp_fit->SetParLimits(2,0,1);
  amp_fit->SetParLimits(3,0.8*spe_amp_estimate,1.2*spe_amp_estimate);
  amp_fit->SetParLimits(4,0,20);
  amp_fit->SetParLimits(5,0,50);
  amp_fit->SetParLimits(6,0,1e10);




  TCanvas *c1 = new TCanvas("c1","",10,10,800,600);
  c1->cd();
  c1->SetLogy();
  amp_fit_histo->Fit("amp_fit","R");

  
  std::vector<double> output;

  //  std::string status;
  //status = gMinuit->fCstatu;
 
  //  if ( status.compare("FAILED") == 0) {

  //  double singl = spe_amp_estimate;
  //  double sigma = 0;

  //}

  //else {

    double singl = amp_fit->GetParameter(3);
    double sigma = amp_fit->GetParError(3);

    //}


    //  output.push_back(singl);
  output.push_back(spe_amp_estimate);
  output.push_back(sigma);

  return output;
}
