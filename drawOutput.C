// Define fitting functions.. Gausian and a quadratic background
Double_t fitGaus(Double_t *x,Double_t *par) {
  Double_t arg = 0;
  if(par[2]!=0){arg = (x[0] - par[1])/par[2];}
  Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg)/par[2]/pow(2*TMath::Pi(),0.5);
  return fitval;
}
Double_t fitBackground(Double_t *x, Double_t *par){
  Double_t arg = 0;
  Double_t fitval = par[0]+par[1]*x[0];
  return fitval;
}

Double_t fitf(Double_t *x, Double_t *par){
  Double_t fitval = fitBackground(x,par) + fitGaus(x,&par[2]);
  return fitval;
}

double calcScaleFactor(double a0,double a1,double a,double b,int num_background_events, double bin_width_Dalitz){
    //double scale_factor = 0.5*(b-a)(2*a0+a1*(a+b))/num_background_events;
    double scale_factor = a0*(b-a)+0.5*a1*(b*b-a*a);
    scale_factor = scale_factor/bin_width_Dalitz/num_background_events;
    return scale_factor;
}
double calcAsymError(double n_Pos,double n_Neg, double asymmetry){
  double error;
  double error = pow((1-asymmetry*asymmetry)/(n_Neg+n_Pos),0.5);
  return error;
  
}

/*
double dalitzFunctionMax(double m12, double m_B, double m_k){
  double E2 = 0.5*m12;
  double E3 = 0.5(m_B*m_B-m12*m12-m_k*m_k)/m12;
  double m23 = pow((E2+E3),2)-pow(
*/
void drawOutput() {
    //////////////////////
    // Example Root Macro for third year B->hhh Lab
    ////////////////////// 

    // Open the root file that was produced by running the example program
    TFile *f = new TFile("outputDataMagnetUp.root");
    
   //Define any new histograms needed 
    TH2F           *h_SignalBgdNeg;
    TH2F           *h_SignalBgdPos;
    TH2F           *h_Asymmetry;
    TH2F           *h_AsymmetryErr;
    TH2F           *h_WeightedAsymmetry;
    
  
    // Get pointers to the example histograms that were made 
    TH1F *hx = (TH1F*)f->Get("h_PX");
    TH1F *hy = (TH1F*)f->Get("h_PY");
    TH2F *hxy = (TH2F*)f->Get("h_TXTY");
    TH1F *HMB = (TH1F*)f->Get("h_MB");
    TH1F *hMBPos = (TH1F*)f->Get("h_MBPos");
    TH1F *hMBNeg = (TH1F*)f->Get("h_MBNeg");
    TH1F *hMResLower = (TH1F*)f->Get("h_MResLower");
    TH1F *hMResHigher = (TH1F*)f->Get("h_MResHigher");
    TH1F *hMResMuLower = (TH1F*)f->Get("h_MResMuLower");
    TH1F *hMResMuHigher = (TH1F*)f->Get("h_MResMuHigher");
    TH1F *hTotalCharge = (TH1F*)f->Get("h_TotalCharge");
    TH2F *hDalitzSignalBNeg = (TH2F*)f->Get("h_DalitzSignalBNeg");
    TH2F *hDalitzBgdBNeg = (TH2F*)f->Get("h_DalitzBgdBNeg");
    TH2F *hDalitzSignalBPos = (TH2F*)f->Get("h_DalitzSignalBPos");
    TH2F *hDalitzBgdBPos = (TH2F*)f->Get("h_DalitzBgdBPos");
    
    int num_bins_dalitz = hDalitzBgdBNeg->GetXaxis()->GetNbins();
    hAsymmetry = new  TH2F("h_Asymmetry","",num_bins_dalitz,0,35,num_bins_dalitz,0,35);
    hSignalBgdNeg = new TH2F("h_SignalBgdNeg","",num_bins_dalitz,0,35,num_bins_dalitz,0,35);
    hSignalBgdPos = new TH2F("h_SignalBgdPos","",num_bins_dalitz,0,35,num_bins_dalitz,0,35);
    hAsymmetryErr = new TH2F("h_AsymmetryErr","",num_bins_dalitz,0,35,num_bins_dalitz,0,35);
    hWeightedAsymmetry = new TH2F("h_WeightedAsymmetry","",num_bins_dalitz,0,35,num_bins_dalitz,0,35);
    hDalitzSignal = (TH2F*)f->Get("h_DalitzSignal");

/*
    // Create a canvas onto which the histograms are plotted and which can be saved
    TCanvas *c1 = new TCanvas("c1","",600,400);
    // Draw the first histogram with a blue line, and an x-axis title
    hx->SetLineColor(kBlue);
    hx->GetXaxis()->SetTitle("Momentum [MeV/c^{2}]");
    hx->Draw();

    // Draw the second histogram on the same plot with a red line
    hy->SetLineColor(kRed);
    //  hy->Draw("same");

    // Save the canvas as pdf file. Other possible formats include root (for later editing) eps, png.
    c1->SaveAs("hx_hy.pdf");

    // Repeatx_max_hDalitzBgdNeg the above for the 2D histogram
    TCanvas *c2 = new TCanvas("c2","",600,400);
    hxy->SetStats(0);                           // remove the statistics box
    hxy->GetXaxis()->SetTitle("Slope in x");    // add axis titles
    hxy->GetYaxis()->SetTitle("Slope in y");
    // hxy->Draw("colz");                          // draw with a colour scale
    c2->SaveAs("hxy.pdf");

    hx->SetLineColor(kBlue);drawOutput.C:69:
    //Mass of B plot, magnet up
    TCanvas *c3 = new TCanvas("c3","",800,400);
    hMB->SetTitle("B+/B- Invariant Mass histogram");
    hMB->GetXaxis()->SetTitle("Mass [MeV/c^{2}]");
    hMB->GetYaxis()->SetTitle("Number of Events");
    //hMB->Draw();
  num_bins_hDalitzBgNeg
    c3->SaveAs("hMB.pdf");
    
     double num_bins_hMBPos = hMBPos->GetNbinsX();
//   calculating the bin width of Dalitz plots
    double bin_min = hDaltizBgdNeg.GetXaxis()->GetXmin();
    double bin_max = hDalitzBgdNeg.GetXaxis()->GetXMax();
    
    double bin_width_Dalitz = (bin_max-bin_min)/num_bins_dalitz;
    */
    TCanvas *c4 = new TCanvas("c4","",800,400);
    hMBPos->SetTitle("B+ Invariant Mass histogram");
    hMBPos->GetXaxis()->SetTitle("Mass [MeV/c^{2}]");
    hMBPos->GetYaxis()->SetTitle("Number of Events");
   // hMBPos->Draw();
   double bin_width_B = hMBPos->GetXaxis()->GetBinWidth(1);
    TF1 *funcPos = new TF1("fitPos",fitf,5.2e3,5.4e3,5);// set the parameters to the mean and RMS of the histogram
    funcPos->SetParameters(0,1e-3,2e2,hMBPos->GetMean(),hMBPos->GetRMS());// give the parameters meaningful names
     funcPos->SetParNames ("a0","a1","Constant","Mean_value","Sigma");// call TH1::Fit with the name of the TF1 object
    hMBPos->Fit("fitPos");
    
    c4->SaveAs("B+ mass histogram.pdf");

    TCanvas *c5 = new TCanvas("c5","",800,400);
    hMBNeg->SetTitle("B- Invariant Mass histogram");
    hMBNeg->GetXaxis()->SetTitle("Mass [MeV/c^{2}]");
    hMBNeg->GetYaxis()->SetTitle("Number of Events");
    //hMBNeg->Draw();
    TF1 *funcNeg = new TF1("fitNeg",fitf,5.2e3,5.4e3,5);// set the parameters to the mean and RMS of the histogram
     funcNeg->SetParameters(0,1e-3,2e2,hMBNeg->GetMean(),hMBNeg->GetRMS());// give the parameters meaningful names
     funcNeg->SetParNames ("a0","a1","Constant","Mean_value","Sigma");// call TH1::Fit with the name of the TF1 object
     hMBNeg->Fit("fitNeg");
     /*
     TF1 *backgroundNeg = new TF1("fitBackground",fitBackground,5.2e3,5.4e3,2);
     TF1 *fitNeg = hMBNeg->GetFunction("fitNeg");
     double a0 = funcNeg->GetParameter(0);
     double a1 = funcNeg->GetParameter(1);
     fitNeg->SetParameters(a0,a1);
     fitNeg->Draw();
     */
    c5->SaveAs("B- mass histogram.pdf");
    
    //Resonaces plot
    TCanvas *c6 = new TCanvas("c6","", 800, 400);
    hMResLower->SetTitle("Resonances Histogram");
    hMResLower->GetXaxis()->SetTitle("[MeV/c^{2}]");
    hMResLower->GetYaxis()->SetTitle("Number of Events");
    hMResLower->SetLineColor(kRed);
    
    
    c6->SaveAs("Resonances histogram.pdf");
    
//     TCanvas *c7 = new TCanvas("c7", "", 800, 400);
    //hMResHigher->SetTitle("Resonances Histogram Higher mass");
    //hMResHigher->GetXaxis()->SetTitle("[MeV/c^{2}]");
   // hMResHigher->GetYaxis()->SetTitle("Number of Events");
    hMResHigher->SetLineColor(kBlue);
    hMResLower->Draw();
    hMResHigher->Draw("SAME");
    
    TCanvas *c19 = new TCanvas("c19", "", 800, 400);
    
    hMResMuLower->SetLineColor(kGreen);
    hMResMuLower->Draw();


    hMResMuHigher->SetLineColor(kYellow);
    hMResMuHigher->Draw("SAME");
    /*
    TCanvas *c8 = new TCanvas("c8", "", 800, 400);
    hTotalCharge->SetTitle("Total Charge");
Bool_t MyAnalysis::CutResonances(double MRes1, double MRes2){
  //Returns true if either of the resonances are in the viscinity of MDMeson
  //Used to deselect decay channels which proceeded via the D meson
  bool valid = false;
  double MDMeson = 1864.84;//Mass of D meson
  double  MRange = 30;//The range in which we will cut the decay chanel
  if(MRes1>MDMeson-
    hTotalCharge->GetXaxis()->SetTitle("[MeV/c^{2}]");
    hTotalCharge->GetYaxis()->SetTitle("Number of Events");
    hTotalCharge->SetLineColor(kBlue);
    hTotalCharge->Draw();
    */
    
    //############    B- Dalitz plots     ##################################################
    
    
   
    
  
  
  
  
  
    TCanvas *c9 = new TCanvas("c9", "", 1200, 1200);
    hDalitzBgdBNeg->SetTitle("Background Dalitz Plot B-");
    hDalitzBgdBNeg->GetXaxis()->SetTitle("mlower^{2}[GeV^{2}/c^{4}]");
    hDalitzBgdBNeg->GetYaxis()->SetTitle("mhigher^{2}[GeV^{2}/c^{4}]");
    hDalitzBgdBNeg->Draw("colz");
    c9->SaveAs("hDalitzBackgroundB-.pdf");
    
    double bin_width_Dalitz = hDalitzBgdBNeg->GetXaxis()->GetBinWidth(1);
    double num_background_events_neg = hDalitzBgdBNeg->GetEntries();
    TF1 *fitNeg = hMBNeg->GetFunction("fitNeg");
    double a0 = fitNeg->GetParameter(0);
    double a1 = fitNeg->GetParameter(1);
    double constant = fitNeg->GetParameter(2);
    double a = 5200;
    double b = 5400;      
    double scale_factor_B_Neg = calcScaleFactor(a0,a1,a,b,num_background_events_neg, bin_width_B);
    
    double num_background = a0*(b-a)+0.5*a1*(b*b-a*a);
    num_background = num_background/bin_width_B;
    double num_signal = constant/bin_width_B;
    double maximization = num_signal/pow(num_signal+num_background,0.5);
    std::cout<<"\n\n\n maximization parameter: "<<maximization<<"\n\n\n";
    std::cout<<"Neg scale factor: "<<scale_factor_B_Neg<<"\n";
    TCanvas *c10 = new TCanvas("c10", "", 1200, 1200);
    hDalitzSignalBNeg->SetTitle("Signal Dalitz Plot B-");
    
    hDalitzSignalBNeg->GetXaxis()->SetTitle("mlower^{2}[GeV^{2}/c^{4}]");
    hDalitzSignalBNeg->GetYaxis()->SetTitle("mhigher^{2}[GeV^{2}/c^{4}]");
    
    
    TH2F *hSignalBgdNeg = (TH2F*)hDalitzSignalBNeg->Clone();//Make a copy of the signal plot to do the subtraction

    hSignalBgdNeg->Add(hDalitzBgdBNeg,-scale_factor_B_Neg);//Perform background subtraction
    int num_bins_dalitz = hSignalBgdNeg->GetXaxis()->GetNbins();
for(int i=1;i<num_bins_dalitz+1;i++){
      for(int j=1;j<num_bins_dalitz+1;j++){
	double content = hSignalBgdNeg->GetBinContent(i,j);
	if(content<0){
	  hSignalBgdNeg->SetBinContent(i,j,0);}
      }
}
    
    hDalitzSignalBNeg->Draw("colz");
    c10->SaveAs("hDalitzSignalBNeg.pdf");

    //#################    B+ Dalitz plots    ###################################
    TCanvas *c11 = new TCanvas("c11", "", 1200, 1200);
    hDalitzBgdBPos->SetTitle("Background Dalitz Plot B+");
    hDalitzBgdBPos->GetXaxis()->SetTitle("mlower^{2}[GeV^{2}/c^{4}]");
    hDalitzBgdBPos->GetYaxis()->SetTitle("mhigher^{2}[GeV^{2}/c^{4}]");
    hDalitzBgdBPos->Draw("colz");
//     c11->SaveAs("hDalitzBackgroundB+.pdf");  
    
    double num_background_events_pos = hDalitzBgdBPos->GetEntries();
    std::cout<<num_background_events_pos;
    TF1 *fitPos = hMBPos->GetFunction("fitPos");
    double a0 = fitPos->GetParameter(0);
    double a1 = fitPos->GetParameter(1);
         
    double scale_factor_B_Pos = calcScaleFactor(a0,a1,a,b,num_background_events_pos, bin_width_B);
    std::cout<<"Pos Scale factor: "<< scale_factor_B_Pos<<"\n";
    TCanvas *c12 = new TCanvas("c12", "", 1200, 1200);
    hDalitzSignalBPos->SetTitle("Signal Dalitz Plot B+");
    hDalitzSignalBPos->GetXaxis()->SetTitle("mlower^{2}[GeV^{2}/c^{4}]");
    hDalitzSignalBPos->GetYaxis()->SetTitle("mhigher^{2}[GeV^{2}/c^{4}]");
    
    TH2F *hSignalBgdPos = (TH2F*)hDalitzSignalBPos->Clone();//Make a copy of the signal plot to do the subtraction

    hSignalBgdPos->Add(hDalitzBgdBPos,-scale_factor_B_Pos);//Perform background subtraction
    for(int i=1;i<num_bins_dalitz+1;i++){
      for(int j=1;j<num_bins_dalitz+1;j++){
	double content = hSignalBgdPos->GetBinContent(i,j);
	if(content<0){
	  hSignalBgdPos->SetBinContent(i,j,0);}
      }
}
    hDalitzSignalBPos->Draw("colz");
    c12->SaveAs("hDalitzSignalBPos.pdf");
        
    //##################    Aslusymetry plots    #############################################
    /*
    TCanvas *c14 = new TCanvas("c14", "", 1200, 1200);
    hSignalBgdNeg->SetTitle("Asymmetry plot");
    hSignalBgdNeg->GetXaxis()->SetTitle("mlower^{2}[GeV^{2}/c^{4}]");
    hSignalBgdNeg->GetYaxis()->SetTitle("mhigher^{2}[GeV^{2}/c^{4}]");
    hSignalBgdNeg = hDalitzSignalBNeg->Add(hDalitzBgdBNeg,-scale_factor_B_Neg);    
    hSignalBgdNeg->Draw("colz");
    c14->SaveAs("hAsymmetry.pdf");
    
    */

    
    
    
   TCanvas *c14 = new TCanvas("c14", "", 1200, 1200);
    hSignalBgdNeg->SetTitle("Signal Dalitz Plot with background subtraction B-");
    hSignalBgdNeg->GetXaxis()->SetTitle("mlower^{2}[GeV^{2}/c^{4}]");
    hSignalBgdNeg->GetYaxis()->SetTitle("mhigher^{2}[GeV^{2}/c^{4}]");
    hSignalBgdNeg->Draw("colz");
    c14->SaveAs("hDalitzSignalBPos.pdf");
    
    TCanvas *c15 = new TCanvas("c15", "", 1200, 1200);
    hSignalBgdPos->SetTitle("Signal Dalitz Plot with background subtraction B+");
    hSignalBgdPos->GetXaxis()->SetTitle("mlower^{2}[GeV^{2}/c^{4}]");
    hSignalBgdPos->GetYaxis()->SetTitle("mhigher^{2}[GeV^{2}/c^{4}]");
    hSignalBgdPos->Draw();
    c15->SaveAs("hDalitzSignalBPos.pdf");
    
    //Asymmetry plot
    TCanvas *c13 = new TCanvas("c13", "", 1200, 1200);
        
   //TH2F *hAsymmetry = (TH2F*)hSignalBgdPos->Clone();

   
    hAsymmetry=(TH2F*) hSignalBgdPos->GetAsymmetry(hSignalBgdNeg);
     hAsymmetry->SetTitle("Asymmetry plot");
    hAsymmetry->GetXaxis()->SetTitle("mlower^{2}[GeV^{2}/c^{4}]");
    hAsymmetry->GetYaxis()->SetTitle("mhigher^{2}[GeV^{2}/c^{4}]");
    hAsymmetry->Draw("colz");
    c13->SaveAs("hAsymmetry.pdf");
    
    
    //Calculate the uncertainty on each bin
    for(int i=1;i<num_bins_dalitz+1;i++){
      for(int j=1;j<num_bins_dalitz+1;j++){
	double n_Pos = hSignalBgdPos->GetBinContent(i,j);
	double n_Neg = hSignalBgdNeg->GetBinContent(i,j);
	double asymmetry = hAsymmetry->GetBinContent(i,j);
	double error = calcAsymError(n_Pos,n_Neg,asymmetry);
	double n_Pos_Sig = hDalitzSignalBPos->GetBinContent(i,j);
	double n_Neg_Sig = hDalitzSignalBNeg->GetBinContent(i,j);
	if(!(n_Pos_Sig+n_Neg_Sig ==0)){
	hAsymmetryErr->Fill(35/num_bins_dalitz*i-0.1,35/num_bins_dalitz*j-0.1,error);
	hWeightedAsymmetry->Fill(35/num_bins_dalitz*i-0.1,35/num_bins_dalitz*j-0.1,asymmetry/error);}
      }
    }
    
     
     TCanvas *c16 = new TCanvas("c16", "", 1200, 1200);
    hAsymmetryErr->SetTitle("Asymmetry error plot");
    hAsymmetryErr->GetXaxis()->SetTitle("mlower^{2}[GeV^{2}/c^{4}]");
    hAsymmetryErr->GetYaxis()->SetTitle("mhigher^{2}[GeV^{2}/c^{4}]");
    hAsymmetryErr->Draw("colz");
    c16->SaveAs("hDalitzSignalBPos.pdf");
    
    TCanvas *c17 = new TCanvas("c17", "", 1200, 1200);
    hWeightedAsymmetry->SetTitle("Weighted Asymmetry");
    hWeightedAsymmetry->GetXaxis()->SetTitle("mlower^{2}[GeV^{2}/c^{4}]");
    hWeightedAsymmetry->GetYaxis()->SetTitle("mhigher^{2}[GeV^{2}/c^{4}]");
    hWeightedAsymmetry->Draw("colz");
    c17->SaveAs("hDalitzSignalBPos.pdf");
    
    TCanvas *c18 = new TCanvas("c18", "", 1200, 1200);
    hDalitzSignal->SetTitle("Dalitz plot");
    hDalitzSignal->GetXaxis()->SetTitle("mlower^{2}[GeV^{2}/c^{4}]");
    hDalitzSignal->GetYaxis()->SetTitle("mhigher^{2}[GeV^{2}/c^{4}]");
    hDalitzSignal->Draw("colz");
    c18->SaveAs("hDalitzSignal.pdf");

    

   
   
    
     
    
     
     
    
    


    
     


}
