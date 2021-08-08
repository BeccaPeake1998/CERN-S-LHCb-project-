#include "Analysis.hpp"

// This is the analysis class, which realises the generic Analysis
// from Analysis.hpp
//
// Look in Analysis.hpp for the event variables available.
class MyAnalysis : public Analysis {
public:
    // Define your histograms here
    TH1F           *h_PX;
    TH1F           *h_PY;
    TH1F           *h_PZ;
    TH2F           *h_TXTY;
    TH1F           *h_ProbK;
    TH1F           *h_MB;
    TH1F           *h_MBNeg;
    TH1F           *h_MBPos;
    TH1F           *h_MResMuLower;
    TH1F           *h_MResMuHigher;
    TH1F           *h_MResLower;//Histogram for lower-mass resonance pairs
    TH1F           *h_MResHigher;//Histogram for higher-mass resonance pairs
    TH1F           *h_TotalCharge;
    TH2F           *h_DalitzSignal;
    TH2F           *h_DalitzSignalBNeg;
    TH2F           *h_DalitzBgdBNeg;
    TH2F           *h_DalitzSignalBPos;
    TH2F           *h_DalitzBgdBPos;
    TH2F           *h_Asymmetry;
    TH2F           *h_SignalBgdNeg;
  

    void     BookHistos();

    Bool_t   Cut();
    void     Execute();
    Bool_t   CutResonances(double,double);
};

void MyAnalysis::BookHistos()
{
    // This function is only called once at the start of the program.
    // Book your histograms here. The format is object_name,
    // histogram_name, number_of_bins, minimum, maximum For a 2D
    // histogram, use TH2F with first the number of bins and limits
    // for the x axis and then for the y axis
    //
    // push_back() adds the histograms to a vector v_Histos.  This
    // will take care of writing out histograms in
    // Analysis::SaveHistos
    int num_bins = 200;
    int num_bins_dalitz = 14;
    v_Histos.push_back( h_PX   = new TH1F("h_PX",  "", 100, -1e4, 1e4) );
    v_Histos.push_back( h_PY   = new TH1F("h_PY",  "", 100, -1e4, 1e4) );
    v_Histos.push_back( h_PZ   = new TH1F("h_PZ",  "", 100, -1e5, 1e5) );
    v_Histos.push_back( h_TXTY = new TH2F("h_TXTY","", 100, -1,1, 100,-1, 1) );
    v_Histos.push_back( h_ProbK = new TH1F("h_ProbK","", 100, 0,1));
    v_Histos.push_back( h_MB = new TH1F("hMB","", 300, 5.2e3, 6.3e3));
    v_Histos.push_back( h_MBNeg = new TH1F("h_MBNeg","",num_bins, 5.2e3,6.3e3));
    v_Histos.push_back( h_MBPos = new TH1F("h_MBPos", "", num_bins, 5.2e3, 6.3e3));
    v_Histos.push_back( h_MResMuLower = new TH1F("h_MResMuLower","",300,987.4,4.79e3));
    v_Histos.push_back( h_MResMuHigher = new TH1F("h_MResMuHigher","",300,987.4,4.79e3));
    v_Histos.push_back( h_MResLower = new TH1F("h_MResLower","",300,987.4,4.79e3));
    v_Histos.push_back( h_MResHigher = new TH1F("h_MResHigher","",300,987.4,4.79e3));    
    v_Histos.push_back( h_TotalCharge = new TH1F("h_TotalCharge","",200,-10,10));
    v_Histos.push_back( h_DalitzSignalBNeg = new TH2F("h_DalitzSignalBNeg","",num_bins_dalitz, 0,35,num_bins_dalitz,0,35));
    v_Histos.push_back( h_DalitzBgdBNeg = new TH2F("h_DalitzBgdBNeg","",num_bins_dalitz, 0,35,num_bins_dalitz,0,35));
    v_Histos.push_back( h_DalitzSignalBPos = new TH2F("h_DalitzSignalBPos","",num_bins_dalitz, 0,35,num_bins_dalitz,0,35));
    v_Histos.push_back( h_DalitzBgdBPos = new TH2F("h_DalitzBgdBPos","",num_bins_dalitz, 0,35,num_bins_dalitz,0,35));
    v_Histos.push_back( h_DalitzSignal = new TH2F("h_DalitzSignal","",num_bins_dalitz, 0,35,num_bins_dalitz,0,35));

}

Bool_t MyAnalysis::Cut(){
    // This function is called for every event from the Execute
    // function to define whether or not to accept this event.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    // This example checks if the PZ component of particle 3 is greater than 0. 
  bool valid=true;
  double probKFloor=0.3;
  if (H1_ProbK<probKFloor||H2_ProbK<probKFloor||H2_ProbK<probKFloor){valid = false;}
  if (H1_ProbK*H2_ProbK*H3_ProbK<pow(probKFloor,3)){valid = false;}
  
  if(H1_isMuon == true||H2_isMuon==true||H3_isMuon==true){
    valid = false;
  }
  return valid;
}

Bool_t MyAnalysis::CutResonances(double MRes1, double MRes2){
  //Returns true if either of the resonances are in the viscinity of MDMeson
  //Used to deselect decay channels which proceeded via the D meson
  bool cut = false;
  double MDMeson = 1864.84;//Mass of D meson
  double  MDRange = 30;//The range in which we will cut the decay chanel
  double MJMeson = 3096.92;
  double MJRange = 50;
//   if(MRes1>MDMeson-MDRange && MRes1<MDMeson+MDRange){cut = true;}
//   if(MRes2>MDMeson-MDRange && MRes2<MDMeson+MDRange){cut = true;}
//   if(MRes1>MJMeson-MJRange&&MRes1<MJMeson+MJRange){cut = true;}
//   if(MRes2>MJMeson-MJRange&&MRes2<MJMeson+MJRange){cut = true;}
  return cut;
}
void MyAnalysis::Execute(){//define some rest masses of the kaons
  //  Double_t H1_P4 [4]; //Particle 1 4-momentum  
  // Double_t H2_P4 [4];
  //  Double_t H3_P4 [4];
  //  Double_t HB_P4 [4];//4 momentum B meson
  Double_t MBMeson;//Mass B meson
  Double_t H1_PT ; //Particle 1 total momentum
  Double_t H2_PT;
  Double_t H3_PT;
  const Double_t RestMKaon=493.7;
  const Double_t MMuon = 105.66;

    // This method gets called on every event.
    // In this example the momentum components are filled into histograms.

    // Call the Cut function to decide whether to plot this event or not
    // it returns if the cut function returns false
    if ( !Cut() )
	return;
    // Fill your histograms below.
    // fill the momentum of all three particles 
    
    h_PX->Fill( H1_PX );
    h_PX->Fill( H2_PX );
    h_PX->Fill( H3_PX );
    // the PY of all three particles
    h_PY->Fill( H1_PY );
    h_PY->Fill( H2_PY );
    h_PY->Fill( H3_PY );
    // the PZ of all three particles
    h_PZ->Fill( H1_PZ );
    h_PZ->Fill( H2_PZ );
    h_PZ->Fill( H3_PZ );
    // 2D histogram of PX/PZ vs PY/PZ
    h_TXTY->Fill( H1_PX / H1_PZ, H1_PY / H1_PZ );
    //1D histogram for prob particle is Kaon
    h_ProbK->Fill( H1_ProbK);
    h_ProbK->Fill( H2_ProbK);
    h_ProbK->Fill( H3_ProbK);
    //Total momentums of the 3 decay products
    H1_PT =pow(H1_PX*H1_PX + H1_PY*H1_PY+H1_PZ*H1_PZ,0.5);
    H2_PT =pow(H2_PX*H2_PX + H2_PY*H2_PY + H2_PZ*H2_PZ,0.5);
    H3_PT =pow( H3_PX*H3_PX + H3_PY*H3_PY + H3_PZ*H3_PZ,0.5);
    //The 4 momentums of the decay products and then the B particle
    Double_t H1_P4[4] = {pow(pow(H1_PT,2)+RestMKaon*RestMKaon,0.5), H1_PX, H1_PY, H1_PZ};
    Double_t H2_P4[4] = {pow(pow(H2_PT,2)+RestMKaon*RestMKaon,0.5), H2_PX, H2_PY, H2_PZ};
    Double_t H3_P4[4] = {pow(pow(H3_PT,2)+RestMKaon*RestMKaon,0.5), H3_PX, H3_PY, H3_PZ};
    //Muon misidentification resonances
    Double_t H1_P4Mu[4] = {pow(pow(H1_PT,2)+MMuon*MMuon,0.5), H1_PX, H1_PY, H1_PZ};
    Double_t H2_P4Mu[4] = {pow(pow(H2_PT,2)+MMuon*MMuon,0.5), H2_PX, H2_PY, H2_PZ};
    Double_t H3_P4Mu[4] = {pow(pow(H3_PT,2)+MMuon*MMuon,0.5), H3_PX, H3_PY, H3_PZ};
    
    Double_t B_P4[4];
    for (int i=0;i<4;i++){
      B_P4[i] = H1_P4[i]+H2_P4[i]+H3_P4[i];
    }
     
    MBMeson=pow(pow(B_P4[0],2)-pow(B_P4[1],2)-pow(B_P4[2],2)-pow(B_P4[3],2),0.5);
    
    //Calc masses of resonances. This is calculated for every event
    Double_t MRes ;
    Double_t MRes1 = 0;
    Double_t MRes2 = 0;
    Double_t Res_P4[4]; 
    
    Double_t MResMu ;
    Double_t MRes1Mu = 0;
    Double_t MRes2Mu = 0;
    Double_t Res_P4Mu[4];
    if(H1_Charge!=H2_Charge){
      for (int i=0;i<4;i++){
	Res_P4[i] = H1_P4[i]+H2_P4[i];
        Res_P4Mu[i] = H1_P4Mu[i]+H2_P4Mu[i];}
      MRes=pow(pow(Res_P4[0],2)-pow(Res_P4[1],2)-pow(Res_P4[2],2)-pow(Res_P4[3],2),0.5);
      MResMu=pow(pow(Res_P4Mu[0],2)-pow(Res_P4Mu[1],2)-pow(Res_P4Mu[2],2)-pow(Res_P4Mu[3],2),0.5);
      MRes1 = MRes;
      MRes1Mu = MResMu;
    }
      if(H2_Charge!=H3_Charge){
      for (int i=0;i<4;i++){
	Res_P4[i] = H2_P4[i]+H3_P4[i];
	Res_P4Mu[i] = H2_P4Mu[i]+H3_P4Mu[i];}
      MRes=pow(pow(Res_P4[0],2)-pow(Res_P4[1],2)-pow(Res_P4[2],2)-pow(Res_P4[3],2),0.5);
      MResMu=pow(pow(Res_P4Mu[0],2)-pow(Res_P4Mu[1],2)-pow(Res_P4Mu[2],2)-pow(Res_P4Mu[3],2),0.5);

      if(MRes1==0){
	MRes1 = MRes;
	MRes1Mu = MResMu;
      }else{
	MRes2=MRes;
	MRes2Mu=MResMu;
      }
    }
      if(H1_Charge!=H3_Charge){
      for (int i=0;i<4;i++){
	Res_P4[i] = H1_P4[i]+H3_P4[i];
      	Res_P4Mu[i] = H1_P4Mu[i]+H3_P4Mu[i];}

      MRes=pow(pow(Res_P4[0],2)-pow(Res_P4[1],2)-pow(Res_P4[2],2)-pow(Res_P4[3],2),0.5);
      MResMu=pow(pow(Res_P4Mu[0],2)-pow(Res_P4Mu[1],2)-pow(Res_P4Mu[2],2)-pow(Res_P4Mu[3],2),0.5);

      MRes2 = MRes;
      MRes2Mu=MResMu;
      }
    if(CutResonances(MRes1, MRes2)){
      return;}

    //####   Fill histograms   ####################################################################
   h_MB->Fill(MBMeson);//B meson histograms
   int totalCharge = H1_Charge+H2_Charge+H3_Charge;

   if(!CutResonances(MRes1, MRes2)){//Cut out D meson decay chanels
    h_TotalCharge->Fill(totalCharge);
    if(totalCharge>=0){      
      h_MBPos->Fill(MBMeson);//B+ meson histogram with D decays cut
    }else{
      h_MBNeg->Fill(MBMeson);}//B- histogram
   
    
   if(MBMeson>5180 && MBMeson<5350){//Now we cut to events around the B meson only
      if(MRes1<=MRes2){
      h_MResLower->Fill(MRes1);//Higher/lower mass resonance histograms
      h_MResHigher->Fill(MRes2);
      h_MResMuLower->Fill(MRes1Mu);
      h_MResMuHigher->Fill(MRes2Mu);
   }else{
      h_MResLower->Fill(MRes2);
      h_MResHigher->Fill(MRes1);
      h_MResMuLower->Fill(MRes2Mu);
      h_MResMuHigher->Fill(MRes1Mu);
   }
   }
   }
   
    //Fill Dalitz plot cutting around the signal event.  
   if(MBMeson>5200 && MBMeson<5400){
     
     if(totalCharge>=0){
      if(MRes1<=MRes2){h_DalitzSignalBPos->Fill(MRes1*MRes1/1e6, MRes2*MRes2/1e6);
      h_DalitzSignal->Fill(MRes1*MRes1/1e6, MRes2*MRes2/1e6);}
      else{h_DalitzSignalBPos->Fill(MRes2*MRes2/1e6, MRes1*MRes1/1e6);
      h_DalitzSignal->Fill(MRes2*MRes2/1e6, MRes1*MRes1/1e6);}
     }
     else{
      if(MRes1<=MRes2){h_DalitzSignalBNeg->Fill(MRes1*MRes1/1e6, MRes2*MRes2/1e6);
      h_DalitzSignal->Fill(MRes1*MRes1/1e6, MRes2*MRes2/1e6);}
      else{h_DalitzSignalBNeg->Fill(MRes2*MRes2/1e6, MRes1*MRes1/1e6);
      h_DalitzSignal->Fill(MRes2*MRes2/1e6, MRes1*MRes1/1e6);}
     }
   }
      
   
   //fill background Dalitz plot choosing only the background events
   if(MBMeson>5600&&MBMeson<6300){
     if(totalCharge>=0){
      if(MRes1<=MRes2){h_DalitzBgdBPos->Fill(MRes1*MRes1/1e6, MRes2*MRes2/1e6);}
      else{h_DalitzBgdBPos->Fill(MRes2*MRes2/1e6, MRes1*MRes1/1e6);}
     }
      else{
	if(MRes1<=MRes2){h_DalitzBgdBNeg->Fill(MRes1*MRes1/1e6, MRes2*MRes2/1e6);
      }else{h_DalitzBgdBNeg->Fill(MRes2*MRes2/1e6, MRes1*MRes1/1e6);}	
     }
}}

    
      


// The main function just calls the generic AnalysisMain function
// with the MyAnalysis class
//
// Normally you don't need to change this
int main(int argc, char* argv[])
{
    MyAnalysis* ana = new MyAnalysis();
    int res = ana->AnalysisMain(argc, argv);
    return res;
}
