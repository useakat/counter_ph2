#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <complex>
#include "TROOT.h"
#include "TApplication.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TMath.h"
#include "TLegend.h"
#include "TPaveText.h"

const double C = 29.9792458;

double GetArg( double r, double i )
{
  std::complex< double > c( r, i );
  return std::arg( c );
}


int Display( char* filename )
{
  gROOT->ProcessLine( "./ ~/.rootlogon.C" );

  gStyle->SetLabelSize(0.07, "xyz");
  gStyle->SetTitleSize(0.07, "xyz");
  gStyle->SetTitleOffset( 0.6, "y" );
  gStyle->SetTitleFontSize(0.07);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadTopMargin(0.12);
  gStyle->SetPadBottomMargin(0.15);


  ifstream ifs( filename );

  // Header info preparation
  double TargetLength = 0;
  double DeltaZ = 0;
  int Xi = 0;
  int EdgeXi = 0;
  double SimulationTime = 0;
  double SimTimePitch = 0;
  double OutputTimePitch = 0;
  int OutputXiPitch = 0;
  double Pressure = 0;
  double Temperature = 0;
  double NumberDensity = 0;
  double Omega_0 = 0;
  double Detuning = 0;
  double Gamma_1g = 0;
  double Gamma_1e = 0;
  double Gamma_2 = 0;
  double ERpIntensity = 0;
  double ELmIntensity = 0;
  double ERpWidth = 0;
  double ELmWidth = 0;
  double ERpDelay = 0;
  double ELmDelay = 0;
  double TrigIntensity = 0;
  double TrigWidth = 0;
  double TrigDelay = 0;
  bool Resonance = false;
  bool Absorption = false;
  bool Relaxation = false;
  bool Square = false;
  int EquationType = 0;
  int SpatialAlgorithm = 0;
  int TimeAlgorithm = 0;
  double AbsError = 0;
  double RelError = 0;
  int openmp = 0;


  std::string line;
  std::getline( ifs, line );
  while( line[0] == '#' ){
    std::istringstream iss(line );
    std::string header;
    std::string dummyString;
    iss >> header >> header;
    if( header == "TargetLength" ){ iss >> TargetLength; }
    else if( header == "DeltaZ"  ){ iss >> DeltaZ; }
    else if( header == "Xi"      ){ iss >> Xi; }
    else if( header == "EdgeXi"      ){ iss >> EdgeXi; }
    else if( header == "SimulationTime" ){ iss >> SimulationTime; }
    else if( header == "SimTimePitch" ){ iss >> SimTimePitch; }
    else if( header == "OutputTimePitch" ){ iss >> OutputTimePitch; }
    else if( header == "OutputXiPitch" ){ iss >> OutputXiPitch; }
    else if( header == "Pressure" ){ iss >> Pressure; }
    else if( header == "Temperature" ){ iss >> Temperature; }
    else if( header == "NumberDensity" ){ iss >> NumberDensity; }
    else if( header == "omega0" ){ iss >> Omega_0; }
    else if( header == "Detuning" ){ iss >> Detuning; }
    else if( header == "Gamma1g" ){ iss >> Gamma_1g; }
    else if( header == "Gamma1e" ){ iss >> Gamma_1e; }
    else if( header == "Gamma2" ){ iss >> Gamma_2; }
    else if( header == "ERpIntensity" ){ iss >> ERpIntensity; }
    else if( header == "ELmIntensity" ){ iss >> ELmIntensity; }
    else if( header == "ERpWidth" ){ iss >> dummyString >> ERpWidth; }
    else if( header == "ELmWidth" ){ iss >> dummyString >> ELmWidth; }
    else if( header == "ERpDelay" ){ iss >> ERpDelay; }
    else if( header == "ELmDelay" ){ iss >> ELmDelay; }
    else if( header == "TrigIntensity" ){ iss >> TrigIntensity; }
    else if( header == "TrigWidth" ){ iss >> dummyString >> TrigWidth; }
    else if( header == "TrigDelay" ){ iss >> TrigDelay; }
    else if( header == "ForceResonance" ){ iss >> Resonance; }
    else if( header == "Absorption" ){ iss >> Absorption; }
    else if( header == "Relaxation" ){ iss >> Relaxation; }
    else if( header == "Square" ){ iss >> Square; }
    else if( header == "Equation" ){ iss >> EquationType; }
    else if( header == "SpatialAlgorithm" ){ iss >> SpatialAlgorithm; }
    else if( header == "TimeAlgorithm" ){ iss >> TimeAlgorithm; }
    else if( header == "AbsoluteError" ){ iss >> AbsError; }
    else if( header == "RelativeError" ){ iss >> RelError; }
    else if( header == "OpenMP" ){ iss >> openmp; }
    else if( header == "Parameter" ){ break; }

    std::getline( ifs, line );
  }
  std::cout << "Target length\t" << TargetLength << " cm" << std::endl;
  std::cout << "Delta Z\t" << DeltaZ << " cm" << std::endl;
  std::cout << "Xi\t" << Xi << std::endl;
  std::cout << "Output Spatial Pitch\t" << OutputXiPitch*DeltaZ << " cm" << std::endl;
  std::cout << "EdgeXi\t" << EdgeXi << std::endl;
  std::cout << "Simulation Time\t" << SimulationTime << " ns" << std::endl;
  std::cout << "Initial Time Pitch\t" << SimTimePitch*1e3 << " ps" << std::endl;
  std::cout << "Output Time Pitch\t" << OutputTimePitch << " ns" << std::endl;
  std::cout << "CFL number\t" << 30*SimTimePitch/DeltaZ << std::endl;
  std::cout << "Pressure\t" << Pressure << " Pa" << std::endl;
  std::cout << "Temperature\t" << Temperature << " degK" << std::endl;
  std::cout << "Number Density\t" << NumberDensity << " /cm3" << std::endl;
  std::cout << "Omega 0\t2pi * " << Omega_0/2/TMath::Pi() << " Grad/s" << std::endl;
  std::cout << "Detuning\t2pi * " << Detuning/2/TMath::Pi() << " Grad/s" << std::endl;
  std::cout << "Gamma1g\t2pi * " << Gamma_1g/2/TMath::Pi() << " Grad/s" << std::endl;
  std::cout << "Gamma1e\t2pi * " << Gamma_1e/2/TMath::Pi() << " Grad/s" << std::endl;
  std::cout << "Gamma2\t2pi * " << Gamma_2/2/TMath::Pi() << " Grad/s" << std::endl;
  std::cout << "ERpIntensity\t" << ERpIntensity << " GW/cm2" << std::endl;
  std::cout << "ELmIntensity\t" << ELmIntensity << " GW/cm2" << std::endl;
  std::cout << "ERpWidth (sigma)\t" << ERpWidth << " ns" << std::endl;
  std::cout << "ELmWidth (sigma)\t" << ELmWidth << " ns" << std::endl;
  std::cout << "ERpDelay\t" << ERpDelay << " ns" << std::endl;
  std::cout << "ELmDelay\t" << ELmDelay << " ns" << std::endl;
  std::cout << "TrigIntensity\t" << TrigIntensity << " GW/cm2" << std::endl;
  std::cout << "TrigWidth\t" << TrigWidth << " ns" << std::endl;
  std::cout << "TrigDelay\t" << TrigDelay << " ns" << std::endl;
  std::cout << "Resonance\t" << Resonance << std::endl;
  std::cout << "Absorption\t" << Absorption << std::endl;
  std::cout << "Relaxation\t" << Relaxation << std::endl;
  std::cout << "AbsError\t" << AbsError << std::endl;
  std::cout << "RelError\t" << RelError << std::endl;
  std::cout << "OpenMP\t" << openmp << std::endl;

  TPaveText* setting[2];
  setting[0] = new TPaveText( 0, 0, 0.5, 1 );
  setting[1] = new TPaveText( 0.5, 0, 1, 1 );
  for( int i=0; i<2; i++ ){
    setting[i]->SetFillStyle(0);
    setting[i]->SetTextAlign(12);
  }
  setting[0]->AddText( Form( "Target length \t %.1e mm", TargetLength*10 ) );
  setting[0]->AddText( Form( "Sim. Spatial Pitch \t %.2e mm", DeltaZ*10 ) );
  setting[0]->AddText( Form( "Output Spatial Pitch \t %.2e mm", OutputXiPitch*DeltaZ*10 ) );
  setting[0]->AddText( Form( "Spatial Divisions \t %d", 2*Xi ) );
  //setting[0]->AddText( Form( "Edge Divisions \t %d", EdgeXi) );
  setting[0]->AddText( Form( "Sim. Time Pitch \t %.2e ps", SimTimePitch*1e3 ) );
  setting[0]->AddText( Form( "Output Time Pitch \t %.2e ns", OutputTimePitch ) );
  setting[0]->AddText( Form( "CFL Number \t %.1f", 30.*SimTimePitch/DeltaZ) );
  setting[0]->AddText( Form( "Resonance \t %d", Resonance ) );
  setting[0]->AddText( Form( "Absorption \t %d", Absorption ) );
  setting[0]->AddText( Form( "Relaxation \t %d", Relaxation ) );
  //setting[0]->AddText( Form( "Square \t %d", Square ) );
  setting[1]->AddText( Form( "Driving peak intensity \t %.1f MW/cm2", ERpIntensity*1000) );
  //setting[0]->AddText( Form( "Driving(L) peak intensity \t %.2f MW/cm2", ELmIntensity*1000) );
  //setting[0]->AddText( Form( "Driving width (sigma) \t %.1f ns", ELmWidth) );
  setting[1]->AddText( Form( "Trig peak intensity \t %.1f MW/cm2", TrigIntensity*1000) );
  //setting[0]->AddText( Form( "Trig width (sigma) \t %.1f ns", TrigWidth) );
  setting[1]->AddText( Form( "Laser width (sigma) \t %.1f ns", ERpWidth) );
  setting[1]->AddText( Form( "Detuning \t 2#pi#times%.1f MHz", Detuning/2/TMath::Pi()*1000 ) );
  setting[1]->AddText( Form( "gamma1 \t 2#pi#times%.1f MHz", Gamma_1g/2/TMath::Pi()*1000 ) );
  setting[1]->AddText( Form( "gamma2 \t 2#pi#times%.1f MHz", Gamma_2/2/TMath::Pi()*1000 ) );
  setting[0]->AddText( Form( "Equation type \t %d", EquationType) );
  if( SpatialAlgorithm == 0 ){
    setting[0]->AddText( "Spatial Algorithm \t Central difference" );
  }else if( SpatialAlgorithm == 1 ){
    setting[0]->AddText( "Spatial Algorithm \t WENO" );
  }
  if( TimeAlgorithm == 0 ){
    setting[0]->AddText( "Time Algorithm \t Dopri5" );
  }else if( TimeAlgorithm == 1 ){
    setting[0]->AddText( "Time Algorithm \t Bulirsch-Stoer" );
  }
  setting[0]->AddText( Form( "Abs. error \t %.1e", AbsError ) );
  setting[0]->AddText( Form( "Rel. error \t %.1e", RelError ) );
  setting[0]->AddText( Form( "OpenMP \t %d" , openmp ) );
  

  // Header Tree branch preparation
  TFile* otf = new TFile( Form( "%s.root", filename), "RECREATE" );
  TTree* otr0 = new TTree( "header", "header" );
  otr0->Branch( "TargetLength", &TargetLength, "TargetLength/D" );
  otr0->Branch( "DeltaZ", &DeltaZ, "DeltaZ/D" );
  otr0->Branch( "Xi", &Xi, "Xi/I" );
  otr0->Branch( "OutputXiPitch", &OutputXiPitch, "OutputXiPitch/I" );
  otr0->Branch( "EdgeXi", &EdgeXi, "EdgeXi/I" );
  otr0->Branch( "SimulationTime", &SimulationTime, "SimulationTime/D" );
  otr0->Branch( "SimTimePitch", &SimTimePitch, "SimTimePitch/D" );
  otr0->Branch( "OutputTimePitch", &OutputTimePitch, "OutputTimePitch/D" );
  otr0->Branch( "Pressure", &Pressure, "Pressure/D" );
  otr0->Branch( "Temperature", &Temperature, "Temperature/D" );
  otr0->Branch( "NumberDensity", &NumberDensity, "NumberDensity/D" );
  otr0->Branch( "Omega_0", &Omega_0, "Omega_0/D" );
  otr0->Branch( "Detuning", &Detuning, "Detuning/D" );
  otr0->Branch( "Gamma_1g", &Gamma_1g, "Gamma_1g/D" );
  otr0->Branch( "Gamma_1e", &Gamma_1e, "Gamma_1e/D" );
  otr0->Branch( "Gamma_2", &Gamma_2, "Gamma_2/D" );
  otr0->Branch( "ERpIntensity", &ERpIntensity, "ERpIntensity/D" );
  otr0->Branch( "ELmIntensity", &ELmIntensity, "ELmIntensity/D" );
  otr0->Branch( "ERpWidth", &ERpWidth, "ERpWidth/D" );
  otr0->Branch( "ELmWidth", &ELmWidth, "ELmWidth/D" );
  otr0->Branch( "ERpDelay", &ERpDelay, "ERpDelay/D" );
  otr0->Branch( "ELmDelay", &ELmDelay, "ELmDelay/D" );
  otr0->Branch( "TrigIntensity", &TrigIntensity, "TrigIntensity/D" );
  otr0->Branch( "TrigWidth", &TrigWidth, "TrigWidth/D" );
  otr0->Branch( "TrigDelay", &TrigDelay, "TrigDelay/D" );
  otr0->Branch( "Resonance", &Resonance, "Resonance/O" );
  otr0->Branch( "Absorption", &Absorption, "Absorption/O" );
  otr0->Branch( "Relaxation", &Relaxation, "Relaxation/O" );
  otr0->Branch( "Square", &Square, "Square/O" );
  otr0->Branch( "EquationType", &EquationType, "EquationType/I" );
  otr0->Branch( "SpatialAlgorithm", &SpatialAlgorithm, "SpatialAlgorithm/I" );
  otr0->Branch( "TimeAlgorithm", &TimeAlgorithm, "TimeAlgorithm/I" );
  otr0->Branch( "AbsError", &AbsError, "AbsError/D" );
  otr0->Branch( "RelError", &RelError, "RelError/D" );
  otr0->Fill();
  otr0->Write();

  
  // Value variables preparation
  int NOutputXi = 2*(Xi/OutputXiPitch) + 3;
  double Time = 0;
  double (*Rhogg)[2] = new double [NOutputXi][2]();
  double (*Rhoee)[2] = new double [NOutputXi][2]();
  double (*Rhoge)[2] = new double [NOutputXi][2]();
  double (*ERp)[2] = new double [NOutputXi][2]();
  double (*ELm)[2] = new double [NOutputXi][2]();
  double (*ETrig)[2] = new double [NOutputXi][2]();
  double (*ESig)[2] = new double [NOutputXi][2]();
  TTree* otr1 = new TTree( "tr", "tr" );
  otr1->Branch( "Time", &Time, "Time/D" );
  otr1->Branch( "Rhogg", Rhogg, Form( "Rhogg[%d][2]/D", NOutputXi) );
  otr1->Branch( "Rhoee", Rhoee, Form( "Rhoee[%d][2]/D", NOutputXi) );
  otr1->Branch( "Rhoge", Rhoge, Form( "Rhoge[%d][2]/D", NOutputXi) );
  otr1->Branch( "ERp", ERp, Form( "ERp[%d][2]/D", NOutputXi) );
  otr1->Branch( "ELm", ELm, Form( "ELm[%d][2]/D", NOutputXi) );
  otr1->Branch( "ETrig", ETrig, Form( "ETrig[%d][2]/D", NOutputXi) );
  otr1->Branch( "ESig", ESig, Form( "ESig[%d][2]/D", NOutputXi) );




  // Histograms preparation for c1
  TH2D* h[12];
  for( int i=0; i<12; i++ ){
    h[i] = new TH2D( Form( "h%d", i), Form( "h%d;z [cm];time [ns];value", i), 201, -TargetLength/2., TargetLength/2., 100, 0-OutputTimePitch/2, SimulationTime-OutputTimePitch/2 );
    h[i]->SetStats( false );
  }
  h[0]->SetTitle( "#rho_ee" );
  h[1]->SetTitle( "Abs( #rho_ge )" );
  h[2]->SetTitle( "Phase( #rho_ge )" );
  h[3]->SetTitle( "Intensity( Trig. )" );
  h[4]->SetTitle( "Intensity( Sig. )" );
  h[5]->SetTitle( "Phase( Trig. )" );
  h[6]->SetTitle( "Intensity( ER+ )" );
  h[7]->SetTitle( "Intensity( EL- )" );
  h[8]->SetTitle( "Phase( Sig. )" );
  h[9]->SetTitle( "Phase( ER+ )" );
  h[10]->SetTitle( "Phase( EL- )" );
  h[11]->SetTitle( "#rho_gg + #rho_ee" );

  // Graphs preparation for c2
  const int nTime = SimulationTime/OutputTimePitch;
  const int nPos = 2*(Xi/OutputXiPitch) + 1;

  TGraph* gField[3][4]; // (left, right, envelop) * ( ERp. ELm. Trig, Sig )
  TGraph* gBloch[2][2];  // (center, coherence max) * ( rho_ee, |rho_ge| )
  for( int i=0; i<4; i++ ){
    gField[0][i] = new TGraph( nTime ); // Left side
    gField[1][i] = new TGraph( nTime ); // Right side
    gField[2][i] = new TGraph( nPos );  // Envelop
    for( int j=0; j<3; j++ ){
      gField[j][i]->SetName( Form( "gField%d%d", j, i ) );
      gField[j][i]->SetLineColor( i+2 );
      gField[j][i]->SetLineWidth(2);
    }
  }
  for( int i=0; i<2; i++ ){
    gBloch[0][i] = new TGraph( nTime ); // Center
    gBloch[1][i] = new TGraph( nPos );  // Coherence max
    for( int j=0; j<2; j++ ){
      gBloch[j][i]->SetTitle( Form("gBloch%d%d", j, i) );
      gBloch[j][i]->SetLineColor( i+6 );
      gBloch[j][i]->SetLineWidth(2);
    }
  }

  TMultiGraph* gMultiField[3];
  gMultiField[0] = new TMultiGraph( "gMultiField0", "Left Side Intensity;time [ns];intensity [GW/cm2]" );
  gMultiField[1] = new TMultiGraph( "gMultiField1", "Right Side Intensity;time [ns];intensity [GW/cm2]" );
  gMultiField[2] = new TMultiGraph( "gMultiField2", "Power;position [cm];power [J/cm2]" );
  TMultiGraph* gMultiBloch[2];
  gMultiBloch[0] = new TMultiGraph( "gMultiBloch0", "Density matrix element at the central;time [ns];" );
  gMultiBloch[1] = new TMultiGraph( "gMultiBloch1", "Max. Density matrix element;position [cm];" );




  // File reading & histogram filling
  const int tInterval = (SimulationTime/OutputTimePitch)/100;
  const int xInterval = Xi/OutputXiPitch/100;
  int tndex = 0;
  double intensities[4] = {0};
  double coherence = 0;
  double pos = 0;
  while( ifs >> Time ){
    for( int xi=0; xi<NOutputXi; xi++ ){
      ifs >> Rhogg[xi][0] >> Rhogg[xi][1] >> Rhoee[xi][0] >> Rhoee[xi][1] >> Rhoge[xi][0] >> Rhoge[xi][1];
      ifs >> ERp[xi][0] >> ERp[xi][1] >> ELm[xi][0] >> ELm[xi][1] >> ETrig[xi][0] >> ETrig[xi][1] >> ESig[xi][0] >> ESig[xi][1];

      intensities[0] = ( pow(ERp[xi][0],2) + pow(ERp[xi][1],2) ) / pow(8.68021098E5,2);
      intensities[1] = ( pow(ELm[xi][0],2) + pow(ELm[xi][1],2) ) / pow(8.68021098E5,2);
      intensities[2] = ( pow(ETrig[xi][0],2) + pow(ETrig[xi][1],2) ) / pow(8.68021098E5,2);
      intensities[3] = ( pow(ESig[xi][0],2) + pow(ESig[xi][1],2) ) / pow(8.68021098E5,2);
      coherence = sqrt( pow(Rhoge[xi][0],2) + pow(Rhoge[xi][1],2) );
      int x = xi-1;

      // Edge region
      if( x < 0 ) continue;
      if( x > 2*Xi/OutputXiPitch ) continue;

      pos = x*DeltaZ*OutputXiPitch - TargetLength/2.;

      
      // Fill
      if( tndex%tInterval == 0 ){
	if( x%xInterval == 0 ){
	  if( h[0]->GetBinContent( h[0]->FindBin( pos, Time) ) == 0 ) {
	    h[0]->Fill( pos, Time, Rhoee[xi][0] );
	    h[1]->Fill( pos, Time, coherence );
	    h[2]->Fill( pos, Time, GetArg( Rhoge[xi][0], Rhoge[xi][1] ) );
	    h[3]->Fill( pos, Time, intensities[2] );
	    h[4]->Fill( pos, Time, intensities[3] );
	    h[5]->Fill( pos, Time, GetArg( ETrig[xi][0],ETrig[xi][1] ) );
	    h[6]->Fill( pos, Time, intensities[0] );
	    h[7]->Fill( pos, Time, intensities[1] );
	    h[8]->Fill( pos, Time, GetArg( ESig[xi][0],ESig[xi][1] ) );
	    h[9]->Fill( pos, Time, GetArg( ERp[xi][0],ERp[xi][1] ) );
	    h[10]->Fill( pos, Time, GetArg( ELm[xi][0],ELm[xi][1] ) );
	    h[11]->Fill( pos, Time, Rhogg[xi][0]+Rhoee[xi][0] );
	  }
	}
      }

      // graph
      if( x == 0 ){ // Left side
	for( int i=0; i<4; i++ ){
	  gField[0][i]->SetPoint( tndex, Time, intensities[i] );
	}
      }
      if( x == 2*Xi/OutputXiPitch ){ // Right side
	for( int i=0; i<4; i++ ){
	  gField[1][i]->SetPoint( tndex, Time, intensities[i] );
	}
      }
      double xp=0, yp=0;
      for( int i=0; i<4; i++ ){ // Integral
	gField[2][i]->GetPoint( x, xp, yp );
	gField[2][i]->SetPoint( x, pos, yp+intensities[i]*OutputTimePitch );
      }
      if( x == Xi/OutputXiPitch ){ // Center
	gBloch[0][0]->SetPoint( tndex, Time, Rhoee[xi][0] );
	gBloch[0][1]->SetPoint( tndex, Time, coherence );
      }
      { // maximum excitation and coherence
	gBloch[1][0]->GetPoint( x, xp, yp );
	if( yp < Rhoee[xi][0] ){
	  gBloch[1][0]->SetPoint( x, pos, Rhoee[xi][0] );
	}
	gBloch[1][1]->GetPoint( x, xp, yp );
	if( yp < coherence ){
	  gBloch[1][1]->SetPoint( x, pos, coherence );
	}
      }

    }
    otr1->Fill();
    tndex++;
  }


  // histogram drawing
  TCanvas* c1 = new TCanvas( "c1", "c1", 1600, 1200 );
  c1->Divide( 3, 4, 0.001, 0.001 );
  for( int ho=0; ho<3; ho++ ){
    for( int v=0; v<4; v++ ){
      c1->cd( ho+v*3+1 );
      h[ho+v*3]->Draw( "colz" );
    }
  }



  // result summary
  double MaxExcitation = TMath::MaxElement( gBloch[1][0]->GetN(), gBloch[1][0]->GetY() );
  double MaxCoherence = TMath::MaxElement( gBloch[1][1]->GetN(), gBloch[1][1]->GetY() );
  double Power[4] = {0}, dummy = 0;
  for( int i=0; i<4; i++ ){
    gField[2][i]->GetPoint( 0, dummy, Power[i] );
  }
  setting[1]->AddText( Form( "Driving Energy Density \t %.1e J/cm2", Power[0]) );
  setting[1]->AddText( Form( "Trigger Energy Density \t %.1e J/cm2", Power[2]) );
  setting[1]->AddText( "" );
  setting[1]->AddText( Form( "Max. excitation \t %.2e", MaxExcitation ) );
  setting[1]->AddText( Form( "Max. coherence \t %.2e", MaxCoherence ) );
  setting[1]->AddText( Form( "Signal Energy Density \t %.2e J/cm2", Power[3] ) );




  // graph scaling
  int scaleOrder[3] = {0};
  scaleOrder[0] = floor(log10(Power[0]/Power[2]));
  scaleOrder[1] = floor(log10(Power[0]/Power[3]));
  double xp=0, yp=0;
  for( int i=2; i<4; i++ ){
    for( int j=0; j<3; j++ ){
      for( int p=0; p<gField[j][i]->GetN(); p++ ){
	gField[j][i]->GetPoint( p, xp, yp );
	gField[j][i]->SetPoint( p, xp, yp*pow(10,scaleOrder[i-2]) );
      }
    }
  }
  scaleOrder[2] = floor(log10(MaxCoherence/MaxExcitation) );
  for( int i=0; i<2; i++ ){
    for( int p=0; p<gBloch[i][0]->GetN(); p++ ){
      gBloch[i][0]->GetPoint( p, xp, yp );
      gBloch[i][0]->SetPoint( p, xp, yp*pow(10,scaleOrder[2]) );
    }
  }



  TLegend* legend[2];
  legend[0] = new TLegend( 0.6, 0.6, 0.9, 0.9, "Field" );
  legend[1] = new TLegend( 0.7, 0.7, 0.9, 0.9, "Density Matrix" );
  for( int i=0; i<2; i++ ){
    legend[i]->SetFillColor(0);
  }
  legend[0]->AddEntry( gField[0][0], "Driving (R+)", "L" );
  legend[0]->AddEntry( gField[0][1], "Driving (L-)", "L" );
  legend[0]->AddEntry( gField[0][2], Form( "Trigger (R-) #times 1E%d", scaleOrder[0]), "L" );
  legend[0]->AddEntry( gField[0][3], Form( "Signal (L+) #times 1E%d", scaleOrder[1]), "L" );
  legend[1]->AddEntry( gBloch[0][0], Form( "rho_ee #times 1E%d", scaleOrder[2]), "L" );
  legend[1]->AddEntry( gBloch[0][1], "|rho_ge|", "L" );






	  

  // graph drawing
  for( int i=0; i<4; i++ ){
    for( int j=0; j<3; j++ ){
      gMultiField[j]->Add( gField[j][i] );
    }
  }
  for( int i=0; i<2; i++ ){
    for( int j=0; j<2; j++ ){
      gMultiBloch[j]->Add( gBloch[j][i] );
    }
  }
  TCanvas* c2 = new TCanvas( "c2", "c2", 1600, 1200 );
  c2->Divide( 2, 3, 0.001, 0.001 );
  for( int i=0; i<2; i++ ){
    c2->cd( i+1 );
    gMultiBloch[i]->Draw( "AL" );
  }
  legend[1]->Draw();
  for( int i=0; i<3; i++ ){
    c2->cd( i+3 );
    gMultiField[i*2-i/2*3]->Draw( "AL" );
  }
  legend[0]->Draw();



  // result calculation
  
  c2->cd(6);
  setting[0]->Draw();
  setting[1]->Draw();


  otr1->Write();
  c1->Write();
  c2->Write();
  for( int i=0; i<12; i++ ){
    h[i]->Write();
  }
  for( int i=0; i<4; i++ ){
    for( int j=0; j<3; j++ ){
      gField[j][i]->Write();
    }
  }
  for( int i=0; i<2; i++ ){
    for( int j=0; j<2; j++ ){
      gBloch[j][i]->Write();
    }
  }


  c1->Print( Form( "%s_0.png", filename) );
  c2->Print( Form( "%s_1.png", filename) );

  //otf->Close();

  return 0;


}

  




