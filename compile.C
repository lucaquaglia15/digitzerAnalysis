//
//TRACKS - MC generation of particle tracks in a detector
//with multiple scattering and noise
//developed by Luca Quaglia and Mattia Ivaldi, 2018
//
//START

void compile(TString myopt="fast"){

  printf("\n\n+++ DIGITIZER START +++\n\n");

  TString opt;

  if(myopt.Contains("force")){
    opt = "kfg";
  }
  else {
    opt = "kg";
  }

  gSystem->CompileMacro("PeakFinder.cpp",opt.Data());//compile classes and libraries
  gSystem->CompileMacro("SGSmooth.cpp",opt.Data());//compile classes and libraries
  gSystem->CompileMacro("utilities.cxx",opt.Data());//compile classes and libraries
  gSystem->CompileMacro("analyzeDigitizer.C",opt.Data());//compile classes and libraries
  gROOT->ProcessLine("analyzeDigitizer(221,false)");//perform analysis of data
  
  printf("\n\n+++ DIGITIZER END +++\n\n");
}  
