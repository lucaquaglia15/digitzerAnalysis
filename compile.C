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

  gSystem->CompileMacro("utilities.cxx",opt.Data());//compile classes and libraries
  gSystem->CompileMacro("analyzeDigitizer.C",opt.Data());//compile classes and libraries
  gROOT->ProcessLine("analyzeDigitizer(182)");//perform generation
  
  printf("\n\n+++ DIGITIZER END +++\n\n");
}  
