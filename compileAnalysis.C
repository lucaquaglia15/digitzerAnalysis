//
//TRACKS - MC generation of particle tracks in a detector
//with multiple scattering and noise
//developed by Luca Quaglia and Mattia Ivaldi, 2018
//
//START

void compileAnalysis(TString myopt="fast"){

  printf("\n\n+++ RUN ANALYSIS START +++\n\n");

  TString opt;

  if(myopt.Contains("force")){
    opt = "kfg";
  }
  else {
    opt = "kg";
  }

  gSystem->CompileMacro("analysis.cxx",opt.Data());//compile classes and libraries
  gSystem->CompileMacro("plotter.C",opt.Data());//compile classes and libraries
  gROOT->ProcessLine("plotter(244)");//perform analysis of data
  
  printf("\n\n+++ RUN ANALYSIS END +++\n\n");
}  
