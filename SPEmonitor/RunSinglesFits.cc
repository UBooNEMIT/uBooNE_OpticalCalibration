#include<ctime>
#include<iostream>
#include<iomanip>
#include<fstream>
#include "TROOT.h"
#include "TRint.h"
#include <vector>

void RunSinglesFits()
{
  gROOT->ProcessLine(".L /uboone/app/users/moon/OpticalStudies/v05_01_01/workdir/SPEmonitor/PoissonFit.C");

  string line;
  std::vector<string> SPE_filepaths;
  ifstream myfile ("/uboone/app/users/moon/OpticalStudies/v05_01_01/workdir/SPEmonitor/SPE_TreePath.txt");
  if (myfile.is_open())
    {
      while ( getline (myfile,line) )
	{
	  SPE_filepaths.push_back(line);
	}
      myfile.close();
    }

  int CHANNELS = 32;

  std::ofstream output_file1("/uboone/app/users/moon/OpticalStudies/v05_01_01/workdir/SPEmonitor/CurrentSPEamps/first_recent_singles.txt");
  std::ofstream output_file2("/uboone/app/users/moon/OpticalStudies/v05_01_01/workdir/SPEmonitor/CurrentSPEamps/second_recent_singles.txt");
  std::ofstream output_file3("/uboone/app/users/moon/OpticalStudies/v05_01_01/workdir/SPEmonitor/CurrentSPEamps/third_recent_singles.txt");
  std::ofstream output_file4("/uboone/app/users/moon/OpticalStudies/v05_01_01/workdir/SPEmonitor/CurrentSPEamps/fourth_recent_singles.txt");

  std::ofstream output_file5("/uboone/app/users/moon/OpticalStudies/v05_01_01/workdir/SPEmonitor/CurrentSPEamps/first_recent_singles_sig.txt");
  std::ofstream output_file6("/uboone/app/users/moon/OpticalStudies/v05_01_01/workdir/SPEmonitor/CurrentSPEamps/second_recent_singles_sig.txt");
  std::ofstream output_file7("/uboone/app/users/moon/OpticalStudies/v05_01_01/workdir/SPEmonitor/CurrentSPEamps/third_recent_singles_sig.txt");
  std::ofstream output_file8("/uboone/app/users/moon/OpticalStudies/v05_01_01/workdir/SPEmonitor/CurrentSPEamps/fourth_recent_singles_sig.txt");

  int current_channel = 0;

  std::vector<double> SPE;

  while (current_channel < CHANNELS) {

    SPE = PoissonFit(SPE_filepaths[0].c_str(),current_channel);

    output_file1 << SPE[0] << "\n";
    output_file5 << SPE[1] << "\n";

    current_channel++;

  }
 
  current_channel = 0;
  while (current_channel < CHANNELS) {

    SPE = PoissonFit(SPE_filepaths[1].c_str(),current_channel);

    output_file2 << SPE[0] << "\n";
    output_file6 << SPE[1] << "\n";

    current_channel++;
  }

  current_channel = 0;
  while (current_channel < CHANNELS) {

    SPE = PoissonFit(SPE_filepaths[2].c_str(),current_channel);

    output_file3 << SPE[0] << "\n";
    output_file7 << SPE[1] << "\n";

    current_channel++;
  }

  current_channel = 0;
  while (current_channel < CHANNELS) {

    SPE = PoissonFit(SPE_filepaths[3].c_str(),current_channel);

    output_file4 << SPE[0] << "\n";
    output_file8 << SPE[1] << "\n";

    current_channel++;
  }

}




