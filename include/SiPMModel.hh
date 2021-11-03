/// \file  SiPMModel.hh
/// \brief Definition of SiPM models charcacteristics

#ifndef SiPMModel_h
#define SiPMModel_h 1

#include "globals.hh"

namespace Model {
	static constexpr int NbPixelsY[3] = {15, 23, 46}; // pitch = 75, 50, 25 mum
	static constexpr int NbPixelsX[3] = {19, 29, 58}; // pitch = 75, 50, 25 mum

	static constexpr double FillFactor[3] = {0.82, 0.74, 0.47}; // pitch = 75, 50, 25 mum

	static constexpr double Gain[3] = {4e6, 1.7e6, 7e5}; // pitch = 75, 50, 25 mum
	static constexpr double OVoltage[3] = {3, 3, 5}; // pitch = 75, 50, 25 mum

	static constexpr double dark_noise_rate[3] = {90*CLHEP::kilohertz, 90*CLHEP::kilohertz, 70*CLHEP::kilohertz}; // pitch = 75, 50, 25 mum
	
	static constexpr double r_index[2] = {1.41, 1.55}; //Window material = Si resin, Epoxy resin
	static constexpr double SiPM_size_Z[2] = {1*CLHEP::mm, .55*CLHEP::mm}; //Window material = Si resin, Epoxy resin
	static constexpr double window_size_Z[2] = {0.5*CLHEP::mm, 0.3*CLHEP::mm}; //Window material = Si resin, Epoxy resin


	static constexpr char eff_name[][23] =  {"SiPM_det_eff_75_CS.txt", 
						 "SiPM_det_eff_50_CS.txt",
						 "SiPM_det_eff_25_CS.txt",
						 "SiPM_det_eff_75_PE.txt", 
						 "SiPM_det_eff_50_PE.txt",
						 "SiPM_det_eff_25_PE.txt"};

	static constexpr char eff_gain_name[][32] =  {"SiPM_photon_det_eff_gain_75.txt", 
						      "SiPM_photon_det_eff_gain_50.txt",
						      "SiPM_photon_det_eff_gain_25.txt"};

	static constexpr char OCT_gain_name[][21] =  {"SiPM_OCT_gain_75.txt", 
					   	      "SiPM_OCT_gain_50.txt",
						      "SiPM_OCT_gain_25.txt"};

	static constexpr double OCT_factor[3] = {3.3904690, 5.2915961, 12.451933};

	static constexpr char photon_gain_name[][26] =  {"SiPM_photon_gain_75.txt", 
					   	         "SiPM_photon_gain_50.txt",
						         "SiPM_photon_gain_25.txt"};
}

#endif


