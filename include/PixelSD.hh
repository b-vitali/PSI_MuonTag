/// \file  PixelSD.hh
/// \brief Definition of the PixelSD class

#ifndef PixelSD_h
#define PixelSD_h 1

#include "PixelHit.hh"
#include "SiPMModel.hh"

#include "G4VSensitiveDetector.hh"
#include "G4OpticalPhoton.hh"

class G4Step;
class G4HCofThisEvent;
class G4VLogicalVolume;
class PixelSDMessenger;

class PixelSD : public G4VSensitiveDetector{
	public:
		PixelSD(G4String name, G4String model);
		virtual ~PixelSD();

		virtual void Initialize(G4HCofThisEvent*);
		virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);
		virtual void EndOfEvent(G4HCofThisEvent*);
		virtual void clear();
		virtual void DrawAll();
		virtual void PrintAll();

		G4double GetAbsProbability(G4double);
		void SetVoltage(G4double val){fVoltage = val;}
		void SetDetEffGain();
		void SetPhotonGain();
		void SetOCT();

		
		inline G4String GetFileName(G4int i){return filename[i];}
		inline void SetNbOfPixels(G4int val){fNbOfPixels = val;}
		inline G4int GetNbOfPixels(){return fNbOfPixels;}

		void DefineProperties();
		void SetModel(G4String name){fModel = name;}

	private:
		PixelSDMessenger* fPixelMessenger;

		PixelHitsCollection* fPixelCollection;
		PixelHitsCollection* fPixelCollectionDraw;
		G4int fNCells, fNbOfPixels;
		G4double fNPhotoElectrons;
		std::vector<G4int> fCells;
		std::vector<G4double> fCellTime;
		G4int fOCTflag;
		std::vector<G4int> fOCTflagvec;
		std::vector<G4int> fDNflagvec;

		G4double fVoltage;
		std::vector<G4double> fDetEff;
		std::vector<G4double> fDetEffx;
		G4double fDetEffGain;
		G4double fPhotonGain;
		G4double fOCT;

		G4double fFillFactor;
		
		std::vector<G4int> firstStep;
		std::vector<G4double> activationTime;
		std::vector<G4int> isParent;

		G4bool fCmdOCT;
		G4bool fCmdDN;

		G4String fModel;
		G4String filename[4];
};

#endif


