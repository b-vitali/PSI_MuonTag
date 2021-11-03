/// \file  PixelHit.hh
/// \brief Definition of the PixelHit class

#ifndef PixelHit_h
#define PixelHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

#include "tls.hh"

#include<vector>

class PixelHit : public G4VHit{
	public:
		PixelHit();
		virtual ~PixelHit();

		PixelHit(const PixelHit &right);
		const PixelHit& operator=(const PixelHit &right);
		G4bool operator==(const PixelHit &right) const;

		inline void *operator new(size_t);
		inline void operator delete(void *aHit);

		virtual void Draw();
		virtual void Print();

		inline void SetNCells(G4int val){fNCells = val;}
		inline void SetNPhotoElectrons(G4double val){fNPhotoElectrons = val;}
		inline void SetCells(std::vector<G4int> val){fCells = val;}
		inline void SetCellTime(std::vector<G4double> val){fCellTime = val;}
		inline void SetOCTFlag(std::vector<G4int> val){fOCTflag = val;}
		inline void SetDNFlag(std::vector<G4int> val){fDNflag = val;}

		inline G4int GetNCells(){return fNCells;}
		inline G4double GetNPhotoElectrons(){return fNPhotoElectrons;}
		inline std::vector<G4int> GetCells(){return fCells;}
		inline std::vector<G4double> GetCellTime(){return fCellTime;}
		inline std::vector<G4int> GetOCTFlag(){return fOCTflag;}
		inline std::vector<G4int> GetDNFlag(){return fDNflag;}

		inline void Clear(){fNCells = 0; fNPhotoElectrons = 0; fCells.clear(); 
			fCellTime.clear(); fPhysVol = nullptr; fDrawit = false; fPixelNumber = -1; 
			fOCTflag.clear(); fDNflag.clear();
		}

		inline void SetPixelPhysVol(G4VPhysicalVolume* physVol){this->fPhysVol = physVol;}
		inline G4VPhysicalVolume* GetPixelPhysVol(){return fPhysVol;}

		inline void SetPixelPhysVolMother(G4VPhysicalVolume* physVol){this->fPhysVolMother = physVol;}
		inline G4VPhysicalVolume* GetPixelPhysVolMother(){return fPhysVolMother;}

		inline void SetPixelPhysVolGMother(G4VPhysicalVolume* physVol){this->fPhysVolGMother = physVol;}
		inline G4VPhysicalVolume* GetPixelPhysVolGMother(){return fPhysVolGMother;}

		inline void SetDrawit(G4bool b){fDrawit=b;}
		inline G4bool GetDrawit(){return fDrawit;}

		inline void SetPixelNumber(G4int val){fPixelNumber = val;}
		inline G4int GetPixelNumber(){return fPixelNumber;}

	private:
		G4int fNCells;
		G4double fNPhotoElectrons;
		std::vector<G4int> fCells;
		std::vector<G4double> fCellTime;
		std::vector<G4int> fOCTflag;
		std::vector<G4int> fDNflag;
		G4VPhysicalVolume* fPhysVol;
		G4VPhysicalVolume* fPhysVolMother;
		G4VPhysicalVolume* fPhysVolGMother;
		G4bool fDrawit;
		G4int fPixelNumber;
};

typedef G4THitsCollection<PixelHit> PixelHitsCollection;

extern G4ThreadLocal G4Allocator<PixelHit>* PixelHitAllocator;

inline void* PixelHit::operator new(size_t){
	if(!PixelHitAllocator) PixelHitAllocator = new G4Allocator<PixelHit>;
	return (void*) PixelHitAllocator->MallocSingle();
}

inline void PixelHit::operator delete(void* aHit){
	PixelHitAllocator->FreeSingle((PixelHit*) aHit);
}

#endif


