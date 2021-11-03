void ordering(TString name, TString file, bool darkNoise){
	TFile* F = TFile::Open(name);
	TTree* T = (TTree*) F->Get("T");

	vector<int>* Cells = 0;
	vector<int>* DNflag = 0;
	vector<int>* Channel = 0;
	vector<double>* CellTime = 0;

	vector<int> CellsAll;
	vector<int> DNflagAll;
	vector<int> ChannelAll;
	vector<double> CellTimeAll;
	vector<double> CellTimeEvAll;
	vector<double> TrackLengthAll;
	vector<double> ThetaInAll;
	vector<int> BounceAll;
	vector<int> SecondaryIDAll;
	vector<int> eventIDAll;
	vector<int> SurfInAll;
	
	double ein, GunTime, TrackLength, ThetaIn;
	int eventID = 0, SurfIn = 0, NCells = 0, Bounce = 0, SecondaryID = 0;

	T->SetBranchAddress("ein", &ein);
	T->SetBranchAddress("Cells", &Cells);
	T->SetBranchAddress("DNflag", &DNflag);
	T->SetBranchAddress("CellTime", &CellTime);
	T->SetBranchAddress("eventID", &eventID);
	T->SetBranchAddress("SurfIn", &SurfIn);
	T->SetBranchAddress("GunTime", &GunTime);
	T->SetBranchAddress("TrackLength", &TrackLength);
	T->SetBranchAddress("ThetaIn", &ThetaIn);
	T->SetBranchAddress("SecondaryID", &SecondaryID);

	int eventsPerCycle = 100;
	int Ni = (T->GetEntries() / eventsPerCycle);
	if(T->GetEntries()%eventsPerCycle == 0) Ni -= 1;
	int nHalf = 0;
	vector<int> tempCells;
	vector<int> tempDNflag;
	vector<int> tempChannel;
	vector<double> tempCellTime;
	vector<double> tempCellTimeEv;
	vector<double> tempTrackLength;
	vector<double> tempThetaIn;
	vector<int> tempSecondaryID;
	vector<int> tempBounce;
	vector<int> tempeventID;
	vector<int> tempSurfIn;

	TFile* N = TFile::Open(file + ".root", "RECREATE");
	TTree* Tnew = new TTree("T", "signals");
	N->cd();
	int sCells = 0, sDNflag = 0, sChannel = 0, seventID = 0, sSurfIn = 0, sBounce = 0, sSecondaryID = 0;
	double sCellTime = 0, sCellTimeEv = 0, sTrackLength = 0, sThetaIn = 0;

	Tnew->Branch("Cells", &sCells);
	Tnew->Branch("DNflag", &sDNflag);
	Tnew->Branch("CellTime", &sCellTime);
	Tnew->Branch("CellTimeEv", &sCellTimeEv); //Time since GunTime
	Tnew->Branch("TrackLength", &sTrackLength);
	Tnew->Branch("ThetaIn", &sThetaIn);
	Tnew->Branch("SecondaryID", &sSecondaryID);
	Tnew->Branch("Bounce", &sBounce);
	Tnew->Branch("eventID", &seventID);
	Tnew->Branch("SurfIn", &sSurfIn);

	for(int i = 0; i < Ni; i++){
		for(int k = i * eventsPerCycle; k < eventsPerCycle * (i + 1); k++){
			T->GetEntry(k);
			int n = Cells->size();
			if(n > 0){
				for(int j = 0; j < n; j++){
					if((darkNoise || DNflag->at(j) == 0)){
						CellsAll.push_back(Cells->at(j));
						DNflagAll.push_back(DNflag->at(j));
						CellTimeAll.push_back(CellTime->at(j));
						CellTimeEvAll.push_back(CellTime->at(j) - GunTime);
						TrackLengthAll.push_back(TrackLength);
						ThetaInAll.push_back(ThetaIn);
						SecondaryIDAll.push_back(SecondaryID);
						BounceAll.push_back(Bounce);
						eventIDAll.push_back(eventID);
						SurfInAll.push_back(SurfIn);
					}
				}
			}
		}
		Cells->clear();
		DNflag->clear();
		CellTime->clear();

		int n = CellsAll.size();
		double CellTimeL[n];
		int I[n];

		for(int i = 0; i < n; i++){
			CellTimeL[i] = CellTimeAll.at(i);
		}

		CellTimeAll.clear();
		TMath::Sort(n, CellTimeL, I, false);
		
		int tempnHalf = n/2;
		if(i == 0){
			for(int j = 0; j < tempnHalf; j++){
				sCells = CellsAll.at(I[j]);
				sDNflag = DNflagAll.at(I[j]);
				sCellTime = CellTimeL[I[j]];
				sCellTimeEv = CellTimeEvAll.at(I[j]);
				sTrackLength = TrackLengthAll.at(I[j]);
				sThetaIn = ThetaInAll.at(I[j]);
				sSecondaryID = SecondaryIDAll.at(I[j]);
				sBounce = BounceAll.at(I[j]);
				seventID = eventIDAll.at(I[j]);
				sSurfIn = SurfInAll.at(I[j]);
				Tnew->Fill();
			}

			for(int j = tempnHalf; j < n; j++){
				tempCells.push_back(CellsAll.at(I[j]));
				tempDNflag.push_back(DNflagAll.at(I[j]));
				tempCellTime.push_back(CellTimeL[I[j]]);
				tempCellTimeEv.push_back(CellTimeEvAll.at(I[j]));
				tempTrackLength.push_back(TrackLengthAll.at(I[j]));
				tempThetaIn.push_back(ThetaInAll.at(I[j]));
				tempSecondaryID.push_back(SecondaryIDAll.at(I[j]));
				tempBounce.push_back(BounceAll.at(I[j]));
				tempeventID.push_back(eventIDAll.at(I[j]));
				tempSurfIn.push_back(SurfInAll.at(I[j]));
			}
			nHalf = n - tempnHalf;
		}

		else{
			for(int j = 0; j < tempnHalf; j++){
				tempCells.push_back(CellsAll.at(I[j]));
				tempDNflag.push_back(DNflagAll.at(I[j]));
				tempCellTime.push_back(CellTimeL[I[j]]);
				tempCellTimeEv.push_back(CellTimeEvAll.at(I[j]));
				tempTrackLength.push_back(TrackLengthAll.at(I[j]));
				tempThetaIn.push_back(ThetaInAll.at(I[j]));
				tempSecondaryID.push_back(SecondaryIDAll.at(I[j]));
				tempBounce.push_back(BounceAll.at(I[j]));
				tempeventID.push_back(eventIDAll.at(I[j]));
				tempSurfIn.push_back(SurfInAll.at(I[j]));
			}

			int ntemp = tempCellTime.size();
			double tempCellTimeL[ntemp];
			int tempI[ntemp];
			for(int j = 0; j < ntemp; j++){
				tempCellTimeL[j] = tempCellTime.at(j);
			}
			tempCellTime.clear();

			TMath::Sort(ntemp, tempCellTimeL, tempI, false);

			for(int j = 0; j < ntemp; j++){
				sCells = tempCells.at(tempI[j]);
				sDNflag = tempDNflag.at(tempI[j]);
				sCellTime = tempCellTimeL[tempI[j]];
				sCellTimeEv = tempCellTimeEv.at(tempI[j]);
				sTrackLength = tempTrackLength.at(tempI[j]);
				sThetaIn = tempThetaIn.at(tempI[j]);
				sSecondaryID = tempSecondaryID.at(tempI[j]);
				sBounce = tempBounce.at(tempI[j]);
				sSurfIn = tempSurfIn.at(tempI[j]);
				seventID = tempeventID.at(tempI[j]);
				Tnew->Fill();
			}
			tempCells.clear();
			tempDNflag.clear();
			tempCellTimeEv.clear();
			tempTrackLength.clear();
			tempThetaIn.clear();
			tempSecondaryID.clear();
			tempBounce.clear();
			tempSurfIn.clear();
			tempeventID.clear();

			nHalf = n - tempnHalf;
			for(int j = tempnHalf; j < n; j++){
				tempCells.push_back(CellsAll.at(I[j]));
				tempDNflag.push_back(DNflagAll.at(I[j]));
				tempCellTime.push_back(CellTimeL[I[j]]);
				tempCellTimeEv.push_back(CellTimeEvAll.at(I[j]));
				tempTrackLength.push_back(TrackLengthAll.at(I[j]));
				tempThetaIn.push_back(ThetaInAll.at(I[j]));
				tempSecondaryID.push_back(SecondaryIDAll.at(I[j]));
				tempBounce.push_back(BounceAll.at(I[j]));
				tempeventID.push_back(eventIDAll.at(I[j]));
				tempSurfIn.push_back(SurfInAll.at(I[j]));
			}
		}
		CellsAll.clear();
		DNflagAll.clear();
		CellTimeEvAll.clear();
		TrackLengthAll.clear();
		ThetaInAll.clear();
		SecondaryIDAll.clear();
		BounceAll.clear();
		eventIDAll.clear();
		SurfInAll.clear();
	}

	for(int i = Ni*eventsPerCycle; i < T->GetEntries(); ++i){
		T->GetEntry(i);
		int n = Cells->size();
		if(n > 0){
			for(int j = 0; j < n; j++){
				if((darkNoise || DNflag->at(j) == 0) && ein > 0){
					CellsAll.push_back(Cells->at(j));
					DNflagAll.push_back(DNflag->at(j));
					CellTimeAll.push_back(CellTime->at(j));
					CellTimeEvAll.push_back(CellTime->at(j) - GunTime);
					TrackLengthAll.push_back(TrackLength);
					ThetaInAll.push_back(ThetaIn);
					SecondaryIDAll.push_back(SecondaryID);
					BounceAll.push_back(Bounce);
					eventIDAll.push_back(eventID);
					SurfInAll.push_back(SurfIn);
				}
			}
		}
	}

	Cells->clear();
	DNflag->clear();
	CellTime->clear();

	int n = CellsAll.size();
	double CellTimeL[n];
	int I[n];

	for(int i = 0; i < n; i++){
		CellTimeL[i] = CellTimeAll.at(i);
	}

	CellTimeAll.clear();
	TMath::Sort(n, CellTimeL, I, false);
	
	for(int j = 0; j < n; j++){
		tempCells.push_back(CellsAll.at(I[j]));
		tempDNflag.push_back(DNflagAll.at(I[j]));
		tempCellTime.push_back(CellTimeL[I[j]]);
		tempCellTimeEv.push_back(CellTimeEvAll.at(I[j]));
		tempTrackLength.push_back(TrackLengthAll.at(I[j]));
		tempThetaIn.push_back(ThetaInAll.at(I[j]));
		tempSecondaryID.push_back(SecondaryIDAll.at(I[j]));
		tempBounce.push_back(BounceAll.at(I[j]));
		tempeventID.push_back(eventIDAll.at(I[j]));
		tempSurfIn.push_back(SurfInAll.at(I[j]));
	}

	CellsAll.clear();
	DNflagAll.clear();
	CellTimeEvAll.clear();
	TrackLengthAll.clear();
	ThetaInAll.clear();
	SecondaryIDAll.clear();
	BounceAll.clear();
	eventIDAll.clear();
	SurfInAll.clear();

	int ntemp = tempCellTime.size();
	double tempCellTimeL[ntemp];
	int tempI[ntemp];
	for(int j = 0; j < ntemp; j++){
		tempCellTimeL[j] = tempCellTime.at(j);
	}
	tempCellTime.clear();

	TMath::Sort(ntemp, tempCellTimeL, tempI, false);

	for(int j = 0; j < ntemp; j++){
		sCells = tempCells.at(tempI[j]);
		sDNflag = tempDNflag.at(tempI[j]);
		sCellTime = tempCellTimeL[tempI[j]];
		sCellTimeEv = tempCellTimeEv.at(tempI[j]);
		sTrackLength = tempTrackLength.at(tempI[j]);
		sThetaIn = tempThetaIn.at(tempI[j]);
		sSecondaryID = tempSecondaryID.at(tempI[j]);
		sBounce = tempBounce.at(tempI[j]);
		seventID = tempeventID.at(tempI[j]);
		sSurfIn = tempSurfIn.at(tempI[j]);
		Tnew->Fill();
	}

	tempCells.clear();
	tempDNflag.clear();
	tempCellTimeEv.clear();
	tempTrackLength.clear();
	tempThetaIn.clear();
	tempSecondaryID.clear();
	tempBounce.clear();
	tempeventID.clear();
	tempSurfIn.clear();

	N->cd();
	Tnew->Write("T", TObject::kOverwrite);
	N->Close();
	F->Close();

}		

void preprocessing(TString file){
	TFile* F = TFile::Open(file + ".root");
	TTree* T = (TTree*) F->Get("T");

	double activationTime[81][2668];
	for(int i = 0; i < 81; i++){
		for(int j = 0; j < 2668; j++){
			activationTime[i][j] = -30; // ns
		}
	}

	vector<TTree*> Tnew;
	int Cell = 0, DNflag = 0, Channel = 0, eventID = 0, SurfIn = 0, Bounce = 0, SecondaryID = 0;
	double Time = 0, TimeEv = 0, TrackLength = 0, ThetaIn = 0;

	T->SetBranchAddress("Cells", &Cell);
	T->SetBranchAddress("DNflag", &DNflag);
	T->SetBranchAddress("CellTime", &Time);
	T->SetBranchAddress("CellTimeEv", &TimeEv);
	T->SetBranchAddress("TrackLength", &TrackLength);
	T->SetBranchAddress("ThetaIn", &ThetaIn);
	T->SetBranchAddress("Bounce", &Bounce);
	T->SetBranchAddress("eventID", &eventID);
	T->SetBranchAddress("SurfIn", &SurfIn);

	TFile* N = TFile::Open(file + ".root", "RECREATE");

	for(int i = 0; i < 81; i++){
		Tnew.push_back(new TTree(TString::Format("Ch%02d", i), TString::Format("Ch%02d", i)));
		Tnew.at(i)->Branch("Cell", &Cell);
		Tnew.at(i)->Branch("DNflag", &DNflag);
		Tnew.at(i)->Branch("Time", &Time);
		Tnew.at(i)->Branch("TimeEv", &TimeEv);
		Tnew.at(i)->Branch("TrackLength", &TrackLength);
		Tnew.at(i)->Branch("ThetaIn", &ThetaIn);
		Tnew.at(i)->Branch("SecondaryID", &SecondaryID);
		Tnew.at(i)->Branch("Bounce", &Bounce);
		Tnew.at(i)->Branch("eventID", &eventID);
		Tnew.at(i)->Branch("SurfIn", &SurfIn);
	}


	int Ntot = T->GetEntries();
	for(int i = 0; i < Ntot; i++){
		T->GetEntry(i);
		if(Time - activationTime[Channel][Cell] > 20){
			activationTime[Channel][Cell] = Time;
			Tnew.at(Channel)->Fill();
		}
	}

	F->Close();
	N->cd();
	for(int i = 0; i < 81; i++){
		Tnew.at(i)->Write(TString::Format("Ch%02d", i), TObject::kOverwrite);
	}
	N->Close();
}

double Gp[5] = {20, 5.94028, 39.2627}; // 
double Signal(double x){
	return (x < Gp[0])*exp(-(x-Gp[0])*(x-Gp[0])/2/Gp[1]/Gp[1]) + (x > Gp[0])*exp(-(x-Gp[0])*(x-Gp[0])/2/(Gp[1]*Gp[1] + 0.5*Gp[2]*x));
}

void processing(double threashold, TString file){
	TF1* f = new TF1("f", "(x < [1])*exp(-(x-[1])*(x-[1])/2/([2]*[2] + 2*x*[3]))*[0] + (x>[1])*gaus(0)", 0, 2);
	f->SetNpx(10000);

	ifstream myfile;
	myfile.open("../../../pars.txt");
	double par = 0;
	for(int i = 0; i < 3; i++){
		myfile >> par;
		f->SetParameter(i + 1,  par);
	}
	f->SetParameter(0, 1);
	myfile.close();

	TFile* F = TFile::Open(file + ".root", "UPDATE");
	TTree* T;
	TTree* Twaves = new TTree("waves", "signals");
	
	double deltaT = 400; // ns
	double pitch  = 0.1; // ns

	vector<double> signal;
	vector<double> signalNoSmearing;
	vector<double> signalT;
	signal.resize(int(deltaT / pitch) + 1);
	signalNoSmearing.resize(int(deltaT / pitch) + 1);
	for(int i = 0; i < int(deltaT / pitch) + 1; i++) signalT.push_back(-deltaT + pitch * i);
	vector<double> activationTime;
	vector<double> smearing;

	int Channel, eventID = 0, signalEventID = 0, SurfIn = 0, signalSurfIn = 0, Bounce = 0, signalBounce = 0, SecondaryID = 0, signalSecondaryID = 0;
	double signalTime, signalTimeEv = 0, Amplitude = 0, AmplitudeNoSmearing = 0, Charge = 0, ChargeNoSmearing = 0, TrackLength = 0, signalTrackLength = 0, DeltaTime = 0, ThetaIn = 0, signalThetaIn = 0;

	int signalDN = 0;
	int DNflag = 0;

	TGraph G, GNoSmearing;
	
	Twaves->Branch("Channel", &Channel);
	Twaves->Branch("DN", &signalDN);
	Twaves->Branch("Amplitude", &Amplitude);
	Twaves->Branch("AmplitudeNoSmearing", &AmplitudeNoSmearing);
	Twaves->Branch("Charge", &Charge);
	Twaves->Branch("ChargeNoSmearing", &ChargeNoSmearing);
	Twaves->Branch("Time", &signalTime);
	Twaves->Branch("TimeEv", &signalTimeEv);
	Twaves->Branch("DeltaTime", &DeltaTime);
	Twaves->Branch("TrackLength", &signalTrackLength);
	Twaves->Branch("ThetaIn", &signalThetaIn);
	Twaves->Branch("SecondaryID", &signalSecondaryID);
	Twaves->Branch("Bounce", &signalBounce);
	Twaves->Branch("eventID", &signalEventID);
	Twaves->Branch("SurfIn", &signalSurfIn);

	int Nsignals = 0, thCheck = 0, timeCheck = - 30;

	double globalTime = 0; // ns, time of the Last Point
	double Time = 0, TimeEv = 0;

	double totalTime = 0;

	vector<double> vecTime, vecTimeEv, vecTrackLength, vecThetaIn;

	vector<int> vecDN, vecEventID, vecSurfIn, vecBounce, vecSecondaryID;
	
	// Buffer
	double bufferTimeEv = 0;
	double bufferTrackLength = 0;
	double bufferThetaIn = 0;
	int bufferDN = 0;
	int bufferEventID = 0;
	int bufferSurfIn = 0;
	int bufferBounce = 0;
	int bufferSecondaryID = 0;
	
	int bufferCheck = 0;
	
	for(int i = 0; i < 81; i++){
		Nsignals = 0;
		thCheck = 0;
		Channel = i;
		T = (TTree*) F->Get(TString::Format("Ch%02d", i));
		if(T->GetEntries() > 0){
			T->SetBranchAddress("Time", &Time);
			T->SetBranchAddress("TimeEv", &TimeEv);
			T->SetBranchAddress("DNflag", &DNflag);
			T->SetBranchAddress("TrackLength", &TrackLength);
			T->SetBranchAddress("ThetaIn", &ThetaIn);
			T->SetBranchAddress("SecondaryID", &SecondaryID);
			T->SetBranchAddress("Bounce", &Bounce);
			T->SetBranchAddress("eventID", &eventID);
			T->SetBranchAddress("SurfIn", &SurfIn);
			int last = -1;
			int n_entries = T->GetEntries();
			for(int j = 0; j < T->GetEntries(); j++){
				if(int(100. * j / n_entries)%10 == 0 && int(100. * j / n_entries)/10 > last){
					TString output = TString::Format("Ch. %02d, advancement: [", i);
					last = int(100. * j / n_entries)/10;
					for(int k = 0; k < 10; k++){
						if (k < last) output = output + "=";
						else if (k == last) output = output + ">";
						else output = output + " ";
					}
					output = output + "]";
					cout << output << endl;
				}

				T->GetEntry(j);
				
				int g = -1; // to skip dead time
				bool check = false;
				while(globalTime < Time){
					if(check) g++;
					if(Time - globalTime > 100 && g == -1){
						check = true;
					}
					if(g + 1 == int(deltaT/pitch/2)){
						globalTime = Time - 10;
						check = false;
						g = -1;
					}
					signal.erase(signal.begin());
					signalNoSmearing.erase(signalNoSmearing.begin());
					signalT.erase(signalT.begin());
					double amplitude = 0;
					double amplitudeNoSmearing = 0;
					int n = activationTime.size();
					int actual = 0;
					for(int k = 0; k < n; k++){
						amplitude += smearing.at(actual)*Signal(globalTime - activationTime.at(actual));
						amplitudeNoSmearing += Signal(globalTime - activationTime.at(actual));
						if(globalTime - activationTime.at(actual) > deltaT/2 - 1){ 
							activationTime.erase(activationTime.begin() + actual);
							smearing.erase(smearing.begin() + actual);
							actual --;
						}
						actual ++;
					}
					signal.push_back(amplitude);
					signalNoSmearing.push_back(amplitudeNoSmearing);
					globalTime += pitch;
					signalT.push_back(globalTime);
					
					// Save in buffer and in arrays signal informations
					if(signal.at(int(deltaT/pitch)) > threashold && bufferCheck == 0){
						bufferCheck = 1;
						vecTimeEv.push_back(bufferTimeEv);
						vecDN.push_back(bufferDN);
						vecTrackLength.push_back(bufferTrackLength);
						vecThetaIn.push_back(bufferThetaIn);
						vecSecondaryID.push_back(bufferSecondaryID);
						vecBounce.push_back(bufferBounce);
						vecEventID.push_back(bufferEventID);
						vecSurfIn.push_back(bufferSurfIn);
					}
					else if(signal.at(int(deltaT/pitch)) < threashold && bufferCheck == 1){
						bufferCheck = 0;
					}
					
					// Register signals
					if(signal.at(int(deltaT/pitch/2)) > threashold && thCheck == 0){
						Nsignals += 1;
						thCheck = 1;
						timeCheck = signalT.at(int(deltaT/pitch/2));
						signalTime = globalTime - deltaT/2;
						G = TGraph(signal.size());
						GNoSmearing = TGraph(signalNoSmearing.size());
						for(int l = 0; l < signal.size(); l++){
							G.SetPoint(l, signalT.at(l), signal.at(l));
							GNoSmearing.SetPoint(l, signalT.at(l), signalNoSmearing.at(l));
						}
						
						signalTimeEv = vecTimeEv.at(0);
						signalDN = vecDN.at(0);
						signalTrackLength = vecTrackLength.at(0);
						signalThetaIn = vecThetaIn.at(0);
						signalSecondaryID = vecSecondaryID.at(0);
						signalBounce = vecBounce.at(0);
						signalEventID = vecEventID.at(0);
						signalSurfIn = vecSurfIn.at(0);
						
						//cout << vecTimeEv.size() << endl,
					
						vecTimeEv.erase(vecTimeEv.begin());
						vecDN.erase(vecDN.begin());
						vecTrackLength.erase(vecTrackLength.begin());
						vecThetaIn.erase(vecThetaIn.begin());
						vecSecondaryID.erase(vecSecondaryID.begin());
						vecBounce.erase(vecBounce.begin());
						vecEventID.erase(vecEventID.begin());
						vecSurfIn.erase(vecSurfIn.begin());
					
						int l = 0;
						double DeltaT = 0;
						while(true){
							if(int(deltaT/pitch/2) + l < int(deltaT / pitch) + 1){
								Charge += (signal.at(int(deltaT/pitch/2) + l)) * pitch/2.6743304;
								ChargeNoSmearing += (signalNoSmearing.at(int(deltaT/pitch/2) + l)) * pitch/2.6743304;
							}
							
							if(int(deltaT/pitch/2) + l >= deltaT/pitch + 1){
								DeltaT = signalT.at(int(deltaT/pitch/2) + l - 1) - signalT.at(int(deltaT/pitch/2));
								break;
							}
							else if(signal.at(int(deltaT/pitch/2) + l) < threashold){
								DeltaT = signalT.at(int(deltaT/pitch/2) + l) - signalT.at(int(deltaT/pitch/2));
								break;
							}
							l++;
						}
						Amplitude = TMath::MaxElement(int(DeltaT/pitch), &G.GetY()[int(deltaT/pitch/2)]);
						AmplitudeNoSmearing = TMath::MaxElement(int(DeltaT/pitch), &GNoSmearing.GetY()[int(deltaT/pitch/2)]);
						DeltaTime -= signalTime;
						Twaves->Fill();
						DeltaTime = signalTime;
						Charge = 0;
						ChargeNoSmearing = 0;
					}
					else if(signal.at(int(deltaT/pitch/2)) < threashold && thCheck == 1){
	// && signalT.at(int(deltaT/pitch/2)) - timeCheck > 20
						thCheck = 0;
					}
				}
				activationTime.push_back(Time);
				smearing.push_back(f->GetRandom());
				
				bufferTimeEv = TimeEv;
				bufferDN = DNflag;
				bufferTrackLength = TrackLength;
				bufferThetaIn = ThetaIn;
				bufferBounce = Bounce;
				bufferEventID = eventID;
				bufferSurfIn = SurfIn;
			}

			for(int j = 0; j < int(deltaT/pitch/2); j++){
				signal.erase(signal.begin());
				signalNoSmearing.erase(signalNoSmearing.begin());
				signalT.erase(signalT.begin());
				double amplitude = 0;
				double amplitudeNoSmearing = 0;
				int n = activationTime.size();
				int actual = 0;
				for(int k = 0; k < n; k++){
					amplitude += smearing.at(actual)*Signal(globalTime - activationTime.at(actual));
					amplitudeNoSmearing += Signal(globalTime - activationTime.at(actual));
					if(globalTime - activationTime.at(actual) > deltaT/2 - 1){ 
						activationTime.erase(activationTime.begin() + actual);
						smearing.erase(smearing.begin() + actual);
						actual --;
					}
					actual ++;
				}
				signal.push_back(amplitude);
				signalNoSmearing.push_back(amplitudeNoSmearing);
				globalTime += pitch;
				signalT.push_back(globalTime);
				
				// Save in buffer and in arrays signal informations
				if(signal.at(int(deltaT/pitch)) > threashold && bufferCheck == 0){
					bufferCheck = 1;
					vecTimeEv.push_back(bufferTimeEv);
					vecDN.push_back(bufferDN);
					vecTrackLength.push_back(bufferTrackLength);
					vecThetaIn.push_back(bufferThetaIn);
					vecSecondaryID.push_back(bufferSecondaryID);
					vecBounce.push_back(bufferBounce);
					vecEventID.push_back(bufferEventID);
					vecSurfIn.push_back(bufferSurfIn);
				}
				else if(signal.at(int(deltaT/pitch)) < threashold && bufferCheck == 1){
					bufferCheck = 0;
				}
				
				// Register signals
				if(signal.at(int(deltaT/pitch/2)) > threashold && thCheck == 0){
					Nsignals += 1;
					thCheck = 1;
					timeCheck = signalT.at(int(deltaT/pitch/2));
					signalTime = globalTime - deltaT/2;
					G = TGraph(signal.size());
					GNoSmearing = TGraph(signalNoSmearing.size());
					for(int l = 0; l < signal.size(); l++){
						G.SetPoint(l, signalT.at(l), signal.at(l));
						GNoSmearing.SetPoint(l, signalT.at(l), signalNoSmearing.at(l));
					}
					
					signalTimeEv = vecTimeEv.at(0);
					signalDN = vecDN.at(0);
					signalTrackLength = vecTrackLength.at(0);
					signalThetaIn = vecThetaIn.at(0);
					signalSecondaryID = vecSecondaryID.at(0);
					signalBounce = vecBounce.at(0);
					signalEventID = vecEventID.at(0);
					signalSurfIn = vecSurfIn.at(0);
				
					vecTimeEv.erase(vecTimeEv.begin());
					vecDN.erase(vecDN.begin());
					vecTrackLength.erase(vecTrackLength.begin());
					vecThetaIn.erase(vecThetaIn.begin());
					vecSecondaryID.erase(vecSecondaryID.begin());
					vecBounce.erase(vecBounce.begin());
					vecEventID.erase(vecEventID.begin());
					vecSurfIn.erase(vecSurfIn.begin());
				
					int l = 0;
					double DeltaT = 0;
					while(true){
						if(int(deltaT/pitch/2) + l < int(deltaT / pitch) + 1){
							Charge += (signal.at(int(deltaT/pitch/2) + l)) * pitch/2.6743304;
							ChargeNoSmearing += (signalNoSmearing.at(int(deltaT/pitch/2) + l)) * pitch/2.6743304;
						}
						
						if(int(deltaT/pitch/2) + l >= int(deltaT / pitch) + 1){
							DeltaT = signalT.at(int(deltaT/pitch/2) + l - 1) - signalT.at(int(deltaT/pitch/2));
							break;
						}
						else if(signal.at(int(deltaT/pitch/2) + l) < threashold){
							DeltaT = signalT.at(int(deltaT/pitch/2) + l) - signalT.at(int(deltaT/pitch/2));
							break;
						}
						l++;
					}
					Amplitude = TMath::MaxElement(int(DeltaT/pitch), &G.GetY()[int(deltaT/pitch/2)]);
					AmplitudeNoSmearing = TMath::MaxElement(int(DeltaT/pitch), &GNoSmearing.GetY()[int(deltaT/pitch/2)]);
					DeltaTime -= signalTime;
					Twaves->Fill();
					DeltaTime = signalTime;
					Charge = 0;
					ChargeNoSmearing = 0;
				}
				else if(signal.at(int(deltaT/pitch/2)) < threashold && thCheck == 1){
	// && signalT.at(int(deltaT/pitch/2)) - timeCheck > 20
					thCheck = 0;
				}
			}
			signal.clear();
			signalNoSmearing.clear();
			signal.resize(int(deltaT / pitch) + 1);
			signalNoSmearing.resize(int(deltaT / pitch) + 1);
			if(globalTime > 0) totalTime = globalTime;
			globalTime = 0;
			timeCheck = -30;
			vecTimeEv.clear();
			vecDN.clear();
			vecTrackLength.clear();
			vecThetaIn.clear();
			vecSecondaryID.clear();
			vecBounce.clear();
			vecEventID.clear();
			vecSurfIn.clear();
			DeltaTime = 0;
		}
	}
	Twaves->Write("waves", TObject::kOverwrite);
	F->Close();
}

void reprocess(double threashold, TString name){
	TFile* F = TFile::Open(name + ".root");
	TTree* T = (TTree*) F->Get("waves");

	double Amplitude = 0;
	int Channel = 0;
	int Counts[81] = {0};
	double Time = 0, deltaT = 0, Time0 = 0;

	T->SetBranchAddress("Amplitude", &Amplitude);
	T->SetBranchAddress("Channel", &Channel);
	T->SetBranchAddress("Time", &Time);

	long n_entries = T->GetEntries();
	for(long i = 0; i < n_entries; i++){
		if(i%(n_entries/100) == 0) cout << 100*i/n_entries << " %" << endl;
		T->GetEntry(i);
		if(i == 0){deltaT = Time; Time0 = Time;}
		else if (Time > deltaT) deltaT = Time;
		if(Time < Time0) Time0 = Time;
		if(Amplitude >= threashold) Counts[Channel]++;
	}
	deltaT -= Time0;

	ofstream myfile;
	myfile.open(name + ".txt");
	for(int i = 0; i < 81; i++){
		myfile << Counts[i] << endl;
	}
	myfile << 0.22 << endl;
	myfile << deltaT/1e9 << endl;

	F->Close();
}


void reprocess(double threashold, TString name, long a1, long a2){
	TFile* F = TFile::Open(name + ".root");
	TTree* T = (TTree*) F->Get("waves");

	double Amplitude = 0;
	int Channel = 0;
	int Counts[81] = {0};
	double Time = 0, deltaT = 0, Time0 = 0;

	T->SetBranchAddress("Amplitude", &Amplitude);
	T->SetBranchAddress("Channel", &Channel);
	T->SetBranchAddress("Time", &Time);

	int n_entries = T->GetEntries();
	TRandom3* random = new TRandom3(299792458);
	if(a1 < n_entries && a2 << n_entries){
		for(long i = a1; i < a2; i++){
			if(i%((a2 - a1)/100) == 0) cout << 100*i/(a2 - a1) << " %" << endl;
			T->GetEntry(int(random->Uniform(0, n_entries)));
			if(Amplitude >= threashold) Counts[Channel]++;
		}
	}
	T->GetEntry(0);
	deltaT = Time;
	T->GetEntry(n_entries - 1);
	deltaT = Time - deltaT;
	deltaT = deltaT*(a2 - a1)/n_entries;
	ofstream myfile;
	myfile.open(name + ".txt");
	for(int i = 0; i < 81; i++){
		myfile << Counts[i] << endl;
	}
	myfile << 0.22 << endl;
	myfile << deltaT/1e9 << endl;

	F->Close();
}



void reprocessOneSide(double threashold, TString name){
	TFile* F = TFile::Open(name + ".root");
	TTree* T = (TTree*) F->Get("waves");

	double Amplitude = 0;
	int Channel = 0, SurfIn = 0;
	int Counts[81] = {0};
	double Time = 0, deltaT = 0, Time0 = 0;

	T->SetBranchAddress("Amplitude", &Amplitude);
	T->SetBranchAddress("Channel", &Channel);
	T->SetBranchAddress("Time", &Time);
	T->SetBranchAddress("SurfIn", &SurfIn);

	int n_entries = T->GetEntries();
	for(int i = 0; i < n_entries; i++){
		if(i%(n_entries/100) == 0) cout << 100*i/n_entries << " %" << endl;
		T->GetEntry(i);
		if(i == 0){deltaT = Time; Time0 = Time;}
		else if (Time > deltaT) deltaT = Time;
		if(Time < Time0) Time0 = Time;
		if(Amplitude >= threashold && SurfIn == 4) Counts[Channel]++;
	}
	deltaT -= Time0;

	ofstream myfile;
	myfile.open(name + ".txt");
	for(int i = 0; i < 81; i++){
		myfile << Counts[i] << endl;
	}
	myfile << 0.22 << endl;
	myfile << deltaT/1e9 << endl;

	F->Close();
}

void reprocessCharge(double threashold, TString name){
	TFile* F = TFile::Open(name + ".root");
	TTree* T = (TTree*) F->Get("waves");

	double Charge = 0;
	int Channel = 0;
	int Counts[81] = {0};
	double Time = 0, deltaT = 0, Time0 = 0;

	T->SetBranchAddress("Charge", &Charge);
	T->SetBranchAddress("Channel", &Channel);
	T->SetBranchAddress("Time", &Time);

	int n_entries = T->GetEntries();
	for(int i = 0; i < n_entries; i++){
		if(i%(n_entries/100) == 0) cout << 100*i/n_entries << " %" << endl;
		T->GetEntry(i);
		if(i == 0){deltaT = Time; Time0 = Time;}
		else if (Time > deltaT) deltaT = Time;
		if(Time < Time0) Time0 = Time;
		if(Charge >= threashold) Counts[Channel]++;
	}
	deltaT -= Time0;

	ofstream myfile;
	myfile.open(name + ".txt");
	for(int i = 0; i < 81; i++){
		myfile << Counts[i] << endl;
	}
	myfile << 0.22 << endl;
	myfile << deltaT/1e9 << endl;

	F->Close();
}

void reprocessDN(double threashold, TString name, double thrDN){
	TFile* F = TFile::Open(name + ".root");
	TTree* T = (TTree*) F->Get("waves");

	double Amplitude = 0;
	double DN = 0;
	int Channel = 0;
	int Counts[81] = {0};
	double Time = 0, deltaT = 0, Time0 = 0;

	T->SetBranchAddress("Amplitude", &Amplitude);
	T->SetBranchAddress("DN", &DN);
	T->SetBranchAddress("Channel", &Channel);
	T->SetBranchAddress("Time", &Time);

	int n_entries = T->GetEntries();
	for(int i = 0; i < n_entries; i++){
		if(i%(n_entries/100) == 0) cout << 100*i/n_entries << " %" << endl;
		T->GetEntry(i);
		if(i == 0){deltaT = Time; Time0 = Time;}
		else if (Time > deltaT) deltaT = Time;
		if(Time < Time0) Time0 = Time;
		if(Amplitude >= threashold && DN > thrDN) Counts[Channel]++;
	}

	deltaT -= Time0;

	ofstream myfile;
	myfile.open(name + ".txt");
	for(int i = 0; i < 81; i++){
		myfile << Counts[i] << endl;
	}
	myfile << 0.22 << endl;
	myfile << deltaT/1e9 << endl;

	F->Close();
}


void signals(double threashold){
	ordering("data.root", "out", true);
	preprocessing("out");
	processing(0.1, "out");
	reprocess(threashold, "out");
}

void signals(TString name){
	ordering(name, "out", true);
	preprocessing("out");
	processing(0.1, "out");
	reprocess(10, "out");
}

void signals(bool darkNoise){
	ordering("data.root", "out", darkNoise);
	preprocessing("out");
	processing(0.1, "out");
	reprocess(10, "out");
}

void signals(TString name, TString out){
	ordering(name, out, true);
	preprocessing(out);
	processing(0.1, out);
	reprocess(10, out);
}

void signals(TString name, TString out, bool darkNoise){
	ordering(name, out, darkNoise);
	preprocessing(out);
	processing(0.1, out);
	reprocess(10, out);
}

void signals(TString name, double threashold, TString out){
	ordering(name, out, true);
	preprocessing(out);
	processing(threashold, out);
	reprocess(threashold, out);
}

void signals(TString name, double threashold, TString out, bool darkNoise){
	ordering(name, out, darkNoise);
	preprocessing(out);
	processing(threashold, out);
	reprocess(threashold, out);
}


void signals(){
/*
	ordering("data.root", "out", true);
	preprocessing("out");
	processing(0.1, "out");
	reprocess(10, "out");
//
	signals("full/mergelooseMu037MeV1e7AlMask.root", "out/mergelooseMu037MeV1e7AlMaskout");
	signals("full/mergelooseMu037MeV1e7NoMask.root", "out/mergelooseMu037MeV1e7NoMaskout");
	signals("full/mergelooseMu037MeV1e7.root", "out/mergelooseMu037MeV1e7out"); 
	signals("full/mergeloosePos528MeV1e7AlMask.root", "out/mergeloosePos528MeV1e7AlMaskout");
	signals("full/mergeloosePos528MeV1e7NoMaskBounce.root", "out/mergeloosePos528MeV1e7NoMaskBounceout");
 	signals("full/mergeloosePos528MeV1e7NoMask.root", "out/mergeloosePos528MeV1e7NoMaskout");
	signals("full/mergeloosePos528MeV1e7.root", "out/mergeloosePos528MeV1e7out");

//*/

	processing(0.1, "out/mergelooseMu037MeV1e7NoMaskout");
	processing(0.1, "out/mergelooseMu037MeV1e7out");
	processing(0.1, "out/mergelooseMu037MeV1e7AlMaskout");
	processing(0.1, "out/mergeloosePos528MeV1e7AlMaskout");
	processing(0.1, "out/mergeloosePos528MeV1e7NoMaskBounceout");
	processing(0.1, "out/mergeloosePos528MeV1e7NoMaskout");
	processing(0.1, "out/mergeloosePos528MeV1e7out");
}

