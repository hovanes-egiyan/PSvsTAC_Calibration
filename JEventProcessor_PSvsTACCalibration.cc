/*
 * JEventProcessorPScsTACCalibration.cc
 *
 *  Created on: May 22, 2017
 *      Author: Hovanes Egiyan
 */

#include <iostream>
#include <map>
#include <vector>
#include <sstream>

#include "TApplication.h"  // needed to display canvas
#include "TSystem.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TF1.h"
#include "TF2.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"

#include <JANA/JApplication.h>
#include <DANA/ReadWriteLock.h>
#include <TRIGGER/DL1Trigger.h>
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250PulseData.h>
#include <DAQ/Df250WindowRawData.h>
#include <DAQ/DF1TDCHit.h>
#include <DAQ/DEPICSvalue.h>
#include <TAC/DTACDigiHit.h>
#include <TAC/DTACTDCDigiHit.h>
#include <TAC/DTACHit.h>
#include <TAGGER/DTAGHDigiHit.h>
#include <TAGGER/DTAGHTDCDigiHit.h>
#include <TAGGER/DTAGMDigiHit.h>
#include <TAGGER/DTAGMTDCDigiHit.h>
#include <TTAB/DTTabUtilities.h>

#include "JEventProcessor_PSvsTACCalibration.h"

using namespace jana;
using namespace std;

// Routine used to create our JEventProcessor
extern "C" {
void InitPlugin(JApplication *app) {
	cout << "Here I go" << endl;
	InitJANAPlugin(app);
	cout << "Plugin initialized" << endl;
	app->AddProcessor(new JEventProcessor_PSvsTAC_Calibration());
	cout << "Processor added" << endl;
}
}

// Mask that specifies the bits of interest for TAC runs
uint32_t JEventProcessor_PSvsTAC_Calibration::tacTriggerMask = 0b00000010;
// Mask that specifies the bits of interest for TAC runs
uint32_t JEventProcessor_PSvsTAC_Calibration::psTriggerMask = 0b00000001;
// Maximum numbr of trigger bits considered in this plugin
uint32_t JEventProcessor_PSvsTAC_Calibration::numberOfTriggerBits = 16;

// Threshold that will define when the TAC hit occurred
unsigned JEventProcessor_PSvsTAC_Calibration::tacThreshold = 500;

string JEventProcessor_PSvsTAC_Calibration::tacRebuildFunctor = "";

// Timing cut value between the TAGH and TAC coincidence in ns
double JEventProcessor_PSvsTAC_Calibration::timeCutValue_TAGH = 100.0;
// Timing cut width between the TAGH and TAC coincidence in ns
double JEventProcessor_PSvsTAC_Calibration::timeCutWidth_TAGH = 20;

// Timing cut value between the TAGM and TAC coincidence in ns
double JEventProcessor_PSvsTAC_Calibration::timeCutValue_TAGM = 90.0;
// Timing cut width between the TAGM and TAC coincidence in ns
double JEventProcessor_PSvsTAC_Calibration::timeCutWidth_TAGM = 20;

jerror_t JEventProcessor_PSvsTAC_Calibration::init(void) {
	cout << "Executing JEventProcessor_PSvsTAC_Calibration::init()" << endl;
	volatile WriteLock rootRWLock(
			*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());

	cout << "lock is taken" << endl;
	// Create parameters and assign values
	gPARMS->SetDefaultParameter<string, double>("TAC:TAGH_FADC_MEAN_TIME",
			timeCutValue_TAGH);
	gPARMS->GetParameter("TAC:TAGH_FADC_MEAN_TIME")->GetValue(
			timeCutValue_TAGH);
	gPARMS->SetDefaultParameter<string, double>("TAC:TAGM_FADC_MEAN_TIME",
			timeCutValue_TAGM);
	gPARMS->GetParameter("TAC:TAGM_FADC_MEAN_TIME")->GetValue(
			timeCutValue_TAGM);

	gPARMS->SetDefaultParameter<string,unsigned>( "TAC:THRESHOLD", tacThreshold );
	gPARMS->GetParameter( "TAC:THRESHOLD" )->GetValue( tacThreshold );

	gPARMS->SetDefaultParameter<string,string>( "TAC:REBUILD_FUNC", tacRebuildFunctor );
	gPARMS->GetParameter( "TAC:REBUILD_FUNC" )->GetValue( tacRebuildFunctor );

	cout << "Parameters are created " << endl;

	// Create TAC directory and the histograms
	TDirectory *mainDir = gDirectory;
	rootDir = gDirectory->mkdir("TAC");
	rootDir->cd();
	createHistograms();
	mainDir->cd();
	cout << "Done executing JEventProcessor_PSvsTAC_Calibration::init()"
			<< endl;
	return NOERROR;
}

jerror_t JEventProcessor_PSvsTAC_Calibration::brun(jana::JEventLoop* eventLoop,
		int32_t runNumber) {
	stringstream fileNameStream;
	fileNameStream << "ps_vs_tac_calib_" << runNumber << ".root";
	rootFileName = fileNameStream.str();
	return NOERROR;
}

jerror_t JEventProcessor_PSvsTAC_Calibration::evnt(jana::JEventLoop* eventLoop,
		uint64_t eventNumber) {
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// loop->Get(mydataclasses);
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();

	// do not do anything if the application is not DApplication
	if (dynamic_cast<DApplication*>(japp) == nullptr)
		return NOERROR;

	// Get First Trigger Type
	const DL1Trigger *trigWords = nullptr;
	try {
		eventLoop->GetSingle(trigWords);
	} catch (...) {
	};

	// Decide if to continue considering this event based on the trigger bit pattern
	if (!triggerIsUseful(trigWords))
		return NOERROR;

	// Check the trigger bits and fill histograms
	for (unsigned trigBit = 0; trigBit < numberOfTriggerBits; trigBit++) {
		unsigned singleBit = 1 << trigBit;
		// This is a TAC trigger, fill TAC-trigger-related histograms
		if (tacTriggerMask & singleBit) {
			fillHistosTAC(eventLoop, trigBit);
		}
		// This is a PS trigger, fill PS-trigger-related histograms
		if (psTriggerMask & singleBit) {
			fillHistosPS(eventLoop, trigBit);
		}
	}

	// Write histograms into ROOT file once in a while
	if (eventNumber % 200000 == 0) {
		this->writeHistograms();
	}

	return NOERROR;
}

jerror_t JEventProcessor_PSvsTAC_Calibration::erun(void) {
	this->writeHistograms();
	return NOERROR;
}

jerror_t JEventProcessor_PSvsTAC_Calibration::fini(void) {
	return NOERROR;
}

jerror_t JEventProcessor_PSvsTAC_Calibration::createHistograms() {
	cout << "Creating TAC histos" << endl;
	for (unsigned trigBit = 0; trigBit < numberOfTriggerBits; trigBit++) {
		unsigned trigPattern = 1 << trigBit;
		if (triggerIsUseful(trigPattern)) {
		}

		if (triggerIsUsefulForTAC(trigPattern)) {
			createHistogramsForTAC(trigBit);
		}
		if (triggerIsUsefulForPS(trigPattern)) {
			createHistogramsForPS(trigBit);
		}
	}
	return NOERROR;
}

jerror_t JEventProcessor_PSvsTAC_Calibration::createHistogramsForTAC(
		unsigned trigBit) {
	// Create TAC number of hits histogram
	createHisto<TH1D>(trigBit, "TAC_NHITS", "Number of hits in TAC",
			"number of hits from FADC FPGA [#]", 7, 0., 7.);

	// Create TAC hit from ADC
	createHisto<TH1D>(trigBit, "TAC_TIME", "TAC time",
			"TAC time [ns]", 1200, -30., 30.);

	// Create TAC hit time wrt RF
	createHisto<TH1D>(trigBit, "TAC_RF_TIME", "TAC time wrt RF",
			"TAC - RF time [ns]", 1200, -30., 30.);

	// Create TAC amplitude histos
	createHisto<TH1D>(trigBit, "TAC_ADCAMP",
			"TAC Signal Amplitude for Trigger ", "TAC Amplitude", 500, 0.,
			5000.);
	// Create TAC signal time histo
	createHisto<TH1D>(trigBit, "TAC_ADCTIME", "TAC Signal time for Trigger ",
			"FlashADC peak time (ns)", 400, 0., 400.);

	// Create TAC hit time vs TAC amplitude
	createHisto<TH2D>(trigBit, "TAC_TIME_VS_E", "TAC time vs TAC energy",
			 "TAC energy", "TAC time [ns]", 5000, 0., 10000. , 500, -20., 20. );

	// Create TAC hit time wrt RF vs TAC amplitude
	createHisto<TH2D>(trigBit, "TAC_RF_TIME_VS_E", "TAC time wrt RF vs TAC energy",
			 "TAC energy", "TAC - RF time [ns]", 5000, 0., 10000. , 500, -20., 20. );

	// Create TAGH Hits detector ID
	createHisto<TH1D>(trigBit, "TAC_TAGH_ENERGY",
			"TAC TAGH Hits Energy for Trigger ",
			"Tagger Hodoscope Energy (GeV)", 500, 3., 12.0);
	// Create TAGH Hits detector ID
	createHisto<TH1D>(trigBit, "TAC_TAGH_ENERGY_MATCHED",
			"Matched to TAC TAGH Hits Energy for Trigger ",
			"Tagger Hodoscope Energy (GeV)", 500, 3., 12.0);
	// Create TAGH Hits detector ID
	createHisto<TH1D>(trigBit, "TAC_TAGH_ENERGY_UNMATCHED",
			"Accidental to TAC TAGH Hits Energy for Trigger ",
			"Tagger Hodoscope Energy (GeV)", 500, 3., 12.0);

	// Create TAGH signal time histo
	createHisto<TH1D>(trigBit, "TAC_TAGH_TIME",
			"TAGH Signal time relative to TAC for Trigger ", "TAGH time (ns)",
			600, -300, 300.);
	createHisto<TH1D>(trigBit, "TAC_TAGH_TIME_MATCHED",
			"TAGH Signal time relative to TAC for matched hits for Trigger ",
			"TAGH time (ns)", 600, -300, 300.);
	createHisto<TH1D>(trigBit, "TAC_TAGH_TIME_UNMATCHED",
			"TAGH Signal time relative to TAC for unmatched hits for Trigger ",
			"TAGH time (ns)", 600, -300, 300.);

	// Create TAC amplitude vs TAGH ID histo
	createHisto<TH2D>(trigBit, "TACAMPvsTAGHID",
			"TAC FADC Amplitude vs TAGH ID ",
			"Tagger Hodoscope Det. Number [#]", "FlashADC peak for TAC", 320,
			0., 320., 1000, 10., 5000.);
	// Create TAGH time vs TAGH ID histo
	createHisto<TH2D>(trigBit, "TAC_TAGHTIMEvsTAGHID",
			"TAGH Time reltive to TAC vs TAGH ID ",
			"Tagger Hodoscope Det. Number [#]", "TAGH time", 320, 0., 320., 400,
			0., 400.);

	// Create TAGM Hits detector ID
	createHisto<TH1D>(trigBit, "TAC_TAGM_ID_UNMATCHED",
			"Accidental to TAC TAGM Hits Detector ID for Trigger ",
			"Tagger Microscope Det. Number [#]", 110, 0., 110.);
	// Create TAGM Hits detector ID
	createHisto<TH1D>(trigBit, "TAC_TAGM_ID_MATCHED",
			"Matched to TAC TAGM Hits Detector ID for Trigger ",
			"Tagger Microscope Det. Number [#]", 110, 0., 110.);
	// Create TAGM signal time histo
	createHisto<TH1D>(trigBit, "TAC_TAGM_TIME",
			"TAGM Signal time reltive to TAC for Trigger ",
			"FlashADC peak time (ns)", 400, 0., 400.);
	// Create TAC amplitude vs TAGM ID histo
	createHisto<TH2D>(trigBit, "TACAMPvsTAGMID",
			"TAC FADC Amplitude vs TAGM ID ",
			"Tagger Microscope Det. Number [#]", "FlashADC peak for TAC", 110,
			0., 110., 1000, 10., 5000.);
	// Create TAGH time vs TAGH ID histo
	createHisto<TH2D>(trigBit, "TAC_TAGMTIMEvsTAGMID",
			"TAGM Time relative to TAC vs TAGM ID ",
			"Tagger Microscope Det. Number [#]", "TAGM time", 110, 0., 110.,
			400, 0., 400.);

	return NOERROR;
}

jerror_t JEventProcessor_PSvsTAC_Calibration::createHistogramsForPS(
		unsigned trigBit) {

	// Create PSC TDC hit time vs RF
	createHisto<TH1D>(trigBit, "PSC_TIME", "PSC time",
			"PSC time [ns]", 10000, -300., 300.);

	// Create PSC TDC hit time vs RF
	createHisto<TH1D>(trigBit, "PSC_RF_TIME", "PSC time wrt RF",
			"PSC - RF time [ns]", 10000, -300., 300.);

	// Create TAGH signal time histo
	createHisto<TH1D>(trigBit, "PSC_TAGH_TIME",
			"TAGH Signal time relative to PSC for Trigger ", "TAGH time (ns)",
			600, -300, 300.);
	createHisto<TH1D>(trigBit, "PSC_TAGH_TIME_MATCHED",
			"TAGH Signal time relative to PSC for matched hits for Trigger ",
			"TAGH time (ns)", 10000, -300, 300.);
	createHisto<TH1D>(trigBit, "PSC_TAGH_TIME_UNMATCHED",
			"TAGH Signal time relative to PSC for unmatched hits for Trigger ",
			"TAGH time (ns)", 10000, -300, 300.);


	// Create TAGH Hits detector ID
	createHisto<TH1D>(trigBit, "PSC_TAGH_ENERGY",
			"PSC TAGH Hits Energy for Trigger ",
			"Tagger Hodoscope Energy (GeV)", 500, 3., 12.0);
	// Create TAGH Hits detector ID
	createHisto<TH1D>(trigBit, "PSC_TAGH_ENERGY_MATCHED",
			"Matched to PSC TAGH Hits Energy for Trigger ",
			"Tagger Hodoscope Energy (GeV)", 500, 3., 12.0);
	// Create TAGH Hits detector ID
	createHisto<TH1D>(trigBit, "PSC_TAGH_ENERGY_UNMATCHED",
			"Accidental to PSC TAGH Hits Energy for Trigger ",
			"Tagger Hodoscope Energy (GeV)", 500, 3., 12.0);


	// Create TAGH Hits detector ID
	createHisto<TH1D>(trigBit, "PSC_TAGH_ID_UNMATCHED",
			"Accidental to PS TAGH Hits Detector ID for Trigger ",
			"Tagger Hodoscope Det. Number [#]", 320, 0., 320.);
	// Create TAGH Hits detector ID
	createHisto<TH1D>(trigBit, "PSC_TAGH_ID_MATCHED",
			"Matched to PS TAGH Hits Detector ID for Trigger ",
			"Tagger Hodoscope Det. Number [#]", 320, 0., 320.);
	// Create TAGH signal time histo
	createHisto<TH1D>(trigBit, "PSC_TAGH_TIME",
			"TAGH Signal time relative to PS for Trigger ",
			"FlashADC peak time (ns)", 400, 0., 400.);
	return NOERROR;
}

jerror_t JEventProcessor_PSvsTAC_Calibration::fillHistosTAC(
		jana::JEventLoop* eventLoop, uint32_t trigBit) {

//	vector<const DTACHit*> tacHitOriginalVector;
//	eventLoop->Get(tacHitOriginalVector);
	vector<const DTACHit*> tacRebuildHitVector;
	eventLoop->Get( tacRebuildHitVector, tacRebuildFunctor.c_str() );
//	eventLoop->Get( tacRebuildHitVector, "REBUILD" );
//	eventLoop->Get( tacRebuildHitVector, "REBUILD_SPIKE" );
//	eventLoop->Get( tacRebuildHitVector, "REBUILD_ERFC" );
//	if( tacRebuildHitVector.size() > 0 )
//		cout << "Rebuild time is " << tacRebuildHitVector[0]->getT() << endl;
	// Here we chose if we want to use the original or rebuild TAC hits.
	vector<const DTACHit*>& tacHitVector = tacRebuildHitVector;
//	cout << "Found " << tacHitVector.size() << " TACHit objects" << endl;
	{
		volatile WriteLock rootRWLock(
				*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
		histoMap["TAC_NHITS"][trigBit]->Fill((double) tacHitVector.size());

	}
	if( tacHitVector.size() < 1 ) return NOERROR;
	for (auto& tacHit : tacHitVector) {

		// Make sure that there is a valid TAC hit and the energy is above some
		// reasonable threshold
		if (tacHit == nullptr || tacHit->getE() < tacThreshold )
			continue;
//		const DRFTime* rfTimeObjectBest;
//		eventLoop->GetSingle(rfTimeObjectBest, "", true);

		const DRFTime* rfTimeObjectTOF;
		eventLoop->GetSingle(rfTimeObjectTOF, "TOF", true);
//		const DRFTime* rfTimeObjectPSC;
//		eventLoop->GetSingle(rfTimeObjectPSC, "PSC", true);
//		const DRFTime* rfTimeObjectTAGH;
//		eventLoop->GetSingle(rfTimeObjectTAGH, "TAGH", true);


//		cout << "TAC time in " << tacHit << " is " << tacHit->getT()  << endl;
//		vector<pair<string,string>> outputStrings;
//
//		tacHit->toStrings( outputStrings );
//		cout << endl << "from hit at " << tacHit  << endl;
//		for( auto& outputPair :  outputStrings ) {
//			cout << outputPair.first << "  = " << outputPair.second << endl;
//		}
//
//		cout << "TOF RF is " << rfTimeObjectTOF->dTime << " , PSC RF is "
//				<< rfTimeObjectPSC->dTime << " , TAGH RF time is "
//				<< rfTimeObjectTAGH->dTime << endl;

		if (rfTimeObjectTOF != nullptr) {
			volatile WriteLock rootRWLock(
					*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
			histoMap["TAC_TIME"][trigBit]->Fill(tacHit->getT());
			histoMap["TAC_RF_TIME"][trigBit]->Fill(
					tacHit->getT() - rfTimeObjectTOF->dTime);
			histoMap["TAC_TIME_VS_E"][trigBit]->Fill(tacHit->getE(),
					tacHit->getT());
			histoMap["TAC_RF_TIME_VS_E"][trigBit]->Fill(tacHit->getE(),
					tacHit->getT() - rfTimeObjectTOF->dTime);
		}

		// Declare comparison functor for this TAC hits and the TAGH hits.
		auto compareFunctorTAGH =
				[&tacHit](const DTAGHHit* lhs, const DTAGHHit* rhs ) ->
				bool {return( (lhs!=nullptr)&(rhs!=nullptr ) ? fabs(lhs->t-tacHit->getT()) < fabs(rhs->t-tacHit->getT() ) : false );};

		vector<const DTAGHHit*> taghHitVector;
		eventLoop->Get(taghHitVector);
		for (auto& taghHit : taghHitVector) {
			if (taghHit != nullptr) {
				volatile WriteLock rootRWLock(
						*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
				histoMap["TAC_TAGH_TIME"][trigBit]->Fill(
						taghHit->t - tacHit->getT());
				histoMap["TAC_TAGH_ENERGY"][trigBit]->Fill(taghHit->E);
			}
		}
		if (taghHitVector.size() > 0) {
			// Find best match in TACH and make an entry by looping through the TAGH hits
			std::nth_element(taghHitVector.begin(), taghHitVector.end(),
					taghHitVector.end(), compareFunctorTAGH);
			const DTAGHHit* taghBestHit = taghHitVector[0];
			const DTAGHHit* taghWorstMatch = taghHitVector[taghHitVector.size()
					- 1];
			if (taghBestHit != nullptr) {
				volatile WriteLock rootRWLock(
						*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
				double photonEnergy = taghBestHit->E;
				double photonTime = taghBestHit->t;
				histoMap["TAC_TAGH_ENERGY_MATCHED"][trigBit]->Fill(
						photonEnergy);
				histoMap["TAC_TAGH_TIME_MATCHED"][trigBit]->Fill(
						photonTime - tacHit->getT());
			}
			if (taghWorstMatch != nullptr) {
				volatile WriteLock rootRWLock(
						*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
				double photonEnergy = taghWorstMatch->E;
				double photonTime = taghWorstMatch->t;
				histoMap["TAC_TAGH_TIME_UNMATCHED"][trigBit]->Fill(
						photonTime - tacHit->getT());
				histoMap["TAC_TAGH_ENERGY_UNMATCHED"][trigBit]->Fill(
						photonEnergy);
			}
		}
	}

	return NOERROR;
}

jerror_t JEventProcessor_PSvsTAC_Calibration::fillHistosPS(
		jana::JEventLoop* eventLoop, uint32_t trigBit) {
	vector<const DPSCHit*> pscHitVector;
	eventLoop->Get(pscHitVector);

	for (auto& pscHit : pscHitVector) {

		if (pscHit == nullptr)
			continue;
		if( pscHit->has_TDC && pscHit->arm == DPSGeometry::kNorth ) {
			volatile WriteLock rootRWLock(
					*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
			histoMap["PSC_TIME"][trigBit]->Fill( pscHit->t );
		}
		const DRFTime* rfTimeObjectBest;
//		eventLoop->GetSingle(rfTimeObjectBest, "", true);
//
//		const DRFTime* rfTimeObjectTOF;
//		eventLoop->GetSingle(rfTimeObjectTOF, "TOF", true);
		const DRFTime* rfTimeObjectPSC;
		eventLoop->GetSingle(rfTimeObjectPSC, "PSC", true);
//		const DRFTime* rfTimeObjectTAGH;
//		eventLoop->GetSingle(rfTimeObjectTAGH, "TAGH", true);

		if (rfTimeObjectPSC != nullptr && pscHit->has_TDC && pscHit->arm == DPSGeometry::kNorth ) {
			volatile WriteLock rootRWLock(
					*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
			histoMap["PSC_RF_TIME"][trigBit]->Fill(
					pscHit->t - rfTimeObjectPSC->dTime);
		}

		vector<const DTAGHHit*> taghHitVector;
		eventLoop->Get(taghHitVector);
		for (auto& taghHit : taghHitVector) {
			if (taghHit != nullptr) {
				volatile WriteLock rootRWLock(
						*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());
				histoMap["PSC_TAGH_TIME"][trigBit]->Fill(
						taghHit->t - pscHit->t);
				histoMap["PSC_TAGH_ENERGY"][trigBit]->Fill(taghHit->E);
			}
		}

	}
	return NOERROR;
}

jerror_t JEventProcessor_PSvsTAC_Calibration::writeHistograms() {
	volatile WriteLock rootRWLock(
			*dynamic_cast<DApplication*>(japp)->GetRootReadWriteLock());

	TDirectory* oldDir = gDirectory;
	TFile outFile(rootFileName.c_str(), "RECREATE");
	outFile.cd();
	for (auto& histNameIter : histoMap) {
//		auto& histName =  histNameIter.first;
		auto& histMapTrig = histNameIter.second;
		for (auto& histTrigIter : histMapTrig) {
			auto histPointer = histTrigIter.second;
			histPointer->Write();
		}
	}
	outFile.Write();
	outFile.Close();
	oldDir->cd();

	return NOERROR;
}
