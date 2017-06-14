/*
 * JEventProcessorPScsTACCalibration.h
 *
 *  Created on: May 22, 2017
 *      Author: Hovanes Egiyan
 */

#ifndef JEVENTPROCESSOR_PSVSTACCALIBRATION_H_
#define JEVENTPROCESSOR_PSVSTACCALIBRATION_H_

#include <iostream>
#include <map>
#include <vector>
#include <iterator>
#include <algorithm>

#include <TH1.h>
#include <TDirectory.h>

#include <JANA/JEventProcessor.h>
#include <DANA/DApplication.h>
#include <TRIGGER/DL1Trigger.h>
#include <TAGGER/DTAGHHit.h>
#include <TAGGER/DTAGMHit.h>
#include <RF/DRFTime.h>
#include <PAIR_SPECTROMETER/DPSCHit.h>

class JEventProcessor_PSvsTAC_Calibration: public jana::JEventProcessor {
protected:
	// Map of all histograms for this monitoring plugin. the first index identifies the
	// name of the histogram, the second index (inner) identifies the trigger bit.
	std::map<std::string, std::map<unsigned, TH1*> > histoMap;

	// ROOT file name
	std::string rootFileName = "tac_monitor.root";

	// ROOT directory pointer
	TDirectory* rootDir = nullptr;

	// Mask indicating which trigger bits this class cares for.
	static uint32_t tacTriggerMask;
	// Mask indicating which trigger bits this class cares for.
	static uint32_t psTriggerMask;
	// Maximum numbr of trigger bits considered in this plugin
	static uint32_t numberOfTriggerBits;

	// Threshold that will define when the TAC hit occurred
	static unsigned tacThreshold;

	static std::string tacRebuildFunctor;

	// Timing cut value between the TAGH and TAC coincidence
	static double timeCutValue_TAGH;
	// Timing cut width between the TAGH and TAC coincidence
	static double timeCutWidth_TAGH;

	// Timing cut value between the TAGH and TAC coincidence
	static double timeCutValue_TAGM;
	// Timing cut width between the TAGH and TAC coincidence
	static double timeCutWidth_TAGM;

	virtual jerror_t init(void);          ///< Called once at program start.
	virtual jerror_t brun(jana::JEventLoop *eventLoop, int32_t runNumber);          ///< Called everytime a new run number is detected.
	virtual jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventNumber);          ///< Called every event.
	virtual jerror_t erun(void);          ///< Called every time run number changes, provided brun has been called.
	virtual jerror_t fini(void);          ///< Called after last event of last event source has been processed.

	// Fill TAC-related histograms
	virtual jerror_t fillHistosTAC(jana::JEventLoop* eventLoop,
			uint32_t trigBit);
	// Fill PS-related histograms
	virtual jerror_t fillHistosPS(jana::JEventLoop* eventLoop,
			uint32_t trigBit);

	// Method where the histograms are created
	virtual jerror_t createHistograms();
	virtual jerror_t createHistogramsForTAC(unsigned trigBit);
	virtual jerror_t createHistogramsForPS(unsigned trigBit);

	template<typename TH1_TYPE>
	jerror_t createHisto(unsigned trigBit, std::string key,
			std::string titlePrefix, std::string xTitle, int nBins, double xMin,
			double xMax);
	template<typename TH2_TYPE>
	jerror_t createHisto(unsigned trigBit, std::string histKey,
			std::string xTitlePrefix, std::string xTitle, std::string yTitle,
			int nBinsX, double xMin, double xMax, int nBinsY, double yMin,
			double yMax);

	// Write histograms into the file
	virtual jerror_t writeHistograms();

	// Check if the trigger bits for the event are useful
	static bool triggerIsUseful(const DL1Trigger* trigWords) {
		if (trigWords) {
			return triggerIsUseful(trigWords->trig_mask);
		}
		return false;
	}
	// Check if the trigger bits for the event are useful for PS
	static bool triggerIsUsefulForPS(const DL1Trigger* trigWords) {
		if (trigWords) {
			return triggerIsUsefulForPS(trigWords->trig_mask);
		}
		return false;
	}
	// Check if the trigger bits for the event are useful for TAC
	static bool triggerIsUsefulForTAC(const DL1Trigger* trigWords) {
		if (trigWords) {
			return triggerIsUsefulForTAC(trigWords->trig_mask);
		}
		return false;
	}

	// Check if the trigger bits for the event are useful
	static bool triggerIsUseful(unsigned trigBits) {
		if ((trigBits & (tacTriggerMask | psTriggerMask)) == 0) {
			return false;
		}
		return true;
	}
	// Check if the trigger bits for the event are useful for PS
	static bool triggerIsUsefulForPS(unsigned trigBits) {
		if ((trigBits & psTriggerMask) == 0) {
			return false;
		}
		return true;
	}
	// Check if the trigger bits for the event are useful for TAC
	static bool triggerIsUsefulForTAC(unsigned trigBits) {
		if ((trigBits & tacTriggerMask) == 0) {
			return false;
		}
		return true;
	}

public:
	JEventProcessor_PSvsTAC_Calibration() {
		// TODO Auto-generated constructor stub

	}
	virtual ~JEventProcessor_PSvsTAC_Calibration() {
		// TODO Auto-generated destructor stub
	}

	TDirectory* getRootDir() {
		return rootDir;
	}

	void setRootDir(TDirectory* rootDir) {
		this->rootDir = rootDir;
	}

	const std::map<std::string, std::map<unsigned, TH1*> >& getHistoMap() const {
		return histoMap;
	}

	void setHistoMap(
			const std::map<std::string, std::map<unsigned, TH1*> >& histoMap) {
		this->histoMap = histoMap;
	}

	const std::string& getRootFileName() const {
		return rootFileName;
	}

	void setRootFileName(const std::string& rootFileName =
			"ps_vs_tac_calib.root") {
		this->rootFileName = rootFileName;
	}

	static double getTimeCutValueTAGH() {
		return timeCutValue_TAGH;
	}

	void setTimeCutValueTAGH(double timeCutValueTagh) {
		JEventProcessor_PSvsTAC_Calibration::timeCutValue_TAGH =
				timeCutValueTagh;
	}

	static double getTimeCutValueTAGM() {
		return timeCutValue_TAGM;
	}

	void setTimeCutValueTAGM(double timeCutValueTagm) {
		JEventProcessor_PSvsTAC_Calibration::timeCutValue_TAGM =
				timeCutValueTagm;
	}

	static double getTimeCutWidthTAGH() {
		return timeCutWidth_TAGH;
	}

	void setTimeCutWidthTAGH(double timeCutWidthTagh) {
		JEventProcessor_PSvsTAC_Calibration::timeCutWidth_TAGH =
				timeCutWidthTagh;
	}

	static double getTimeCutWidthTAGM() {
		return timeCutWidth_TAGM;
	}

	void setTimeCutWidthTAGM(double timeCutWidthTagm) {
		JEventProcessor_PSvsTAC_Calibration::timeCutWidth_TAGM =
				timeCutWidthTagm;
	}
};


// Create a 1D histogram of type TH1_TYPE and assign it to the histogram map based on the argument valeus
// provided in the function call.
template<typename TH1_TYPE>
jerror_t JEventProcessor_PSvsTAC_Calibration::createHisto(unsigned trigBit,
		string histKey, string titlePrefix, string xTitle, int nBins,
		double xMin, double xMax) {
	static_assert(std::is_base_of<TH1, TH1_TYPE>::value,
	              "TH1_TYPE must be derived from TH1");
	stringstream histName;
	stringstream histTitle;
	histName << histKey << "_" << trigBit;
	histTitle << titlePrefix << trigBit;
	if (histoMap.count(histName.str()) == 0)
		histoMap[histName.str()] = map<unsigned, TH1*>();
	histoMap[histKey][trigBit] = new TH1_TYPE(histName.str().c_str(),
			histTitle.str().c_str(), nBins, xMin, xMax);
	histoMap[histKey][trigBit]->GetXaxis()->SetTitle(xTitle.c_str());
	return NOERROR;
}

// Create a 2D histogram of type TH2_TYPE and assign it to the histogram map based on the argument valeus
// provided in the function call.
template<typename TH2_TYPE>
jerror_t JEventProcessor_PSvsTAC_Calibration::createHisto(unsigned trigBit,
		string histKey, string titlePrefix, string xTitle, string yTitle,
		int nBinsX, double xMin, double xMax, int nBinsY, double yMin,
		double yMax) {
	static_assert(std::is_base_of<TH2, TH2_TYPE>::value,
	              "TH2_TYPE must be derived from TH2");
	stringstream histName;
	stringstream histTitle;
	histName << histKey << "_" << trigBit;
	histTitle << titlePrefix << trigBit;
	if (histoMap.count(histName.str()) == 0)
		histoMap[histName.str()] = map<unsigned, TH1*>();
	histoMap[histKey][trigBit] = new TH2_TYPE(histName.str().c_str(),
			histTitle.str().c_str(), nBinsX, xMin, xMax, nBinsY, yMin, yMax);
	histoMap[histKey][trigBit]->GetXaxis()->SetTitle(xTitle.c_str());
	histoMap[histKey][trigBit]->GetYaxis()->SetTitle(yTitle.c_str());
	return NOERROR;
}



#endif /* JEVENTPROCESSOR_PSVSTACCALIBRATION_H_ */
