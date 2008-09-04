#include "RecoParticleFlow/PFClusterTools/interface/PFClusterCalibration.h"
#include "RecoParticleFlow/PFClusterTools/interface/CalibrationProvenance.h"

#include <cmath>
#include <cassert>
#include <TBranch.h>

using namespace pftools;

PFClusterCalibration::PFClusterCalibration() :
	barrelEndcapEtaDiv_(1.0), ecalOnlyDiv_(0.3), hcalOnlyDiv_(0.5),
			doCorrection_(1), globalP0_(0.0), globalP1_(1.0), lowEP0_(0.0),
			lowEP1_(1.0), allowNegativeEnergy_(0), correction_("correction",
					"((x-[0])/[1])*(x>[4])+((x-[2])/[3])*(x<[4])") {
	std::cout << __PRETTY_FUNCTION__ << std::endl;
	correction_.FixParameter(0, globalP0_);
	correction_.FixParameter(1, globalP1_);
	correction_.FixParameter(2, lowEP0_);
	correction_.FixParameter(3, lowEP1_);
	correction_.FixParameter(4, correctionLowLimit_);

	std::string eoeb("ecalOnlyEcalBarrel");
	names_.push_back(eoeb);
	std::string eoee("ecalOnlyEcalEndcap");
	names_.push_back(eoee);
	std::string hohb("hcalOnlyHcalBarrel");
	names_.push_back(hohb);
	std::string hohe("hcalOnlyHcalEndcap");
	names_.push_back(hohe);

	std::string eheb("ecalHcalEcalBarrel");
	names_.push_back(eheb);
	std::string ehee("ecalHcalEcalEndcap");
	names_.push_back(ehee);
	std::string ehhb("ecalHcalHcalBarrel");
	names_.push_back(ehhb);
	std::string ehhe("ecalHcalHcalEndcap");
	names_.push_back(ehhe);

	char
			* funcString("([0]*[5]*x*([1]-[5]*x)/pow(([2]+[5]*x),3)+[3]*pow([5]*x, 0.1))*([5]*x<[8] && [5]*x>[7])+[4]*([5]*x>[8])+([6]*[5]*x)*([5]*x<[7])");

	//Create functions for each sector
	for (std::vector<std::string>::const_iterator cit = names_.begin(); cit
			!= names_.end(); ++cit) {
		std::string name = *cit;
		TF1 func(name.c_str(), funcString);
		//some sensible defaults
		func.FixParameter(5, 1);
		func.FixParameter(6, 1);
		func.FixParameter(7, 0);
		func.SetMinimum(0);
		//Store in map
		namesAndFunctions_[name] = func;

	}
	std::cout
			<< "WARNING! Evolution functions have not yet been initialised - ensure this is done.\n";
	std::cout << "PFClusterCalibration construction complete."<< std::endl;

}

void PFClusterCalibration::setEvolutionParameters(const std::string& sector,
		std::vector<double> params) {
	TF1* func = &(namesAndFunctions_.find(sector)->second);
	unsigned count(0);
	std::cout << "Fixing for "<< sector << "\n";
	for (std::vector<double>::const_iterator dit = params.begin(); dit
			!= params.end(); ++dit) {
		func->FixParameter(count, *dit);
		std::cout << "\t"<< count << ": "<< *dit;
		++count;
	}
	std::cout << std::endl;
	func->SetMinimum(0);
}

void PFClusterCalibration::setCorrections(const double& lowEP0,
		const double& lowEP1, const double& globalP0, const double& globalP1) {
	globalP0_ = globalP0;
	globalP1_ = globalP1;
	lowEP0_ = lowEP0;
	lowEP1_ = lowEP1;
	correctionLowLimit_ = (lowEP0_ - globalP0_)/(globalP1_ - lowEP1_);

	correction_.FixParameter(0, globalP0_);
	correction_.FixParameter(1, globalP1_);
	correction_.FixParameter(2, lowEP0_);
	correction_.FixParameter(3, lowEP1_);
	correction_.FixParameter(4, correctionLowLimit_);

	std::cout << __PRETTY_FUNCTION__ << ": setting correctionLowLimit_ = "
			<< correctionLowLimit_ << "\n";
}

double PFClusterCalibration::getCalibratedEcalEnergy(const double& ecalE,
		const double& hcalE, const double& eta, const double& phi) const {
	const TF1* theFunction(0);
	if (ecalE > ecalOnlyDiv_ && hcalE > hcalOnlyDiv_) {
		//ecalHcal class
		if (fabs(eta) < barrelEndcapEtaDiv_) {
			//barrel
			theFunction = &(namesAndFunctions_.find("ecalHcalEcalBarrel")->second);
		} else {
			//endcap
			theFunction = &(namesAndFunctions_.find("ecalHcalEcalEndcap")->second);
		}
	} else if (ecalE > ecalOnlyDiv_ && hcalE < hcalOnlyDiv_) {
		//ecalOnly class
		if (fabs(eta) < barrelEndcapEtaDiv_)
			theFunction = &(namesAndFunctions_.find("ecalOnlyEcalBarrel")->second);
		else
			theFunction = &(namesAndFunctions_.find("ecalOnlyEcalEndcap")->second);
	} else {
		//either hcal only or too litte energy, in any case,
		return ecalE;
	}
	assert(theFunction != 0);
	double totalE(ecalE + hcalE);
	double bCoeff = theFunction->Eval(totalE);
	return ecalE * bCoeff;
}

double PFClusterCalibration::getCalibratedHcalEnergy(const double& ecalE,
		const double& hcalE, const double& eta, const double& phi) const {
	const TF1* theFunction(0);
	if (ecalE > ecalOnlyDiv_ && hcalE > hcalOnlyDiv_) {
		//ecalHcal class
		if (fabs(eta) < barrelEndcapEtaDiv_) {
			//barrel
			theFunction = &(namesAndFunctions_.find("ecalHcalHcalBarrel")->second);
		} else {
			//endcap
			theFunction = &(namesAndFunctions_.find("ecalHcalHcalEndcap")->second);
		}
	} else if (ecalE < ecalOnlyDiv_ && hcalE> hcalOnlyDiv_) {
		//hcalOnly class
		if (fabs(eta) < barrelEndcapEtaDiv_)
		theFunction = &(namesAndFunctions_.find("hcalOnlyHcalBarrel")->second);
		else
		theFunction = &(namesAndFunctions_.find("hcalOnlyHcalEndcap")->second);
	} else {
		//either ecal only or too litte energy, in any case,
		return hcalE;
	}
	double totalE(ecalE + hcalE);
	assert(theFunction != 0);
	double cCoeff = theFunction->Eval(totalE);
	return hcalE * cCoeff;
}

double PFClusterCalibration::getCalibratedEnergy(const double& ecalE,
		const double& hcalE, const double& eta, const double& phi) const {
	double totalE(ecalE + hcalE);
	double answer(totalE);

	answer = getCalibratedEcalEnergy(ecalE, hcalE, eta, phi)
			+ getCalibratedHcalEnergy(ecalE, hcalE, eta, phi);

	//apply correction?
	if (!doCorrection_)
		return answer;
	if (!allowNegativeEnergy_ && correction_.Eval(answer) < 0)
		return 0;
	return correction_.Eval(answer);

}

const void PFClusterCalibration::calibrate(Calibratable& c) {
	CalibrationResultWrapper crw;
	getCalibrationResultWrapper(c, crw);
	c.calibrations_.push_back(crw);

}

const void PFClusterCalibration::getCalibrationResultWrapper(
		const Calibratable& c, CalibrationResultWrapper& crw) {

	crw.ecalEnergy_ = getCalibratedEcalEnergy(c.cluster_energyEcal_,
			c.cluster_energyHcal_, fabs(c.cluster_meanEcal_.eta_),
			fabs(c.cluster_meanEcal_.phi_));

	crw.hcalEnergy_ = getCalibratedHcalEnergy(c.cluster_energyEcal_,
			c.cluster_energyHcal_, fabs(c.cluster_meanEcal_.eta_),
			fabs(c.cluster_meanEcal_.phi_));

	crw.particleEnergy_ = getCalibratedEnergy(c.cluster_energyEcal_,
			c.cluster_energyHcal_, fabs(c.cluster_meanEcal_.eta_),
			fabs(c.cluster_meanEcal_.phi_));

	crw.provenance_ = LINEARCORR;
	crw.b_ = crw.ecalEnergy_ / c.cluster_energyEcal_;
	crw.c_ = crw.hcalEnergy_ / c.cluster_energyHcal_;
	crw.truthEnergy_ = c.sim_energyEvent_;

}

void PFClusterCalibration::calibrateTree(TTree* input) {
	std::cout << __PRETTY_FUNCTION__
			<< ": WARNING! This isn't working properly yet!\n";
	TBranch* calibBr = input->GetBranch("Calibratable");
	Calibratable* calib_ptr = new Calibratable();
	calibBr->SetAddress(&calib_ptr);

//	TBranch* newBranch = input->Branch("NewCalibratable",
//			"pftools::Calibratable", &calib_ptr, 32000, 2);

	std::cout << "Looping over tree's "<< input->GetEntries()<< " entries...\n";
	for (unsigned entries(0); entries < 20000; entries++) {
		if (entries % 10000== 0)
			std::cout << "\tProcessing entry "<< entries << "\n";
		input->GetEntry(entries);
		calibrate(*calib_ptr);
		input->Fill();
	}
	//input.Write("",TObject::kOverwrite);
}

std::ostream& pftools::operator<<(std::ostream& s, const PFClusterCalibration& cc) {
	s << "PFClusterCalibration: dump...\n";
	s << "barrelEndcapEtaDiv:\t" << cc.barrelEndcapEtaDiv_ << ", ecalOnlyDiv:\t" << cc.ecalOnlyDiv_;
	s << ", \nhcalOnlyDiv:\t" << cc.hcalOnlyDiv_ << ", doCorrection:\t" << cc.doCorrection_;
	s << ", \nallowNegativeEnergy:\t" << cc.allowNegativeEnergy_;
	s << ", \ncorrectionLowLimit:\t" << cc.correctionLowLimit_ << ",parameters:\t" << cc.globalP0_ << ", ";
	s << cc.globalP1_ << ", " << cc.lowEP0_ << ", " << cc.lowEP1_;
	return s;
}

