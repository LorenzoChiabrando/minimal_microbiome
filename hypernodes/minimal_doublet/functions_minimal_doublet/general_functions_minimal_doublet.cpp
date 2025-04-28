
// Import libraries and header files
#include <string.h>
#include <sys/resource.h>
#include <map>
#include <regex>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <string>
#include <sstream>
#include <cmath>
#include <iomanip>    // Include for formatting output

// Bring the FBGLPK namespace into the current scope
using namespace FBGLPK;

double V = 0.001; // (mL) (1 microL)
long long int delta = 1e+10; // (cell/mL)
long long int max_total_bacteria = V*delta;
double total_bacteria = 0; //total number of bacteria actually in the simulation 

// CarbonFlux transition
double SharedMetabolite(double *Value,
                  vector<class FBGLPK::LPprob>& vec_fluxb,
                  map <string,int>& NumTrans,
                  map <string,int>& NumPlaces,
                  const vector<string> & NameTrans,
                  const struct InfTr* Trans,
                  const int T,
                  const double& time
                  ) {
  
  double rate = 0;
  
  return(rate);
  
}

vector<string> bacteria_names = {
    "E_coli", "Clost_buty"
};


vector<string> bacteriaBiomass_names = {
    "E_coli_biomass_e", "Clost_buty_biomass_e"
};


void updateTotalBacteria(double *Value, 
                         const map<string, 
                         int>& NumPlaces, 
                         const vector<class FBGLPK::LPprob>& vec_fluxb) {
    total_bacteria = 0;
    for (size_t i = 0; i < bacteria_names.size(); ++i) {
        auto placeIt = NumPlaces.find(bacteria_names[i]);
        auto biomassIt = NumPlaces.find(bacteriaBiomass_names[i]);
        if (placeIt != NumPlaces.end() && biomassIt != NumPlaces.end()) {
            total_bacteria += floor(Value[placeIt->second]);
        }
    }
}

// Debug rate
static std::map<std::string, std::ofstream> outfileMap;
void WriteRateToFile(const std::string& bacteriumName, double time, double rate) {

    if (outfileMap.find(bacteriumName) == outfileMap.end()) {
        std::string fileName = "../logStarvationRate" + bacteriumName + ".csv";
        outfileMap[bacteriumName].open(fileName, std::ios::out);
        outfileMap[bacteriumName] << "Time,Rate" << std::endl;
    }

    outfileMap[bacteriumName] << std::fixed << std::setprecision(16) << time << "," << rate << std::endl;
}

// Map to store file streams for each bacterium
static std::map<std::string, std::ofstream> deathOutfileMap;
static std::map<std::string, std::ofstream> duplicationOutfileMap;

void WriteDeathRateToFile(const std::string& bacteriumName, double time, double rate) {
    if (deathOutfileMap.find(bacteriumName) == deathOutfileMap.end()) {
        std::string fileName = "../logDeathRate" + bacteriumName + ".csv";
        deathOutfileMap[bacteriumName].open(fileName, std::ios::out);
        deathOutfileMap[bacteriumName] << "Time,Rate" << std::endl;
    }
    deathOutfileMap[bacteriumName] << std::fixed << std::setprecision(16) << time << "," << rate << std::endl;
}

void WriteDuplicationRateToFile(const std::string& bacteriumName, double time, double rate) {
    if (duplicationOutfileMap.find(bacteriumName) == duplicationOutfileMap.end()) {
        std::string fileName = "../logDuplicationRate" + bacteriumName + ".csv";
        duplicationOutfileMap[bacteriumName].open(fileName, std::ios::out);
        duplicationOutfileMap[bacteriumName] << "Time,Rate" << std::endl;
    }
    duplicationOutfileMap[bacteriumName] << std::fixed << std::setprecision(16) << time << "," << rate << std::endl;
}

double Starvation(double *Value,
                  vector<class FBGLPK::LPprob>& vec_fluxb,
                  map <string,int>& NumTrans,
                  map <string,int>& NumPlaces,
                  const vector<string> & NameTrans,
                  const struct InfTr* Trans,
                  const int T,
                  const double& time,
                  const double starvation_rate,
                  long unsigned int bacterium) {

    if (bacterium < 0 || bacterium >= bacteria_names.size()) {
        std::cerr << "Invalid bacterium index: " << bacterium << std::endl;
        return 0; 
    }

    const std::string& bacteriumName = bacteria_names[bacterium];
    
    auto placeIndex = NumPlaces.find(bacteriumName);
    auto biomassIndex = NumPlaces.find(bacteriaBiomass_names[bacterium]);
    if (placeIndex == NumPlaces.end() || biomassIndex == NumPlaces.end()) {
        std::cerr << "Place index not found for bacterium index: " << bacterium << std::endl;
        return 0; 
    }

    double currentBiomass = Value[biomassIndex->second];
    double BioMin = vec_fluxb[bacterium].getBioMin();
    double num_places_specie = floor(Value[placeIndex->second]);
    double rate = 0;
		double effectiveBiomass = (currentBiomass > BioMin) ? (currentBiomass - BioMin) : 0.0;
		double starvRate = starvation_rate * effectiveBiomass;

    if (num_places_specie >= 1) {
        rate = starvRate; // starvation_rate * (currentBiomass * 1.0); 
    }

    
    cout << "Starvation rate for: " << bacteriumName << " is --->  " << rate << " (pgDW/h) "<< "time ---> " << time << endl;
    // Uncomment the following line when debugging is needed
    // WriteRateToFile(bacteriumName, time, rate);

    return rate;
}



double Death(double *Value,
             vector<class FBGLPK::LPprob>& vec_fluxb,
             map <string,int>& NumTrans,
             map <string,int>& NumPlaces,
             const vector<string> & NameTrans,
             const struct InfTr* Trans,
             const int T,
             const double& time,
             const double half_life,
             long unsigned int bacterium) {
             
  
  if (bacterium < 0 || bacterium >= bacteria_names.size()) {
      cerr << "Invalid bacterium index in Death: " << bacterium << endl;
      return 0;
  }

  double rate = 0;
  auto placeIt = NumPlaces.find(bacteria_names[bacterium]);
  auto biomassIt = NumPlaces.find(bacteriaBiomass_names[bacterium]);
  if (placeIt != NumPlaces.end() && biomassIt != NumPlaces.end()) {
    double currentBiomass = Value[biomassIt->second];
    double num_places_specie = floor(Value[placeIt->second]);
    double meanBiomass = vec_fluxb[bacterium].getBioMean();
    
    if(currentBiomass > 0 && num_places_specie >= 1) {
        rate = half_life * num_places_specie * (meanBiomass / currentBiomass);
    }
    
  } else {
     cerr << "Place not found for bacterium index: " << bacterium << endl;
  }

  const std::string& bacteriumName = bacteria_names[bacterium];
  cout << "Death rate for: " << bacteriumName << " is --->  " << rate << " (cell/h) "<< "time ---> " << time << endl;
  // Uncomment the following line when debugging is needed
  // WriteDeathRateToFile(bacteria_names[bacterium], time, rate);
  return rate;
}


double Duplication(double *Value,
						       vector<class FBGLPK::LPprob>& vec_fluxb,
						       map <string,int>& NumTrans,
						       map <string,int>& NumPlaces,
						       const vector<string> & NameTrans,
						       const struct InfTr* Trans,
						       const int T,
						       const double& time,
						       const double duplication_rate,
						       long unsigned int bacterium) {
						       
  if (bacterium < 0 || bacterium >= bacteria_names.size()) {
      cerr << "Invalid bacterium index in Duplication: " << bacterium << endl;
      return 0;
  }

  double rate = 0;
  auto placeIt = NumPlaces.find(bacteria_names[bacterium]);
  auto biomassIt = NumPlaces.find(bacteriaBiomass_names[bacterium]);
  if (placeIt != NumPlaces.end() && biomassIt != NumPlaces.end()) {
    double currentBiomass = Value[biomassIt->second];
    double num_places_specie = floor(Value[placeIt->second]);
    double maxBiomass = vec_fluxb[bacterium].getBioMax();
    //double minBiomass = vec_fluxb[bacterium].getBioMin();
    
    updateTotalBacteria(Value, NumPlaces, vec_fluxb);

    if (num_places_specie >= 1 && total_bacteria < max_total_bacteria /*&& currentBiomass >= minBiomass*/) {
        rate = num_places_specie * duplication_rate * (currentBiomass / maxBiomass) * (1 - (total_bacteria / static_cast<double>(max_total_bacteria)));
    }
  } else {
      cerr << "Place not found for bacterium index: " << bacterium << endl;
  }
  
  const std::string& bacteriumName = bacteria_names[bacterium];
  cout << "Duplication rate for: " << bacteriumName << " is --->  " << rate << " (cell/h) "<< "time ---> " << time << endl;
  // Uncomment the following line when debugging is needed
  // WriteDuplicationRateToFile(bacteria_names[bacterium], time, rate);
  
  return rate;
}
