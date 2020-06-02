#include "hicReader.hpp"
#include "FreeHiC.hpp"
#include "Spikein.hpp"

// [[Rcpp::plugins(cpp11)]]

/* Utilitis */

Rcpp::IntegerMatrix contactVectorToMatrix(const std::vector<FreeHiC::contactRecord> &records) {
  
  size_t N = records.size();
  Rcpp::IntegerMatrix ans(N, 3);
  for (size_t i = 0; i < N; i++) {
    FreeHiC::contactRecord record = records[i];
    if (record.counts < 0.5) continue;
    ans(i, 0) = record.binX;
    ans(i, 1) = record.binY;
    ans(i, 2) = record.counts;
  }
  return ans;
}

std::vector<FreeHiC::contactRecord> matrixToContactVector(const Rcpp::IntegerMatrix &records) {
  
  size_t N = records.nrow();
  std::vector<FreeHiC::contactRecord> ans(N);
  for (size_t i = 0; i < N; i++) {
    int binX = records(i, 0);
    int binY = records(i, 1);
    float count = records(i, 2);
    ans[i] = FreeHiC::contactRecord(binX, binY, count);
  }
  return ans;
}

std::unordered_map<std::string, std::vector<FreeHiC::contactRecord>> listToMap(const Rcpp::List &records, 
                                                                               const std::vector<std::string> &names) {
  std::unordered_map<std::string, std::vector<FreeHiC::contactRecord>> ans;
  int N = names.size();
  for (int i = 0; i < N; ++ i) {
    std::string name = names[i];
    ans[name] = matrixToContactVector(records[name]);
  }
  return ans;
}

std::pair<std::string, std::string> split(std::string s, std::string delimiter) {
  size_t pos = 0;
  std::string t1, t2;
  std::pair<std::string, std::string> ans;
  pos = s.find(delimiter);
  t1 = s.substr(0, pos);
  s.erase(0, pos + delimiter.length());
  pos = s.find(delimiter);
  t2 = s.substr(0, pos);
  s.erase(0, pos + delimiter.length());
  ans = std::make_pair(t1, t2);
  return ans;
}


// READ data

// [[Rcpp::export]]
Rcpp::IntegerMatrix hicData(std::string fileName, std::string chr1, 
                            std::string chr2, std::string unit, 
                            int resolution) {
  FreeHiC::Juicer::hicReader reader(fileName,  unit, resolution);
  
  const std::vector<FreeHiC::contactRecord> records = reader.readData(chr1, chr2);
  
  size_t N = records.size();
  Rcpp::IntegerMatrix ans(N, 3);
  for (size_t i = 0; i < N; i++) {
    FreeHiC::contactRecord record = records[i];
    ans(i, 0) = record.binX;
    ans(i, 1) = record.binY;
    ans(i, 2) = record.counts;
  }
  return ans;
}
  

// [[Rcpp::export]]
Rcpp::IntegerMatrix hicDataHttp(std::string fileName, std::string chr1, 
                            std::string chr2, std::string unit, 
                            int resolution) {
  FreeHiC::Juicer::hicReaderHttp reader(fileName, unit, resolution);
  
  std::vector<FreeHiC::contactRecord> records = reader.readData(chr1, chr2);
  
  size_t N = records.size();
  Rcpp::IntegerMatrix ans(N, 3);
  for (size_t i = 0; i < N; i++) {
    FreeHiC::contactRecord record = records[i];
    ans(i, 0) = record.binX;
    ans(i, 1) = record.binY;
    ans(i, 2) = record.counts;
  }
  return ans;
}

// [[Rcpp::export]]
Rcpp::List hicDataExtra(const std::string &fileName, bool isHttp,
                        const std::vector<std::string> &pair, 
                        const std::string &unit, int resolution) {
  
  FreeHiC::Juicer::hicReader *reader;
  if (isHttp) {
    reader = new FreeHiC::Juicer::hicReaderHttp(fileName, unit, resolution);
  } else {
    reader = new FreeHiC::Juicer::hicReader(fileName, unit, resolution);
  }
  
  
  Rcpp::List contact;
  
  for (size_t i = 0; i < pair.size(); i++) {
      std::pair<std::string, std::string> chrPair = split(pair[i], "_");
      std::vector<FreeHiC::contactRecord> res = reader->readData(chrPair.first, chrPair.second); 
      contact.push_back(contactVectorToMatrix(res), pair[i]);
  }

  Rcpp::List resolutions;
  for (const auto & it : reader->getResolutions()) {
    resolutions.push_back(it.second, it.first);
  }
  
  std::vector<std::string> chroName;
  std::vector<int> chroSize;
  Rcpp::DataFrame chromosomeSize;
  for (const auto & it: reader->getChromosomeSize()) {
    chroName.push_back(it.first);
    chroSize.push_back(it.second);
  }
  
  Rcpp::List ans = Rcpp::List::create(Rcpp::Named("contact") = contact,
                                      Rcpp::Named("information") = Rcpp::List::create(
                                        Rcpp::Named("genomId") = reader->getGenome(),
                                        Rcpp::Named("resolution") = resolutions,
                                        Rcpp::Named("chromosomeSizes") = Rcpp::DataFrame::create(
                                          Rcpp::Named("chromosome") = chroName,
                                          Rcpp::Named("size") = chroSize)
                                      ));
  return ans;
}


// [[Rcpp::export]]
Rcpp::List hicDataSimu(const std::string &fileName, 
                       const std::vector<std::string> &chromosomes, 
                       const std::string &unit, int resolution,
                       const int &sequenceDepth, const double &countScale, 
                       const double &noiseRate, const double &neighborZeroRate) {
  FreeHiC::Juicer::hicReader reader(fileName, unit, resolution);
  FreeHiC::FreeHiC simulator(resolution, sequenceDepth, countScale,
                             noiseRate, neighborZeroRate);

  std::unordered_map<std::string, std::vector<FreeHiC::contactRecord>> resMap;
  for (int i = 0; i < chromosomes.size(); i++) {
    for (int j = i; j < chromosomes.size(); j++) {
      std::stringstream ss;
      ss << chromosomes[i] << "_" << chromosomes[j];
      std::vector<FreeHiC::contactRecord> res = reader.readData(chromosomes[i], chromosomes[j]);
    
      resMap[ss.str()] = res;
    }
  }
  simulator.simulate(resMap);
  std::unordered_map<std::string, std::vector<FreeHiC::contactRecord>> res = simulator.getData();
  
  Rcpp::List ans;
  for (const auto &it : res) {
    ans.push_back(contactVectorToMatrix(it.second), it.first);
  }
  return ans;
}

// [[Rcpp::export]]
Rcpp::List hicDataSimuPure(const Rcpp::List &contactRecords, 
                           const std::vector<std::string> &names,
                           int resolution, const int &sequenceDepth, const double &countScale, 
                           const double &noiseRate, const double &neighborZeroRate) {
  
  FreeHiC::FreeHiC simulator(resolution, sequenceDepth, countScale,
                             noiseRate, neighborZeroRate);
  
  std::unordered_map<std::string, std::vector<FreeHiC::contactRecord>> resMap = listToMap(contactRecords, names);
  simulator.simulate(resMap);
  std::unordered_map<std::string, std::vector<FreeHiC::contactRecord>> res = simulator.getData();
  
  Rcpp::List ans;
  for (const auto &it : res) {
    ans.push_back(contactVectorToMatrix(it.second), it.first);
  }
  return ans;
}


// [[Rcpp::export]]
Rcpp::IntegerMatrix spikein(const Rcpp::IntegerMatrix &background, 
                            const Rcpp::IntegerMatrix &SpikeInSignal, 
                            int resolution, bool smooth) {
  std::vector<FreeHiC::contactRecord> backgound_ = matrixToContactVector(background);
  std::vector<FreeHiC::contactRecord> SpikeInSignal_ = matrixToContactVector(SpikeInSignal);
  
  FreeHiC::Spikein sol(backgound_, resolution, smooth);
  
  sol.smoothSignal(SpikeInSignal_);
  return contactVectorToMatrix(sol.getData());
}
  
  
  
  
  