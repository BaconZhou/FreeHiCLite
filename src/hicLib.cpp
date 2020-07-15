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
  Rcpp::colnames(ans) = Rcpp::CharacterVector({"x", "y", "counts"});
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
  std::sort(ans.begin(), ans.end());
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

Rcpp::List hicDataInformation(FreeHiC::Juicer::hicReader *reader) {
  reader->init();
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
  Rcpp::List info = Rcpp::List::create(
    Rcpp::Named("genomeID") = reader->getGenome(),
    Rcpp::Named("resolution") = resolutions,
    Rcpp::Named("pairs") = reader->getPairs(),
    Rcpp::Named("chromosomeSizes") = Rcpp::DataFrame::create(
      Rcpp::Named("chromosome") = chroName,
      Rcpp::Named("size") = chroSize));
  return info;
}

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
  
  Rcpp::List ans = Rcpp::List::create(Rcpp::Named("contact") = contact,
                                      Rcpp::Named("information") = hicDataInformation(reader));
  return ans;
}

// [[Rcpp::export]]
Rcpp::List hicDataFragSites(const std::string &fileName, bool isHttp,
                            const std::vector<std::string> &chromosomes) {
  
  FreeHiC::Juicer::hicReader *reader;
  if (isHttp) {
    reader = new FreeHiC::Juicer::hicReaderHttp(fileName, "FRAG", 1);
  } else {
    reader = new FreeHiC::Juicer::hicReader(fileName, "FRAG", 1);
  }
  
  Rcpp::List fragmentSites;
  for (const std::string &chr : chromosomes) {
    fragmentSites.push_back(reader->getChromosomeSites(chr), chr);
  }
  return fragmentSites;
}


// [[Rcpp::export]]
Rcpp::List hicDataInformation(const std::string &fileName, bool isHttp) {
  FreeHiC::Juicer::hicReader *reader;
  if (isHttp) {
    reader = new FreeHiC::Juicer::hicReaderHttp(fileName, "FRAG", 1);
  } else {
    reader = new FreeHiC::Juicer::hicReader(fileName, "FRAG", 1);
  }
  return hicDataInformation(reader);
}

// [[Rcpp::export]]
Rcpp::List hicDataSimuFromFile(const std::string &fileName,  bool isHttp,
                               const std::vector<std::string> &pair,  
                               const std::string &unit, int resolution,
                               const int &sequenceDepth, const double &countScale, 
                               const double &noiseRate, const double &neighborZeroRate) {
  FreeHiC::Juicer::hicReader *reader;
  if (isHttp) {
    reader = new FreeHiC::Juicer::hicReaderHttp(fileName, unit, resolution);
  } else {
    reader = new FreeHiC::Juicer::hicReader(fileName, unit, resolution);
  }
  
  FreeHiC::FreeHiC simulator(resolution, sequenceDepth, countScale,
                             noiseRate, neighborZeroRate);

  std::unordered_map<std::string, std::vector<FreeHiC::contactRecord>> resMap;
  
  for (size_t i = 0; i < pair.size(); i++) {
    std::pair<std::string, std::string> chrPair = split(pair[i], "_");
    std::vector<FreeHiC::contactRecord> res = reader->readData(chrPair.first, chrPair.second); 
    resMap[pair[i]] = res;
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
Rcpp::List hicDataSimuList(const Rcpp::List &contactRecords, 
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
Rcpp::IntegerMatrix hicDataSimuMatrix(const Rcpp::IntegerMatrix &contactRecords, 
                       int resolution, const int &sequenceDepth, const double &countScale, 
                       const double &noiseRate, const double &neighborZeroRate) {
  
  FreeHiC::FreeHiC simulator(resolution, sequenceDepth, countScale,
                             noiseRate, neighborZeroRate);
  
  std::unordered_map<std::string, std::vector<FreeHiC::contactRecord>> resMap;
  resMap["1_1"] = matrixToContactVector(contactRecords);
  simulator.simulate(resMap);
  // cout << "check random number: " << simulator.rng.uniform() << endl;
  std::unordered_map<std::string, std::vector<FreeHiC::contactRecord>> res = simulator.getData();
  return contactVectorToMatrix(res["1_1"]);
}


// [[Rcpp::export]]
Rcpp::IntegerMatrix spikein(const Rcpp::IntegerMatrix &background, 
                            const Rcpp::IntegerMatrix &SpikeInSignal, 
                            int bandwidth, bool smooth) {
  std::vector<FreeHiC::contactRecord> backgound_ = matrixToContactVector(background);
  std::vector<FreeHiC::contactRecord> SpikeInSignal_ = matrixToContactVector(SpikeInSignal);
  
  FreeHiC::Spikein sol(backgound_, bandwidth, smooth);
  
  sol.smoothSignal(SpikeInSignal_);
  return contactVectorToMatrix(sol.getData());
}
  
  
  
  
  