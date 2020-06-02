//
// Created by Peigen Zhou on 5/23/20.
//

#include "FreeHiC.hpp"

namespace FreeHiC {

double FreeContact::getBinomial(double n, double p) {
    if (p < 0 || p > 1)
        return 0;
    double log_q = std::log(1.0 - p);
    double x = 0;
    double sum = 0;
    for (;;) {
        sum += std::log(this->rng->uniform()) / (n - x);
        if (sum < log_q) {
            return x;
        }
        x = x + 1.0;
    }
}

void FreeContact::probIncrement() {
    double r = rng->uniform();
    if (r < 0.8) {
        this->sampleCounts += 1.0;
    } else if (r < 0.95) {
        this->sampleCounts += 2.0;
    } else {
        this->sampleCounts += 3.0;
    }
}

void FreeContact::probDecrement() {
    double r = rng->uniform();
    if (r < 0.8) {
        this->sampleCounts -= 1.0;
    } else if (r < 0.95) {
        this->sampleCounts -= 2.0;
    } else {
        this->sampleCounts -= 3.0;
    }
    this->sampleCounts = std::max(this->sampleCounts, 0.0);
}

void FreeContact::updateNeighborZero(const int &prevX, const int &prevY,
                                     const double &prevCounts, const int &maxX,
                                     const int &maxY, const int &resolution,
                                     const double &totalScaledCounts) {
    int yLength = this->binY - prevY;
    int xLength = (this->binX + prevX) / 2;
    this->area = (double)xLength * yLength / (resolution * resolution);
    this->neighborAverageCounts =
        (this->counts + prevCounts) / (this->area + 1.0);
    this->neighborZeroProb = 2.0 / (this->area + 1.0);
    int counts = this->sampleCounts;

    while (counts > 0) {
#ifdef RVERSION
        Rcpp::checkUserInterrupt();
#endif
        counts -= 1.0;
        if (rng->uniform() > neighborZeroProb)
            continue;
        this->sampleCounts -= 1.0;

        contactRecord tmp;

        int x = (rng->uniform() * maxX) / resolution;
        int y = (rng->uniform() * maxY) / resolution;
        if (maxX == maxY) {
            if (x < y) {
                tmp.binX = x * resolution;
                tmp.binY = y * resolution;
            } else {
                tmp.binX = y * resolution;
                tmp.binY = x * resolution;
            }
        } else {
            tmp.binX = x;
            tmp.binY = y;
        }
        tmp.counts = 1.0;
        this->neighborRecords.emplace_back(tmp);
    }
}

void FreeContact::sample(const double &totalCounts,
                         const double &totalScaleCounts) {
    double prob = this->counts / totalCounts;
    this->sampleCounts = getBinomial(totalScaleCounts, prob);
}

void FreeContact::assignRng(RandomNumber *rng_) { this->rng = rng_; }

void FreeHiC::simulate(
    const std::unordered_map<std::string, std::vector<contactRecord>>
        &dataMap) {
    this->getBasicInformation(dataMap);

    double totalScaledCounts =
        std::max(this->totalCounts * this->countScale_, this->sequenceDepth_);
    for (auto &it : recordsMap) {
        std::string key = it.first;
        std::vector<FreeContact> &data = it.second;
        const std::pair<int, int> &maxBin = pairMaxBin[key];
        // int maxBin = pairMaxBin[key];
        int prevX = 0, prevY = 0;
        double prevCounts = 0;

        for (FreeContact &contact : data) {
#ifdef RVERSION
            Rcpp::checkUserInterrupt();
#endif
            contact.assignRng(&rng);
            if (rng.uniform() < this->noiseRate_)
                contact.noise();
            contact.sample(this->totalCounts, totalScaledCounts);

            if (rng.uniform() < this->neighborZeroRate_) {
                contact.updateNeighborZero(
                    prevX, prevY, prevCounts, maxBin.first, maxBin.second,
                    this->resolution_, totalScaledCounts);
            }

            prevX = contact.getBinX();
            prevY = contact.getBinY();
            prevCounts = contact.getCounts();
        }
    }
}

void FreeHiC::getBasicInformation(
    const std::unordered_map<std::string, std::vector<contactRecord>>
        &dataMap) {

    this->totalCounts = 0.0;
    this->recordsMap.clear();
    this->pairMaxBin.clear();
    int counter = 0;
    for (const auto &it : dataMap) {
        std::string key = it.first;
        const std::vector<contactRecord> &data = it.second;
        int maxBinX = 0;
        int maxBinY = 0;
        std::vector<FreeContact> records;
        records.reserve(data.size());
        for (const contactRecord &contact : data) {
            records.emplace_back(contact);
            totalCounts += contact.counts;
            maxBinX = std::max(maxBinX, contact.binX);
            maxBinY = std::max(maxBinY, contact.binY);
        }
        this->recordsMap[key] = records;
        this->pairMaxBin[key] = std::make_pair(maxBinX, maxBinY);
        counter++;
    }
}

std::unordered_map<std::string, std::vector<contactRecord>> FreeHiC::getData() {
    std::unordered_map<std::string, std::vector<contactRecord>> ans;
    for (const auto &it : this->recordsMap) {
        std::string key = it.first;
        const std::vector<FreeContact> &data = it.second;
        std::vector<contactRecord> tmp;
        tmp.reserve(data.size() * 3);
        for (const auto &fc : data) {
            const std::vector<contactRecord> neighbor = fc.getNeighborRecords();
            if (neighbor.size() > 0) {
                for (auto &cr : neighbor) {
                    if (cr.counts < 0.5)
                        continue;
                    tmp.emplace_back(cr);
                }
            }
            if (fc.getSampleCounts() < 0.5)
                continue;
            contactRecord cr;
            cr.binX = fc.getBinX();
            cr.binY = fc.getBinY();
            cr.counts = fc.getSampleCounts();
            tmp.emplace_back(cr);
        }
        std::sort(tmp.begin(), tmp.end());
        std::vector<contactRecord>::iterator uit =
            std::unique(tmp.begin(), tmp.end());
        tmp.resize(std::distance(tmp.begin(), uit));
        ans[key] = tmp;
    }
    return ans;
}
} // namespace FreeHiC
