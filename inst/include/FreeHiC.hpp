//
// Created by Peigen Zhou on 5/23/20.
//

#ifndef HICLIB_FREEHIC_H
#define HICLIB_FREEHIC_H

#include "common.hpp"
// #include "hicReader.hpp"

namespace FreeHiC {

//    static std::default_random_engine generator;
//    static std::bernoulli_distribution bernoulli(0.5);
//    static std::uniform_real_distribution<double> uniform(0, 1);
    

    struct FreeInfo {
        int totalCounts = 0.0;
        int maxBinX = 0;
        int maxBinY = 0;
        std::string key;
    };

    class FreeContact : public contactRecord {
    private:
        // count information
        double sampleCounts;

        double neighborAverageCounts;
        double neighborZeroProb;
        double area;

        void probIncrement();

        void probDecrement();

        double getBinomial(double n, double p);

        RandomNumber *rng;
        std::vector<contactRecord> neighborRecords;
    public:
        FreeContact();

        ~FreeContact() = default;

        FreeContact(const int &x, const int &y, const double &counts) : contactRecord(x, y, counts),
                                                                        sampleCounts(0.0) {}

        FreeContact(const contactRecord &contact) : contactRecord(contact.binX, contact.binY, contact.counts),
                                                    sampleCounts(0.0) {}

        void noise() {
            if (rng->bernoulli()) probIncrement();
            else probDecrement();
        }

        void updateNeighborZero(const int &prevX, const int &prevY, const double &prevCounts,
                                const int &maxX, const int &maxY,
                                const int &resolution, const double &totalScaledCounts);

        void sample(const double &totalCounts, const double &totalScaleCounts);

        const double getCounts() const { return this->counts; }

        const double getSampleCounts() const { return this->sampleCounts; };

        const int getBinX() const { return this->binX; };

        const int getBinY() const { return this->binY; };

        const double getNeighborCounts() const { return neighborAverageCounts; }

        const double getNeighbor() const { return area; }

        const std::vector<contactRecord> getNeighborRecords() const { return this->neighborRecords; };

        void assignRng(RandomNumber *rng_);

    };

    class FreeHiC {
    public:
        FreeHiC();

        ~FreeHiC() = default;

        FreeHiC(const int &resolution, const double &sequenceDepth,
                const double &countScale, const double &noiseRate,
                const double &neighborZeroRate) : resolution_(resolution), sequenceDepth_(sequenceDepth),
                                                  countScale_(countScale), noiseRate_(noiseRate),
                                                  neighborZeroRate_(neighborZeroRate) {};

        void simulate(const std::unordered_map<std::string, std::vector<contactRecord>> &data);
        void simulate(const std::vector<contactRecord> &data);

        std::unordered_map<std::string, std::vector<contactRecord>> getData();

        RandomNumber rng;
    private:
        // Data setting
        int resolution_;

        // simulation setting
        double sequenceDepth_;
        double countScale_;

        double noiseRate_;
        double neighborZeroRate_;

        // data information
        std::unordered_map<std::string, std::pair<int, int>> pairMaxBin;
        // std::unordered_map<std::string, double> pairArea;
        // std::unordered_map<std::string, double> pairAverageCounts;

        std::unordered_map<std::string, std::vector<FreeContact>> recordsMap;

        double totalCounts = 0.0;

        void getBasicInformation(const std::unordered_map<std::string, std::vector<contactRecord>> &dataMap);
        std::vector<FreeContact> getBasicInformation(const std::vector<contactRecord> &data, const std::string& key, FreeInfo &info);
    };
}


#endif //HICLIB_FREEHIC_H
