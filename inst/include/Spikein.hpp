#ifndef HICLIB_SPIKEIN_HPP
#define HICLIB_SPIKEIN_HPP

#include "common.hpp"
#include "kdtree.hpp"

namespace FreeHiC {

using kdtree2::KDTreeResultVector;

class KernelSmooth {
  public:
    KernelSmooth(){};
    ~KernelSmooth() = default;

    inline int abs(const int &x) { return x > 0 ? x : -x; }

    inline int dist(const contactRecord &x1, const contactRecord &x2) {
        return abs(x1.binX - x2.binX) + abs(x1.binY - x2.binY);
    }

    std::vector<float> fit(const std::vector<contactRecord> &records,
                           const int &resolution) {

        size_t N = records.size();
        std::vector<float> f(N, 0.0);

        for (size_t i = 0; i < N; i++) {
            double sum = 0.0;
            for (size_t k = 0; k < N; k++) {
                int d = dist(records[i], records[k]) / resolution;
                if (d > thresh)
                    continue;

                f[i] += weights[d] * records[k].counts;
                sum += weights[d];
            }

            f[i] /= sum;
        }
        return f;
    }

    std::vector<float> fit(const std::vector<contactRecord> &background,
                           const KDTreeResultVector &recordIdx,
                           const int &resolution) {

        size_t N = recordIdx.size();
        std::vector<float> f(N, 0.0);

        for (size_t i = 0; i < N; i++) {
            double sum = 0.0;
            for (size_t k = 0; k < N; k++) {
                int d = dist(background[recordIdx[i].idx],
                             background[recordIdx[k].idx]) /
                        resolution;
                if (d > thresh)
                    continue;

                f[i] += weights[d] * background[recordIdx[k].idx].counts;
                sum += weights[d];
            }

            f[i] /= sum;
        }
        return f;
    }

  private:
    const double weights[5] = {0.5703485, 0.3459338, 0.07718827, 0.06335999,
                              0.001933};
    const int thresh = 5;
};

class Spikein {
  public:
    Spikein(std::vector<contactRecord> &background, int resolution, bool smooth)
        : background_(background), resolution_(resolution),
          kdt(background, true) {
        kdt.sort_results = true;
    }

    ~Spikein() = default;

    const std::vector<contactRecord> getData() const {return background_;}

  private:
    std::vector<contactRecord> background_;
    int resolution_;
    KernelSmooth ks;
    kdtree2::KDTree kdt;

    KDTreeResultVector findIndex(contactRecord &trueSignal) {
        KDTreeResultVector signalIndex;
        signalIndex.clear();
        kdt.n_nearest(trueSignal, 1, signalIndex); // find the cloest one
        return signalIndex;
    }

    KDTreeResultVector findSignalNeighbor(contactRecord &trueSignal) {
        KDTreeResultVector signalNeighborIndex;
        signalNeighborIndex.clear();
        kdt.r_nearest(trueSignal, resolution_ * 8, signalNeighborIndex);
        return signalNeighborIndex;
    }

  public:
    void smoothSignal(std::vector<contactRecord> &trueSignal) {
        size_t N = trueSignal.size();
        double avg = 0.0;
        for (size_t i = 0; i < N; i++) {
            KDTreeResultVector selfIndex = findIndex(trueSignal[i]);
            KDTreeResultVector neighborIndex =
                findSignalNeighbor(trueSignal[i]);
            std::vector<contactRecord> block(neighborIndex.size());
            avg += neighborIndex.size();
            for (size_t b = 0; b < neighborIndex.size(); b++) {
                if (neighborIndex[b].idx == selfIndex[b].idx) {
                    block[b] = trueSignal[b];
                } else {
                    block[b] = this->background_[neighborIndex[b].idx];
                }
            }
            std::vector<float> smoothCount = ks.fit(block, this->resolution_ / 2);
            for (size_t b = 0; b < neighborIndex.size(); b++) {
                if (neighborIndex[b].idx == selfIndex[b].idx) smoothCount[b] = trueSignal[b].counts;
                this->background_[neighborIndex[b].idx].counts = smoothCount[b];
            }
        }
        if (debug) {        
            cout << "average neighbor: " << avg / (double) N << endl;
        }
    }
};

} // namespace FreeHiC

#endif