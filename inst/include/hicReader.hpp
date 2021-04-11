#include <climits>

#ifndef HICLIB_HICREADER_HPP
#define HICLIB_HICREADER_HPP

#include "common.hpp"

// First write class for local file

namespace FreeHiC {
// sparse matrix entry

namespace Juicer {

    typedef long long ll;

    /*
    The MIT License (MIT)

    Copyright (c) 2016 Aiden Lab

    Permission is hereby granted, free of charge, to any person obtaining a copy
            of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
            to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
            copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
            copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
            AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
    */
// pointer structure for reading blocks or matrices,
// holds the size and position
struct indexEntry {
    int size;
    ll position;
};

struct FragIndexEntry {
    int nSites;
    ll position;
    FragIndexEntry() = default;
    FragIndexEntry(int nSites_, int position_)
        : nSites(nSites_), position(position_) {}
};

struct membuf : std::streambuf {
    membuf(char *begin, char *end) { this->setg(begin, begin, end); }
};

// for holding data from URL call
struct MemoryStruct {
    char *memory;
    size_t size;
};

// chromosome
struct chromosome {
    std::string name;
    int index;
    int length;
};

class hicReader {

  public:
    hicReader();

    hicReader(const std::string &fileName, const std::string &unit,
              const int resolution)
        : fileName_(fileName), unit_(unit), resolution_(resolution) {}

    ~hicReader() = default;

    std::vector<contactRecord> readData(const std::string &chr1,
                                        const std::string &chr2) {
        this->chr1 = cleanUpName(chr1);
        this->chr2 = cleanUpName(chr2);
        if (!this->initialized) {
            this->init();
        }
        if (!this->checked)
            return this->records;
        else {
            if (!this->readData_(this->chr1, this->chr2)) {
                cerr << "Please check the error message." << endl;
            }
        }
        return this->records;
    }

    const std::string getGenome() const { return this->genome; }

    const std::map<std::string, std::vector<int>> getResolutions() const {
        return this->fileResolutions;
    };

    std::vector<std::string> getPairs() {
        std::vector<std::string> ans;
        for (const auto &it : pairFilePositions) {
            if (it.first == "0_0") continue;
            ans.push_back(it.first);
        }
        return ans;
    }

    std::map<std::string, int> getChromosomeSize() {
        std::map<std::string, int> ans;
        for (const auto &it : this->chromosomeMap) {
            ans[it.first] = it.second.length;
        }
        return ans;
    };

    std::vector<int> getChromosomeSites(std::string chr) {
        if (!this->initialized) {
            this->init();
        }
        std::string chrClean = cleanUpName(chr);
        if (this->fragmentSitesCache.find(chrClean) ==
            this->fragmentSitesCache.end()) {
            getChromosomeSites_(chrClean);
        }
        std::sort(this->fragmentSitesCache[chrClean].begin(),
                  this->fragmentSitesCache[chrClean].end());
        return this->fragmentSitesCache[chrClean];
    }

    bool init() {
        if (this->initialized)
            return true;
        this->initialized = true;
        this->checked = this->prepareData();
        this->records.clear();
        if (!this->checked)
            cerr << "please checked the input and error message." << endl;
        return this->checked;
    }

    bool closeFile();
  protected:
    // hic file
    std::string fileName_;
    std::string chr1, chr2;
    std::string unit_ = "BP";
    int resolution_ = 5000;
    bool useRegionIndex = false;
    ll filePos = 0;

    std::ifstream fileInStream;

    // file attributes
    int version = 0;
    ll master = 0;

    std::map<std::string, std::vector<int>> fileResolutions;
    std::map<std::string, ll> pairFilePositions;

    // fragment
    bool needFragmentSite = true;
    ll fragmentSitePos = 0;
    std::map<std::string, FragIndexEntry> fragmentSitesIndex;
    std::map<std::string, std::vector<int>> fragmentSitesCache;
    std::map<std::string, int> chromosomeNSites; // fragmentCounts

    // hic info
    std::string genome;
    std::map<std::string, chromosome> chromosomeMap;
    std::vector<std::string> chromosomes;
    std::map<int, indexEntry> blockMap;
    std::vector<contactRecord> records;

    bool checked = false;
    bool initialized = false;

    virtual bool readData_(const std::string &chr1, const std::string &chr2);

    virtual bool prepareData();

    bool openFile();

    static bool readMagicString(std::istream &fin);

    bool readHeader(std::istream &fin);

    bool readFooter(std::istream &fin);

    bool readMatrix(std::istream &fin, ll filePosition, int &blockBinCount,
                    int &blockColumnCount);

    bool readMatrixZoomData(std::istream &fin, int &blockBinCount,
                            int &blockColumnCount);

    //            bool MatrixZoomData(int blockBinCount, int blockColumnCount,
    //            const std::vector<int> &chr1Sites, const std::vector<int>
    //            &chr2Sites);
    //
    std::vector<int> readSites(std::istream &fin, int nSites);
    std::vector<contactRecord> readBlock(char *compressedBytes, indexEntry idx);

    std::unordered_set<int>
    getBlockNumbersForRegionFromBinPosition(int *regionIndices,
                                            int blockBinCount,
                                            int blockColumnCount, bool intra);

    std::unordered_set<int> getBlockNumbers();
    bool findResolution();

    virtual void getChromosomeSites_(std::string chr);
    std::string cleanUpName(std::string name) {
        trim(name);
        toLower(name);
        replace(name, "chr", "");
        return name;
    }
};

// Read file from web

class hicReaderHttp : public hicReader {
  public:
    hicReaderHttp();

    hicReaderHttp(const std::string &fileName, const std::string &unit,
                  const int resolution)
        : hicReader(fileName, unit, resolution) {
        this->needFragmentSite = false;
    }

    ~hicReaderHttp() = default;

  private:
    CURL *urlBuffer = nullptr;
    // ll totalBytes;

    bool prepareData() override;

    bool readData_(const std::string &chr1, const std::string &chr2) override;

    static size_t WriteMemoryCallback(void *contents, size_t size, size_t nmemb,
                                      void *userp);

    char *getData(CURL *curl, ll position, int chunksize);

    CURL *initCURL(const char *url);

    bool readMatrix(CURL *curl, ll filePosition, int &blockBinCount,
                        int &blockColumnCount);

    bool readMatrixZoomData(CURL *curl, ll &FilePosition,
                                int &blockBinCount, int &blockColumnCount);
    void getChromosomeSites_(std::string chr) override;
};

} // namespace Juicer

} // namespace FreeHiC
#endif // !HICLIB_HICREADER_HPP
