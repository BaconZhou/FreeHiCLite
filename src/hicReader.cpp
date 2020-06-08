#include "hicReader.hpp"
#include <chrono>
#include <ctime>

#define log2 0.6931472
namespace FreeHiC {
    namespace Juicer {

        char *time() {
            auto end = std::chrono::system_clock::now();
            std::time_t end_time = std::chrono::system_clock::to_time_t(end);
            return std::ctime(&end_time);
        }

        std::unordered_set<int>
        hicReader::getBlockNumbersForRegionFromBinPosition(int *regionIndices, int blockBinCount, int blockColumnCount,
                                                           bool intra) {
            int col1 = regionIndices[0] / blockBinCount;
            int col2 = (regionIndices[1] + 1) / blockBinCount;
            int row1 = regionIndices[2] / blockBinCount;
            int row2 = (regionIndices[3] + 1) / blockBinCount;
            std::unordered_set<int> blocksSet;
            // first check the upper triangular matrix
            for (int r = row1; r <= row2; r++) {
                for (int c = col1; c <= col2; c++) {
#ifdef RVERSION
                    Rcpp::checkUserInterrupt();
#endif
                    int blockNumber = r * blockColumnCount + c;
                    blocksSet.insert(blockNumber);
                }
            }
            // check region part that overlaps with lower left triangle
            // but only if intra chromosomal
            if (intra) {
                for (int r = col1; r <= col2; r++) {
                    for (int c = row1; c <= row2; c++) {
#ifdef RVERSION
                        Rcpp::checkUserInterrupt();
#endif
                        int blockNumber = r * blockColumnCount + c;
                        blocksSet.insert(blockNumber);
                    }
                }
            }
            // std::vector<int> vec( blocksSet.begin(), blocksSet.end() );
            // std::sort( vec.begin(), vec.end() );

            return blocksSet;
        }

        bool hicReader::findResolution() {
            const std::vector<int> &resolution = this->fileResolutions[this->unit_];
            auto it = std::find(resolution.cbegin(), resolution.cend(), this->resolution_);
            return it != resolution.cend();
        }

        bool hicReader::prepareData() {
            if (!this->openFile()) return false;
            this->readHeader(this->fileInStream);

            this->fileInStream.seekg(master, std::ios::beg);
            bool foundFooter = this->readFooter(this->fileInStream);
            return foundFooter;
        }

        bool hicReader::readData_(const std::string &chr1, const std::string &chr2) {

            this->records.clear();
            if (!this->findResolution()) {
                cerr << "Not find this resolution: (" << this->unit_ << ", " << this->resolution_ << ")" << endl;
                cerr << "Available resolutions:" << endl;
                cerr << "BP: ";
                for (auto x : this->fileResolutions["BP"]) cerr << x << " ";
                cerr << endl;
                cerr << "FRAG: ";
                for (auto x : this->fileResolutions["FRAG"]) cerr << x << " ";
                cerr << endl;
                return false;
            }
            int c1 = std::min(this->chromosomeMap[chr1].index,
                              this->chromosomeMap[chr2].index);
            int c2 = std::max(this->chromosomeMap[chr1].index,
                              this->chromosomeMap[chr2].index);

            int binSize = this->resolution_;

            int origRegionIndices[4] = {
                    0, this->chromosomeMap[chr1].length, 0,
                    this->chromosomeMap[chr2].length}; // as given by user
            int regionIndices[4];                        //

            for (int i = 0; i < 4; i++)
                regionIndices[i] = origRegionIndices[i] / binSize;

            std::stringstream ss;
            ss << c1 << "_" << c2;
            std::string key = ss.str();

            if (this->pairFilePositions.find(key) == this->pairFilePositions.end()) {
                cerr << "Can not find this pair " << key << endl;
                return false;
            }
            ll filePosition = this->pairFilePositions[key];

            int blockBinCount, blockColumnCount;
            bool foundBlock = this->readMatrix(this->fileInStream, filePosition,
                                               blockBinCount, blockColumnCount);
            if (!foundBlock)
                return false;

            std::unordered_set<int> blockNumbers;
            if (useRegionIndex) {
                blockNumbers = this->getBlockNumbersForRegionFromBinPosition(
                        regionIndices, blockBinCount, blockColumnCount, c1 == c2);
            } else {
                blockNumbers = this->getBlockNumbers();
            }

            std::vector<int> blockVec(blockNumbers.begin(), blockNumbers.end());
            std::sort(blockVec.begin(), blockVec.end());

            std::vector<contactRecord> tmp_records;
            for (const int & blockNumber : blockVec) {
                // get contacts in this block
#ifdef RVERSION
                Rcpp::checkUserInterrupt();
#endif
                char *compressedBytes = new char[blockMap[blockNumber].size];
                this->fileInStream.seekg(blockMap[blockNumber].position, std::ios::beg);
                this->fileInStream.read(compressedBytes, blockMap[blockNumber].size);
                tmp_records = this->readBlock(compressedBytes, blockMap[blockNumber]);
                for (auto rec : tmp_records) {
                    int x = rec.binX * this->resolution_;
                    int y = rec.binY * this->resolution_;
                    float c = rec.counts;
//            if (this->norm != "NONE") {
//                 c = c / (c1Norm[rec.binX] * c2Norm[rec.binY]);
//            }
                    if ((x >= origRegionIndices[0] && x <= origRegionIndices[1] &&
                         y >= origRegionIndices[2] && y <= origRegionIndices[3]) ||
                        // or check regions that overlap with lower left
                        ((c1 == c2) && y >= origRegionIndices[0] &&
                         y <= origRegionIndices[1] && x >= origRegionIndices[2] &&
                         x <= origRegionIndices[3])) {
                        contactRecord record;
                        record.binX = x;
                        record.binY = y;
                        record.counts = c;
                        this->records.push_back(record);
                    }
                }
            }
            std::sort(this->records.begin(), this->records.end());
            return true;
        }

        bool hicReader::openFile() {
            this->fileInStream.open(this->fileName_, std::fstream::in | std::fstream::binary);
            if (!this->fileInStream) {
                cerr << "File " << this->fileName_ << " cannot be opened for reading"
                     << endl;
                return false;
            }
            return true;
        }

        bool hicReader::readMagicString(std::istream &fin) {
            std::string str;
            getline(fin, str, '\0');
            return str[0] == 'H' && str[1] == 'I' && str[2] == 'C';
        }

        bool hicReader::readHeader(std::istream &fin) {
#ifdef WINTEST
            DEBUG(this->filePos);
#endif
            if (!this->readMagicString(fin)) {
                cerr << "Hi-C magic string is missing, does not appear to be a hic file"
                     << endl;
                this->master = -1;
                return false;
            }
#ifdef WINTEST
            DEBUG(this->filePos);
#endif

            this->filePos += 4;
            fin.read((char *) &this->version, sizeof(int));
#ifdef WINTEST
            DEBUG(this->version);
            DEBUG(this->filePos);
#endif

            if (version < 6) {
                cerr << "Version " << version << " no ller supported" << endl;
                this->master = -1;
                return false;
            }
            this->filePos += 4;

            fin.read((char *) &this->master, sizeof(ll));
#ifdef WINTEST
            DEBUG(this->master);
            DEBUG(sizeof(ll));
            DEBUG(this->filePos);
#endif
            this->filePos += 8;
            getline(fin, this->genome, '\0');

            this->filePos += this->genome.size() + 1;
#ifdef WINTEST
            DEBUG(this->genome);
            DEBUG(this->genome.size());
            DEBUG(this->filePos);
#endif
            int nAttributes;
            fin.read((char *) &nAttributes, sizeof(int));
#ifdef WINTEST
            DEBUG(nAttributes);
            DEBUG(this->filePos);
#endif

            this->filePos += 4;

            for (int i = 0; i < nAttributes; i++) {
                std::string key, value;
                getline(fin, key, '\0');
                this->filePos += key.size() + 1;
                getline(fin, value, '\0');
                this->filePos += value.size() + 1;
#ifdef WINTEST
                DEBUG(key.size());
                DEBUG(value.size());
                DEBUG(this->filePos);
#endif
            }

            int nChrs;
            fin.read((char *) &nChrs, sizeof(int));
            this->filePos += 4;
#ifdef WINTEST
            DEBUG(nChrs);
            DEBUG(this->filePos);
#endif
            if (nChrs > 1000) {
                cerr << "It may contains some bug. Please create an issue in the github. Thanks!" << endl;
#ifdef RVERSION
                Rcpp::stop("wired nchrs");
#else
                return false;
#endif
            }
            for (int i = 0; i < nChrs; i++) {
                std::string name;
                int length;
                getline(fin, name, '\0');

                this->filePos += name.size()+1;
                name = cleanUpName(name);
                fin.read((char *) &length, sizeof(int));

                this->filePos += 4;

                chromosome chr;
                chr.index = i;
                chr.name = name;
                chr.length = length;
                this->chromosomeMap[name] = chr;
                this->chromosomes.push_back(name);
            }

            int nBpResolution;
            fin.read((char *) &nBpResolution, sizeof(int));
            this->filePos += 4;
            if (nBpResolution > 10) {
                cerr << "Current number of bp resolution is " << nBpResolution
                    <<  ". It should smaller than 10" << endl;
                cerr << "It may contains some bug. Please create an issue in the github. Thanks!" << endl;
#ifdef RVERSION
                Rcpp::stop("wired bp resolution");
#else
                return false;
#endif
            }
            std::vector<int> bpResolution;
            for (int i = 0; i < nBpResolution; i++) {
                int bpResolution_;
                fin.read((char *) &bpResolution_, sizeof(int));
                this->filePos += 4;
                bpResolution.push_back(bpResolution_);
            }
            this->fileResolutions["BP"] = bpResolution;

            int nFragResolution;
            fin.read((char *) &nFragResolution, sizeof(int));
            this->filePos += 4;
            if (nFragResolution > 9) {
                cerr << "Current number of FRAG resolution is " << nFragResolution
                     <<  ". It should smaller than 9" << endl;
                cerr << "It may contains some bug. Please create an issue in the github. Thanks!" << endl;
#ifdef RVERSION
                Rcpp::stop("wired frag resolution");
#else
                return false;
#endif
            }

            std::vector<int> FragResolution;
            for (int i = 0; i < nFragResolution; i++) {
                int FragResolution_;
                fin.read((char *) &FragResolution_, sizeof(int));
                this->filePos += 4;
                FragResolution.push_back(FragResolution_);
            }
            this->fileResolutions["FRAG"] = FragResolution;
            this->fragmentSitePos = this->filePos;
            if (nFragResolution > 0 && this->needFragmentSite) {
                for (int i = 0; i < nChrs; i++) {
                    std::string chr = chromosomes[i];
                    int nSite_;
                    fin.read((char *) &nSite_, sizeof(int));
                    this->fragmentSitesIndex[chr] = FragIndexEntry(nSite_, this->filePos);
#ifdef HTTPTEST
                    cerr << "chr: " << chr << ", nSites_: " << nSite_ << ", pos: " << this->filePos << endl;
#endif
                    this->chromosomeNSites[chr] = nSite_;
                    fin.ignore(nSite_ * 4);
                    this->filePos += nSite_ * 4;
                }
            }
            return true;
        }

        bool hicReader::readFooter(std::istream &fin) {
            int nBytes;
            fin.read((char *) &nBytes, sizeof(int));
            int nEntries;
            fin.read((char *) &nEntries, sizeof(int));
            bool found = false;
            for (int i = 0; i < nEntries; i++) {
#ifdef RVERSION
                Rcpp::checkUserInterrupt();
#endif
                std::string str;
                getline(fin, str, '\0');
                ll fpos;
                fin.read((char *) &fpos, sizeof(ll));
                int sizeInBytes;
                fin.read((char *) &sizeInBytes, sizeof(int));
                this->pairFilePositions[str] = fpos;
                found = true;
            }
            if (!found) {
                cerr << "File doesn't have any contact pair map." << endl;
            }
            return found;
        }

        bool hicReader::readMatrix(std::istream &fin, ll filePosition,
                                   int &blockBinCount, int &blockColumnCount) {
            fin.seekg(filePosition, std::ios::beg);
            int c1, c2;
            fin.read((char *) &c1, sizeof(int)); // chr1
            fin.read((char *) &c2, sizeof(int)); // chr2
            int nRes;
            fin.read((char *) &nRes, sizeof(int));


            int i = 0;
            bool found = false;
            //while (i < nRes && !found) {
            while(i < nRes) {
                bool found_ = this->readMatrixZoomData(fin, blockBinCount, blockColumnCount);
                if (found_) found = true;
                i++;
            }
            if (!found) {
                cerr << "Error finding block data. The given chromosomes pair not find: (" << chr1 << ", " << chr2 << ")" << endl;
                cerr << "Available chromosomes are: \n";
                for (const std::string &chr : chromosomes) {
                    if (chr == "ALL") continue;
                    cerr << chr << " ";
                }
                cerr << endl;
            }
            return found;
        }

        bool hicReader::readMatrixZoomData(std::istream &fin, int &myBlockBinCount,
                                           int &myBlockColumnCount) {
            std::string unit;
            getline(fin, unit, '\0'); // unit
            int tmp;
            fin.read((char *) &tmp, sizeof(int)); // Old "zoom" index -- not used
            float tmp2;
            fin.read((char *) &tmp2, sizeof(float)); // sumCounts
            fin.read((char *) &tmp2, sizeof(float)); // occupiedCellCount
            fin.read((char *) &tmp2, sizeof(float)); // stdDev
            fin.read((char *) &tmp2, sizeof(float)); // percent95
            int binSize;
            fin.read((char *) &binSize, sizeof(int));
            int blockBinCount;
            fin.read((char *) &blockBinCount, sizeof(int));
            int blockColumnCount;
            fin.read((char *) &blockColumnCount, sizeof(int));

            bool found = false;

            if (this->unit_ == unit && this->resolution_ == binSize) {
                myBlockBinCount = blockBinCount;
                myBlockColumnCount = blockColumnCount;
                found = true;
            }

            int nBlocks;
            fin.read((char *) &nBlocks, sizeof(int));

            std::unordered_map<int, indexEntry> blockMap_;
            for (int b = 0; b < nBlocks; b++) {

#ifdef RVERSION
                Rcpp::checkUserInterrupt();
#endif
                int blockNumber;
                fin.read((char *) &blockNumber, sizeof(int));
                ll filePosition;
                fin.read((char *) &filePosition, sizeof(ll));
                int blockSizeInBytes;
                fin.read((char *) &blockSizeInBytes, sizeof(int));
                indexEntry entry{};
                entry.size = blockSizeInBytes;
                entry.position = filePosition;
                blockMap_[blockNumber] = entry;
                if (found) blockMap[blockNumber] = entry;
            }
            return found;
        }

        std::vector<contactRecord> hicReader::readBlock(char *compressedBytes,
                                                        indexEntry idx) {
            if (idx.size == 0) {
                std::vector<contactRecord> v;
                return v;
            }
            char *uncompressedBytes =
                    new char[idx.size * 10]; // biggest seen so far is 3

            // Decompress the block
            // zlib struct
            z_stream infstream;
            infstream.zalloc = Z_NULL;
            infstream.zfree = Z_NULL;
            infstream.opaque = Z_NULL;
            infstream.avail_in = (uInt) (idx.size);           // size of input
            infstream.next_in = (Bytef *) compressedBytes;    // input char array
            infstream.avail_out = (uInt) idx.size * 10;       // size of output
            infstream.next_out = (Bytef *) uncompressedBytes; // output char array
            // the actual decompression work.
            inflateInit(&infstream);
            inflate(&infstream, Z_NO_FLUSH);
            inflateEnd(&infstream);
            int uncompressedSize = infstream.total_out;

            // create stream from buffer for ease of use
            membuf sbuf(uncompressedBytes, uncompressedBytes + uncompressedSize);
            std::istream bufferin(&sbuf);
            int nRecords;
            bufferin.read((char *) &nRecords, sizeof(int));
            std::vector<contactRecord> v(nRecords);
            // different versions have different specific formats
            if (this->version < 7) {
                for (int i = 0; i < nRecords; i++) {
                    int binX, binY;
                    bufferin.read((char *) &binX, sizeof(int));
                    bufferin.read((char *) &binY, sizeof(int));
                    float counts;
                    bufferin.read((char *) &counts, sizeof(float));
                    contactRecord record;
                    record.binX = binX;
                    record.binY = binY;
                    record.counts = counts;
                    v[i] = record;
                }
            } else {
                int binXOffset, binYOffset;
                bufferin.read((char *) &binXOffset, sizeof(int));
                bufferin.read((char *) &binYOffset, sizeof(int));
                char useShort;
                bufferin.read((char *) &useShort, sizeof(char));
                char type;
                bufferin.read((char *) &type, sizeof(char));
                int index = 0;
                if (type == 1) {
                    // List-of-rows representation
                    short rowCount;
                    bufferin.read((char *) &rowCount, sizeof(short));
                    for (int i = 0; i < rowCount; i++) {
#ifdef RVERSION
                        Rcpp::checkUserInterrupt();
#endif
                        short y;
                        bufferin.read((char *) &y, sizeof(short));
                        int binY = y + binYOffset;
                        short colCount;
                        bufferin.read((char *) &colCount, sizeof(short));
                        for (int j = 0; j < colCount; j++) {
#ifdef RVERSION
                            Rcpp::checkUserInterrupt();
#endif
                            short x;
                            bufferin.read((char *) &x, sizeof(short));
                            int binX = binXOffset + x;
                            float counts;
                            if (useShort == 0) { // yes this is opposite of usual
                                short c;
                                bufferin.read((char *) &c, sizeof(short));
                                counts = c;
                            } else {
                                bufferin.read((char *) &counts, sizeof(float));
                            }
                            contactRecord record;
                            record.binX = binX;
                            record.binY = binY;
                            record.counts = counts;
                            v[index] = record;
                            index++;
                        }
                    }
                } else if (type == 2) {
                    int nPts;
                    bufferin.read((char *) &nPts, sizeof(int));
                    short w;
                    bufferin.read((char *) &w, sizeof(short));

                    for (int i = 0; i < nPts; i++) {
#ifdef RVERSION
                        Rcpp::checkUserInterrupt();
#endif
                        // int idx = (p.y - binOffset2) * w + (p.x - binOffset1);
                        int row = i / w;
                        int col = i - row * w;
                        int bin1 = binXOffset + col;
                        int bin2 = binYOffset + row;

                        float counts;
                        if (useShort == 0) { // yes this is opposite of the usual
                            short c;
                            bufferin.read((char *) &c, sizeof(short));
                            if (c != -32768) {
                                contactRecord record;
                                record.binX = bin1;
                                record.binY = bin2;
                                record.counts = c;
                                v[index] = record;
                                index++;
                            }
                        } else {
                            bufferin.read((char *) &counts, sizeof(float));
                            if (!isnan(counts)) {
                                contactRecord record;
                                record.binX = bin1;
                                record.binY = bin2;
                                record.counts = counts;
                                v[index] = record;
                                index++;
                            }
                        }
                    }
                }
            }
            delete[] compressedBytes;
            delete[] uncompressedBytes; // don't forget to delete your heap arrays in
            // C++!
            return v;
        }

//        int HiCFragmentAxis(int binSize, std::vector<int> sites, int chrLength) {
//            double averageBinSizeInBP = ((double) chrLength) / (sites.size() + 1) * binSize;
//            int igvZoom = (int) ((std::log((double) chrLength / 700.0) / averageBinSizeInBP) / log2);
//            return igvZoom;
//        }
//
//        bool hicReader::MatrixZoomData(int blockBinCount, int blockColumnCount, const std::vector<int> &chr1Sites, const std::vector<int> &chr2Sites) {
//            int correctedBinCount = blockBinCount;
//            if (this->version < 8 && this->chromosomeMap[chr1].length < this->chromosomeMap[chr2].length) {
//                bool isFrag = unit_ == "FRAG";
//                int len1 = this->chromosomeMap[chr1].length;
//                int len2 = this->chromosomeMap[chr2].length;
//
//                if (chr1Sites.size() > 0 && chr2Sites.size() > 0 && isFrag) {
//                    len1 = chr1Sites.size() + 1;
//                    len2 = chr2Sites.size() + 1;
//                }
//                int nBinsX = std::max(len1, len2) / resolution_ + 1;
//                correctedBinCount = nBinsX / blockColumnCount + 1;
//            }
//
//            return false;
//        }
//
        std::vector<int> hicReader::readSites(std::istream &fin, int nSites) {
            std::vector<int> sites(nSites);
            for (int i = 0; i < nSites; i++) {
                int site;
                fin.read((char *) &site, sizeof(int));
                sites[i] = site;
            }
            return sites;
        }

        // ll total_bytes;

        size_t hdf(char *b, size_t size, size_t nitems, void *userdata) {
            size_t numbytes = size * nitems;
            b[numbytes + 1] = '\0';
            std::string s(b);
            int found = s.find("Content-Range");
            if (found != std::string::npos) {
                int found2 = s.find("/");
                //Content-Range: bytes 0-100000/891471462
                if (found2 != std::string::npos) {
                    std::string total = s.substr(found2 + 1);
                    // total_bytes = stol(total);
                }
            }
            return numbytes;
        }

        size_t hicReaderHttp::WriteMemoryCallback(void *contents, size_t size, size_t nmemb, void *userp) {
            size_t realsize = size * nmemb;
            struct MemoryStruct *mem = (struct MemoryStruct *) userp;

            mem->memory = static_cast<char *>(realloc(mem->memory, mem->size + realsize + 1));
            if (mem->memory == nullptr) {
                /* out of memory! */
                cerr << "not enough memory (realloc returned NULL)\n";
                return 0;
            }

            std::memcpy(&(mem->memory[mem->size]), contents, realsize);
            mem->size += realsize;
            mem->memory[mem->size] = 0;

            return realsize;
        }

        char *hicReaderHttp::getData(CURL *curl, ll position, int chunksize) {
            std::ostringstream oss;
            struct MemoryStruct chunk;

            chunk.memory = static_cast<char *>(malloc(1));
            chunk.size = 0;    /* no data at this point */
            oss << position << "-" << position + chunksize;
            curl_easy_setopt(curl, CURLOPT_WRITEDATA, (void *) &chunk);
            curl_easy_setopt(curl, CURLOPT_RANGE, oss.str().c_str());
            CURLcode res = curl_easy_perform(curl);
            if (res != CURLE_OK) {
                cerr << "curl_easy_perform() failed: " << curl_easy_strerror(res) << endl;
            }
            return chunk.memory;
        }

        CURL *hicReaderHttp::initCURL(const char *url) {
            CURL *curl = curl_easy_init();
            if (curl) {
                curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, this->WriteMemoryCallback);
                curl_easy_setopt(curl, CURLOPT_URL, url);
                //curl_easy_setopt (curl, CURLOPT_VERBOSE, 1L);
                curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
                curl_easy_setopt(curl, CURLOPT_HEADERFUNCTION, hdf);
                curl_easy_setopt(curl, CURLOPT_USERAGENT, "straw");
            }
            return curl;
        }

        bool hicReaderHttp::readMatrix(CURL *curl, ll filePosition, int &blockBinCount, int &blockColumnCount) {
#ifdef HTTPTEST
            cout << "Read matrix data " << time() << endl;
#endif
            char *buffer;
            int size = sizeof(int) * 3;
            buffer = getData(curl, filePosition, size);
            membuf sbuf(buffer, buffer + size);
            std::istream bufin(&sbuf);

            int c1, c2;
            bufin.read((char *) &c1, sizeof(int)); //chr1
            bufin.read((char *) &c2, sizeof(int)); //chr2
            int nRes;
            bufin.read((char *) &nRes, sizeof(int));
            int i = 0;
            bool found = false;
            filePosition = filePosition + size;
            delete buffer;

            while (i < nRes && !found) {
#ifdef RVERSION
                Rcpp::checkUserInterrupt();
#endif
                // myFilePosition gets updated within call
                found = this->readMatrixZoomData(curl, filePosition, blockBinCount, blockColumnCount);
                i++;
            }

            if (!found) {
                cerr << "Error finding block data. The given chromosomes pair not find: (" << chr1 << ", " << chr2 << ")" << endl;
                cerr << "Available chromosomes are: \n";
                for (const std::string &chr : chromosomes) {
                    if (chr == "ALL") continue;
                    cerr << chr << " ";
                }
                cerr << endl;
            }
#ifdef HTTPTEST
        cout << "Finished Read matrix data" << time() << endl;
#endif
            return found;
        }

//
        bool hicReaderHttp::readMatrixZoomData(CURL *curl, ll &FilePosition, int &myBlockBinCount,
                                                   int &myBlockColumnCount) {
            char *buffer;
            int header_size = 5 * sizeof(int) + 4 * sizeof(float);
            char *first;
            first = this->getData(curl, FilePosition, 1);
            if (first[0] == 'B') {
                header_size += 3;
            } else if (first[0] == 'F') {
                header_size += 5;
            } else {
                cerr << "Unit not understood" << endl;
                return false;
            }
            buffer = this->getData(curl, FilePosition, header_size);
            membuf sbuf(buffer, buffer + header_size);
            std::istream fin(&sbuf);

            std::string unit;
            getline(fin, unit, '\0'); // unit
            int tmp;
            fin.read((char *) &tmp, sizeof(int)); // Old "zoom" index -- not used
            float tmp2;
            fin.read((char *) &tmp2, sizeof(float)); // sumCounts
            fin.read((char *) &tmp2, sizeof(float)); // occupiedCellCount
            fin.read((char *) &tmp2, sizeof(float)); // stdDev
            fin.read((char *) &tmp2, sizeof(float)); // percent95
            int binSize;
            fin.read((char *) &binSize, sizeof(int));
            int blockBinCount;
            fin.read((char *) &blockBinCount, sizeof(int));
            int blockColumnCount;
            fin.read((char *) &blockColumnCount, sizeof(int));

            bool found = false;
            if (this->unit_ == unit && this->resolution_ == binSize) {
                myBlockBinCount = blockBinCount;
                myBlockColumnCount = blockColumnCount;
                found = true;
            }

            int nBlocks;
            fin.read((char *) &nBlocks, sizeof(int));

            if (found) {
                buffer = getData(curl, FilePosition + header_size,
                                 nBlocks * (sizeof(int) + sizeof(ll) + sizeof(int)));
                membuf sbuf2(buffer, buffer + nBlocks * (sizeof(int) + sizeof(ll) + sizeof(int)));
                std::istream fin2(&sbuf2);
                for (int b = 0; b < nBlocks; b++) {
#ifdef RVERSION
                    Rcpp::checkUserInterrupt();
#endif
                    int blockNumber;
                    fin2.read((char *) &blockNumber, sizeof(int));
                    ll filePosition;
                    fin2.read((char *) &filePosition, sizeof(ll));
                    int blockSizeInBytes;
                    fin2.read((char *) &blockSizeInBytes, sizeof(int));
                    indexEntry entry;
                    entry.size = blockSizeInBytes;
                    entry.position = filePosition;
                    blockMap[blockNumber] = entry;
                }
            } else {
                FilePosition = FilePosition + header_size + (nBlocks * (sizeof(int) + sizeof(ll) + sizeof(int)));
            }
            delete buffer;
            return found;
        }


        bool hicReaderHttp::prepareData() {
#ifdef HTTPTEST
            cout << "Init curl " << time() <<  endl;
#endif
            char *buffer = nullptr;
            this->urlBuffer = initCURL(this->fileName_.c_str());
#ifdef HTTPTEST
            cout << "Get buffer data " << time() << endl;
#endif
            ll headerReader = this->needFragmentSite ? 50000000 : 100000;
            if (this->urlBuffer) buffer = this->getData(this->urlBuffer, 0, headerReader);

#ifdef HTTPTEST
            cout << "Read head " << time() << endl;
#endif
            membuf sbuf(buffer, buffer + headerReader);
            std::istream bufin(&sbuf);
            this->readHeader(bufin);
            delete buffer;
            // this->totalBytes = total_bytes;
#ifdef HTTPTEST
            cout << "Find resolution " << time() << endl;
#endif
            if (!this->findResolution()) {
                cerr << "Not find this resolution: (" << this->unit_ << ", " << this->resolution_ << ")" << endl;
                cerr << "Available resolutions:" << endl;
                cerr << "BP: ";
                for (auto x : this->fileResolutions["BP"]) cerr << x << " ";
                cerr << endl;
                cerr << "FRAG: ";
                for (auto x : this->fileResolutions["FRAG"]) cerr << x << " ";
                cerr << endl;
                return false;
            }
#ifdef HTTPTEST
            cout << "get Data footer " << time() << endl;
#endif
            char *buffer2;

            ll bytesToRead = 100000; //this->totalBytes - this->master; // this mean download all the dataset
            // bool foundFooter = false;

            buffer2 = getData(this->urlBuffer, this->master, bytesToRead);
#ifdef HTTPTEST
            cout << "Read footer " << time() << endl;
#endif
            membuf sbuf2(buffer2, buffer2 + bytesToRead);
            std::istream bufin2(&sbuf2);
            bool foundFooter = this->readFooter(bufin2);

            delete buffer2;
            if (!foundFooter) {
                cerr << "Did not find footer" << endl;
                return false;
            }
#ifdef HTTPTEST
            cout << "Finished prepared " << time() << endl;
#endif
            return true;
        }

//
        bool hicReaderHttp::readData_(const std::string &chr1, const std::string &chr2) {

#ifdef HTTPTEST
            cout << "start read data " << time() << endl;
#endif
            int c1 = std::min(this->chromosomeMap[chr1].index,
                              this->chromosomeMap[chr2].index);
            int c2 = std::max(this->chromosomeMap[chr1].index,
                              this->chromosomeMap[chr2].index);
            int binSize = this->resolution_;

            int origRegionIndices[4] = {
                    0, this->chromosomeMap[chr1].length, 0,
                    this->chromosomeMap[chr2].length}; // as given by user
            int regionIndices[4];                        //

            for (int i = 0; i < 4; i++)
                regionIndices[i] = origRegionIndices[i] / binSize;

            int blockBinCount, blockColumnCount;
            std::stringstream ss;
            ss << c1 << "_" << c2;
            std::string key = ss.str();

            if (this->pairFilePositions.find(key) == this->pairFilePositions.end()) {
                cerr << "Can not find this pair " << key << endl;
                return false;
            }
            ll filePosition = this->pairFilePositions[key];


            bool foundBlock = this->readMatrix(this->urlBuffer, filePosition,
                                                   blockBinCount, blockColumnCount);
            if (!foundBlock) {
                cerr << "Did not find block" << endl;
                return false;
            }

            std::unordered_set<int> blockNumbers;
            if (useRegionIndex) {
                blockNumbers = this->getBlockNumbersForRegionFromBinPosition(
                        regionIndices, blockBinCount, blockColumnCount, c1 == c2);
            } else {
                blockNumbers = getBlockNumbers();

            }

            std::vector<int> blockVec(blockNumbers.begin(), blockNumbers.end());
            std::sort(blockVec.begin(), blockVec.end());
            std::vector<contactRecord> tmp_records;

            for (const int & blockNumber : blockVec) {
                // get contacts in this block
#ifdef RVERSION
                Rcpp::checkUserInterrupt();
#endif
                char *compressedBytes = new char[blockMap[blockNumber].size];
                compressedBytes = this->getData(this->urlBuffer, blockMap[blockNumber].position, blockMap[blockNumber].size);
                tmp_records = this->readBlock(compressedBytes, blockMap[blockNumber]);
                for (auto rec : tmp_records) {
                    int x = rec.binX * this->resolution_;
                    int y = rec.binY * this->resolution_;
                    float c = rec.counts;
//            if (this->norm != "NONE") {
//                 c = c / (c1Norm[rec.binX] * c2Norm[rec.binY]);
//            }
                    if ((x >= origRegionIndices[0] && x <= origRegionIndices[1] &&
                         y >= origRegionIndices[2] && y <= origRegionIndices[3]) ||
                        // or check regions that overlap with lower left
                        ((c1 == c2) && y >= origRegionIndices[0] &&
                         y <= origRegionIndices[1] && x >= origRegionIndices[2] &&
                         x <= origRegionIndices[3])) {
                        contactRecord record;
                        record.binX = x;
                        record.binY = y;
                        record.counts = c;
                        this->records.push_back(record);
                    }
                }
            }
            // auto end = std::chrono::system_clock::now();
            // std::chrono::duration<double> elapsed_seconds = end-start;
            // cout << "elapsed time: " << elapsed_seconds.count() << endl;
            std::sort(this->records.begin(), this->records.end());
            return true;
        }

        std::unordered_set<int> hicReader::getBlockNumbers() {
            std::unordered_set<int> ans;
            for (auto & it : blockMap) {
                ans.insert(it.first);
            }
            return ans;
        }

        void hicReader::getChromosomeSites_(std::string chr) {
            std::string chrClean = cleanUpName(chr);
            std::vector<int> ans;
            if (this->fragmentSitesCache.find(chrClean) != this->fragmentSitesCache.end()) {
                ans = this->fragmentSitesCache[chrClean];
            }
            if (ans.size() < 1) {
                if (this->fragmentSitesIndex.find(chrClean) != this->fragmentSitesIndex.end()) {
                    FragIndexEntry entry = this->fragmentSitesIndex[chrClean];
                    if (entry.nSites > 0) {
                        this->fileInStream.seekg(entry.position, std::ios::beg);
                        ans = readSites(this->fileInStream, entry.nSites);
                    }
                }
                this->fragmentSitesCache[chrClean] = ans;
            }
        }

        void hicReaderHttp::getChromosomeSites_(std::string chr) {
            std::string chrClean = cleanUpName(chr);
            std::vector<int> ans;
            if (this->fragmentSitesCache.find(chrClean) != this->fragmentSitesCache.end()) {
                ans = this->fragmentSitesCache[chrClean];
            }
            if (ans.size() < 1) {
                ll position = this->fragmentSitePos;
                if (fileResolutions["FRAG"].size() > 0) {
                    for (size_t i = 0; i < chromosomes.size(); i++) {
#ifdef RVERSION
                        Rcpp::checkUserInterrupt();
#endif
                        std::string chr = chromosomes[i];
                        if (this->fragmentSitesCache.find(chr) != this->fragmentSitesCache.end()) continue;
                        char *buffer = this->getData(this->urlBuffer, position, 4);
                        membuf sbuf(buffer, buffer + 4);
                        std::istream fin(&sbuf);
                        int nSite_;
                        fin.read((char *) &nSite_, sizeof(int));
                        delete buffer;

                        this->fragmentSitesIndex[chr] = FragIndexEntry(nSite_, position);
#ifdef HTTPTEST
                        cout << "chr: " << chr << ", nSite: " << nSite_ << ", position: " << position << endl;
#endif
                        position += 4;
                        char *buffer2 = this->getData(this->urlBuffer, position, nSite_ * 4 + 4);
                        membuf sbuf2(buffer2, buffer2 + nSite_ * 4 + 4);
                        std::istream fin2(&sbuf2);

                        this->chromosomeNSites[chr] = nSite_;
                        this->fragmentSitesCache[chr] = readSites(fin2, nSite_);
                        position += nSite_ * 4;
                        this->fragmentSitePos = position;
                        if (chrClean == chr) break;
                        delete buffer2;
                    }
                }
            }
        }

    } // namespace Juicer
} // namespace FreeHiC