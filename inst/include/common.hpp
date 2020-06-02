//
// Created by Peigen Zhou on 5/19/20.
//

#ifndef HICLIB_COMMON_HPP
#define HICLIB_COMMON_HPP

#include "math.h"
#include "zlib.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <curl/curl.h>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <streambuf>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <random>

#ifdef RVERSION
#include <Rcpp.h>
#define cerr Rcpp::Rcerr
#define cout Rcpp::Rcout
#else
#define cerr std::cerr
#define cout std::cout
#define HTTPTEST
#endif

#define debug false

#define DEBUG(x)                                                               \
    do {                                                                       \
        cerr << #x << ": " << x << endl;                                       \
    } while (0)

using std::endl;

namespace FreeHiC {

class RandomNumber {
  public:
    RandomNumber() : uni{0, 1.0}, bern{0.5} {};

    ~RandomNumber() = default;

    int bernoulli() { return this->bern(this->generator); }

    double uniform() { return this->uni(this->generator); }

  private:
    std::default_random_engine generator;
    std::uniform_real_distribution<double> uni;
    std::bernoulli_distribution bern;
};

class contactRecord {
  public:
    int binX;
    int binY;
    float counts;
    contactRecord() = default;
    contactRecord(int x, int y, float c) : binX(x), binY(y), counts(c) {}
    ~contactRecord() = default;

    int &operator[](int c) {
        if (c == 0)
            return binX;
        if (c == 1)
            return binY;
        return binX;
    }

    const int &operator[](int c) const {
        if (c == 0)
            return binX;
        if (c == 1)
            return binY;
        return binY;
    }

    friend bool operator<(const contactRecord &c1, const contactRecord &c2) {
        if (c1.binY < c2.binY)
            return true;
        if (c1.binY > c2.binY)
            return false;
        if (c1.binY == c2.binY) {
            return c1.binX < c2.binX;
        }
        return true;
    }

    friend bool operator==(const contactRecord &c1, const contactRecord &c2) {
        if (c1.binY != c2.binY)
            return false;
        if (c1.binX != c2.binX)
            return false;
        return true;
    }
};

inline void ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int ch) {
            return !std::isspace(ch);
        }));
    }

// trim from end (in place)
inline void rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), [](int ch) {
            return !std::isspace(ch);
        }).base(), s.end());
    }

// trim from both ends (in place)
inline void trim(std::string &s) {
        ltrim(s);
        rtrim(s);
    }

    inline void toLower(std::string &s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {return std::tolower(c);});
}

inline void replace(std::string &s, const std::string &from ,const std::string &to) {
    size_t start_pos = s.find(from);
    if (start_pos != std::string::npos) s.replace(start_pos, from.length(), to);
}

} // namespace FreeHiC

#endif // HICLIB_COMMON_HPP
