#include <vector>
#include <tuple>
#include <algorithm>
#include <iostream>
#include <set>

#include "spectrogram.h"
#include "fingerprint.h"

#include "vendor/MurmurHash3.h"


#define LOC_BINS 40
#define LOC_WINDOWS 20
#define MIN_AMPLITUDE -10
#define MAX_FAN 25


bool is_local_maximum(Spectrogram& spec, int window, int bin) {
    float value = spec[window][bin];

    for (int w_shift = -LOC_WINDOWS; w_shift < LOC_WINDOWS; w_shift++) {
        if (window + w_shift < 0 || window + w_shift >= spec.size()) 
        { continue; }

        for (int b_shift = -LOC_BINS; b_shift < LOC_BINS; b_shift++) {
            if (bin + b_shift < 0 || bin + b_shift >= BINS_AMOUNT) 
            { continue; }

            if (w_shift == 0 && b_shift == 0)
            { continue; }

            if (spec[window+w_shift][bin+b_shift] >= value) 
            { return false; }
        }
    }

    return true;
}

std::vector<Peak> find_peaks(Spectrogram& spec) {
    std::set<int> bitmap_by_windows;
    std::vector<int> bitmap_by_bins;
    
    int cur_max_i;
    int cur_max;
    bool pushed;
    // std::cout << "12312312" << std::endl;
    // Iterate all windows of single bin
    for (int bin = 0; bin < BINS_AMOUNT; bin++) {
        cur_max_i = 0;
        cur_max = spec[0][bin];
        pushed = false;
        if (cur_max >= MIN_AMPLITUDE) {
            pushed = true;
            bitmap_by_windows.insert(bin);
        }
        for (int window = 0; window < spec.size(); window++) {
            auto& value = spec[window][bin];
            if ((window - cur_max_i) >= LOC_WINDOWS) {
                cur_max_i = window;
                cur_max = value;
                if (cur_max >= MIN_AMPLITUDE) {
                    pushed = true;
                    bitmap_by_windows.insert(window*BINS_AMOUNT + bin);
                }
                continue;
            }
            if (value < MIN_AMPLITUDE) continue;

            if (cur_max == value) {
                // Uncheck old max, but do not check new max
                if (pushed) {
                    bitmap_by_windows.erase(cur_max_i*BINS_AMOUNT + bin);
                    pushed = false;
                }
                cur_max_i = window;
            }
            if (cur_max < value) {
                if (pushed) {
                    bitmap_by_windows.erase(cur_max_i*BINS_AMOUNT + bin);
                }
                cur_max_i = window;
                cur_max = value;
                bitmap_by_windows.insert(window*BINS_AMOUNT + bin);
                pushed = true;
            }
        }
    }

    // Iterate all bins of a single window of single bin
    for (int window = 0; window < spec.size(); window++) {
        cur_max_i = 0;
        cur_max = spec[window][0];
        pushed = false;
        if (cur_max >= MIN_AMPLITUDE) {
            pushed = true;
            bitmap_by_bins.push_back(window*BINS_AMOUNT);
        }
        for (int bin = 0; bin < BINS_AMOUNT; bin++) {
            auto& value = spec[window][bin];
            if ((bin - cur_max_i) >= LOC_BINS) {
                cur_max_i = bin;
                cur_max = value;
                if (cur_max >= MIN_AMPLITUDE) {
                    pushed = true;
                    bitmap_by_bins.push_back(window*BINS_AMOUNT + bin);
                }
                continue;
            }
            if (value < MIN_AMPLITUDE) continue;

            if (cur_max == value) {
                // Uncheck old max, but do not check new max
                if (pushed) {
                    bitmap_by_bins.pop_back();
                    pushed = false;
                }
                cur_max_i = bin;
            }
            if (cur_max < value) {
                if (pushed) {
                    bitmap_by_bins.pop_back();
                }
                cur_max_i = bin;
                cur_max = value;
                bitmap_by_bins.push_back(window*BINS_AMOUNT + bin);
                pushed = true;
            }
        }
    }

    std::vector<Peak> peaks;
    // std::cout << bitmap_by_windows.size() << " AND set of " << bitmap_by_bins.size() << std::endl;

    // bitmap2 is already sorted
    for (auto& v: bitmap_by_bins) {
        if (bitmap_by_windows.find(v) != bitmap_by_windows.end()) {
            auto window = v % BINS_AMOUNT;
            auto bin = v - (window * BINS_AMOUNT);
            peaks.emplace_back(window, bin);
        }
    }

    // std::cout << "asdasdasd  " << peaks.size() << std::endl;
    return peaks;

    
}

std::vector<Hash> generate_hashes(std::vector<Peak>& peaks) {
    std::vector<Hash> hashes;

    auto format_key = new int[3];
    for(int i = 0; i < peaks.size(); i++) {
        int c = 0;
        for (int j = i + 1; (j < peaks.size() && c < MAX_FAN); j++) {
            auto distance = std::get<0>(peaks[j]) - std::get<0>(peaks[i]);
            if (distance == 0) continue;

            format_key[0] = std::get<1>(peaks[i]);  // freq 1
            format_key[1] = std::get<1>(peaks[j]);  // freq 2
            format_key[2] = distance;  // window2 - window1
            
            Hash new_hash;
            MurmurHash3_x64_128(
                format_key,
                sizeof(int) * 3,
                0,
                std::get<0>(new_hash).data()
            );

            std::get<1>(new_hash) = std::get<0>(peaks[i]);

            hashes.push_back(new_hash);
            c++;
        }
    }
    
    delete format_key;
    return hashes;
}