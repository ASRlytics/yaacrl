#include <vector>
#include <tuple>
#include <array>

#include "spectrogram.h"


typedef std::tuple<int, int> Peak;  // <window, bin>
typedef std::array<char, 16> Hash;

std::vector<Peak> find_peaks(Spectrogram spec);
std::vector<Hash*> generate_hashes(std::vector<Peak> peaks);
