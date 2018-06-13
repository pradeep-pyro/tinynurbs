#include "util.h"

namespace tinynurbs {
namespace util {

unsigned int binomial(unsigned int n, unsigned int k) {
    unsigned int result = 1;
    if (k > n) {
        return 0;
    }
    for (unsigned int i = 1; i <= k; ++i) {
        result *= (n + 1 - i);
        result /= i;
    }
    return result;
}

} // namespace util
} // namespace tinynurbs
