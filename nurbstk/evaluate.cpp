#include "evaluate.h"

namespace nurbstk {

bool isValidRelation(unsigned int degree, size_t nKnots, size_t nCtrlPts) {
    return nKnots - degree - 1 == nCtrlPts;
}

}
