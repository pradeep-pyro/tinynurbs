#include "evaluate.h"

namespace nurbstk {

bool isValidRelation(unsigned int degree, size_t num_knots, size_t num_ctrl_pts) {
    return (num_knots - degree - 1) == num_ctrl_pts;
}

}
