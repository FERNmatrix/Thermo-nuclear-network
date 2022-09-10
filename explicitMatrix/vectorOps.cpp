#include <vector>
#include <algorithm>
#include <functional>

/**
 * Compare two vectors
 * @param rv1 First vector to compare
 * @param rv2 Second vector to compare
 * @return 0 if the vectors are not equal, 1 if they are the same, 2 if one
 * vector is the negative of the other
 */
int compareVectors(const std::vector<double> & rv1, const std::vector<double> & rv2) {
    
    int retVal = 0;
    // Check for element-wise equality
    auto equal = std::equal(rv1.begin(),rv1.end(),rv2.begin(),rv2.end());
    // If not equal, are they negated images?
    if (!equal) {
       // Flip rv2 by using the handy transformation instead of a for loop
       std::vector<double> flippedRv2(rv2);
       std::transform(rv2.begin(), rv2.end(), flippedRv2.begin(),
               std::bind(std::multiplies<double>(), std::placeholders::_1, -1.0));
       // Check the value and set the flag
       bool flipped = std::equal(rv1.begin(), rv1.end(), flippedRv2.begin(), flippedRv2.end());
       if (flipped) retVal = 2;
    } else {
       // Just flip the flag
       retVal = 1;
    }

    return retVal;
}
