#define BOOST_TEST_MODULE vectorTests
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <vector>
#include <vectorOps.h>

BOOST_AUTO_TEST_CASE(checkCompareVectors) {

    // Define vectors
    std::vector<double> vec1 = {1.0, 2.0, 3.0};
    std::vector<double> vec2(vec1);
    std::vector<double> vec3 = {-1.0,-2.0,-3.0};
    std::vector<double> vec4 = {0.0, 5.0, 9.99};

    // Equal
    BOOST_REQUIRE_EQUAL(1, compareVectors(vec1,vec2));
    // Flipped
    BOOST_REQUIRE_EQUAL(2, compareVectors(vec1,vec3));
    // Not equal
    BOOST_REQUIRE_EQUAL(0, compareVectors(vec1,vec4));

}
