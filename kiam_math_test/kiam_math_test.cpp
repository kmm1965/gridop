#define _USE_MATH_DEFINES
//#include <kiam/math/array_value.hpp>
#include <kiam/math/vector2_value.hpp>
//#include <kiam/math/dual3.hpp>
using namespace _KIAM_MATH;

//int main()
//{
//    array_value<double, 3> av1, av2, av3;
//    av1 = 1.;
//    av2 = 2.;
//    av1 += av2;
//    av1 -= av2;
//    av1 *= 2;
//    av1 /= 2;
//    const double l = av1.length();
//    av3 = av1 + av2;
//    av3 = av1 - av2;
//    av2 = av1 * 2.;
//    av2 = 2. * av1;
//    av2 = av1 / 2.;
//    const double sc = av1 & av2;
//    av3 = av1 ^ av2;
//    math_copy(av1.cbegin(), av1.cend(), av3.begin());
//}

int main()
{
    vector2_value<double> v1, v2, v3;
    v3 = v1 + v2;
    v3 = v1 - v2;
    //v3 = v1 * v2;
    const double sc = v1 & v2;
    v3 = v1 ^ v2;
    v3 = v1 * 2.;
    v3 = 2. * v1;
    v3 = v1 / 2.;
    math_copy(v1.cbegin(), v1.cend(), v3.begin());
}

//int main()
//{
//    dual3<double> d1, d2, d3;
//    d3 = d1 + d2;
//    d3 = d1 - d2;
//    d3 = d1 * d2;
//    d3 = d1 * 2.;
//    d3 = 2. * d1;
//    d3 = d1 / 2.;
//    d3 = func::sin(d1);
//    d3 = func::cos(d1);
//}
