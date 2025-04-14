#ifndef DBFGSP_NEW_RAND_H
#define DBFGSP_NEW_RAND_H
#include <random>
#include <cfloat>
#include <type_traits>  // 引入 std::is_integral
using namespace std;

// 返回一个种子随机数生成器
inline std::mt19937 &rand_generator() {
    static thread_local std::mt19937 gen(std::random_device{}());
    return gen;
}
// 生成[min, max]范围内的整数
template<typename T, std::enable_if_t<std::is_integral<T>::value> * = nullptr>
T wyt_rand(T min, T max) {
    std::uniform_int_distribution<T> dist(min, max);
    return dist(rand_generator());
}

// 生成[0, max-1]范围内的整数
template<typename T, std::enable_if_t<std::is_integral<T>::value> * = nullptr>
T wyt_rand(T max) {
    std::uniform_int_distribution<T> dist(0, max - 1);
    return dist(rand_generator());
}

// 生成[min, max)范围内的浮点数
template<typename T, std::enable_if_t<std::is_floating_point<T>::value> * = nullptr>
T wyt_rand(T min, T max) {
    std::uniform_real_distribution<T> dist(min, max);
    return dist(rand_generator());
}

// 生成[min, max]范围内的浮点数（包括右端点）
template<typename T, std::enable_if_t<std::is_floating_point<T>::value> * = nullptr>
T wyt_rand_include_right(T min, T max) {
    std::uniform_real_distribution<T> dist(min, std::nextafter(max, DBL_MAX));
    return dist(rand_generator());
}

// 生成布尔值
bool wyt_rand(double par = 0.5);

#endif //DBFGSP_NEW_RAND_H
