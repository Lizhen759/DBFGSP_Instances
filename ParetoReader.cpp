#include <algorithm>
#include "ParetoReader.h"
//
//// 读取Pareto解集文件的函数
//void ParetoReader::readParetoSetFromFile(const std::string& filename, std::vector<std::vector<Individual>>& paretoSets)
//{
//    std::ifstream file(filename);
//
//    if (!file.is_open()) {
//        throw std::ios_base::failure("Failed to open file.");
//    }
//
//    std::string line;
//    int instanceCount = 0;  // 用于计数实例
//
//    // 逐行读取文件
//    while (std::getline(file, line)) {
//        std::vector<Individual> currentParetoSet;
//        std::replace(line.begin(), line.end(), ',', ' ');
//
//        std::istringstream stream(line);
//        double ms, tec;
//
//        // 读取 MS 和 TEC
//        while (stream >> ms >> tec) {
//            Individual sol;
//            sol.MS = ms;
//            sol.TEC = tec;
//            currentParetoSet.push_back(sol);
//        }
//
//        // 只有非空解集才添加
//        if (!currentParetoSet.empty()) {
//            paretoSets.push_back(currentParetoSet);
//            instanceCount++;
//        }
//    }
//
//    file.close();
//
//    // 输出读取的解集数量
//    if (instanceCount > 0) {
//       // std::cout << " Successfully read " << instanceCount << " Pareto sets from " << filename << std::endl;
//    } else {
//       // std::cerr << "️ Warning: No valid Pareto sets found in " << filename << std::endl;
//    }
//}
void ParetoReader::readParetoSetFromFile(const std::string& filename, std::vector<std::vector<Individual>>& paretoSets)
{
    paretoSets.clear(); // ✅ 每次读取前先清空

    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::ios_base::failure("Failed to open file.");
    }

    std::string line;
    int instanceCount = 0;

    while (std::getline(file, line)) {
        std::vector<Individual> currentParetoSet;
        std::replace(line.begin(), line.end(), ',', ' ');
        std::istringstream stream(line);
        double ms, tec;

        while (stream >> ms >> tec) {
            Individual sol;
            sol.MS = ms;
            sol.TEC = tec;
            currentParetoSet.push_back(sol);
        }

        if (!currentParetoSet.empty()) {
            paretoSets.push_back(currentParetoSet);
            instanceCount++;
        }
    }

    file.close();
}

