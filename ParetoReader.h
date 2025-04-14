#ifndef PARETOREADER_H
#define PARETOREADER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include "Individual.h"

class ParetoReader {
public:
    // 构造函数
    ParetoReader() = default;

    // 读取Pareto解集文件的函数
    void readParetoSetFromFile(const std::string& filename, std::vector<std::vector<Individual>>& paretoSets);
};

#endif // PARETOREADER_H
