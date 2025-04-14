#include "ParetoOptimizer.h"
#include <iostream>
#include <algorithm>
#include <valarray>
#include "Individual.h"

// 合并多个Pareto集合
void ParetoOptimizer::mergeParetoSets(int ins,
                                      vector<vector<Individual>>& TempTotalParetoSet,  // 非 const 引用：我们要修改这个容器
                                      const vector<vector<Individual>>& MOSAParetoSet, // const 引用：不会修改这个集合
                                      const vector<vector<Individual>>& CMOEAParetoSet,
                                      const vector<vector<Individual>>& GCCEAParetoSet,
                                      const vector<vector<Individual>>& TSCCEAParetoSet,
                                      const vector<vector<Individual>> &QCCEAParetoSet)
{
    // 合并所有Pareto集合
    for (int PS = 0; PS < MOSAParetoSet[ins].size(); PS++) {
        TempTotalParetoSet[ins].push_back(MOSAParetoSet[ins][PS]);
    }
    for (int PS = 0; PS < CMOEAParetoSet[ins].size(); PS++) {
        TempTotalParetoSet[ins].push_back(CMOEAParetoSet[ins][PS]);
    }
    for (int PS = 0; PS < GCCEAParetoSet[ins].size(); PS++) {
        TempTotalParetoSet[ins].push_back(GCCEAParetoSet[ins][PS]);
    }
    for (int PS = 0; PS < TSCCEAParetoSet[ins].size(); PS++) {
        TempTotalParetoSet[ins].push_back(TSCCEAParetoSet[ins][PS]);
    }
    for (int PS = 0; PS < QCCEAParetoSet[ins].size(); PS++) {
        TempTotalParetoSet[ins].push_back(QCCEAParetoSet[ins][PS]);
    }
}

// 合并多个Pareto集合
void ParetoOptimizer::mergeParetoSets(int ins,
                                      vector<vector<Individual>>& TempTotalParetoSet,  // 非 const 引用：我们要修改这个容器
                                      const vector<vector<Individual>>& QCCEAParetoSet, // const 引用：不会修改这个集合
                                      const vector<vector<Individual>>& QCCEANoQParetoSet,
                                      const vector<vector<Individual>>& QCCEANoRapidParetoSet)
{
    TempTotalParetoSet[ins].clear();
    // 合并所有Pareto集合
    TempTotalParetoSet[ins].insert(TempTotalParetoSet[ins].end(),QCCEANoQParetoSet[ins].begin(), QCCEANoQParetoSet[ins].end());
    TempTotalParetoSet[ins].insert(TempTotalParetoSet[ins].end(),QCCEANoRapidParetoSet[ins].begin(), QCCEANoRapidParetoSet[ins].end());
    TempTotalParetoSet[ins].insert(TempTotalParetoSet[ins].end(),QCCEAParetoSet[ins].begin(), QCCEAParetoSet[ins].end());

}

// 合并多个Pareto集合
void ParetoOptimizer::mergeParetoSets(int ins,
                                      vector<vector<Individual>>& TempTotalParetoSet,  // 非 const 引用：我们要修改这个容器
                                      const vector<vector<Individual>>& QCCEAParetoSet, // const 引用：不会修改这个集合
                                      const vector<vector<Individual>>& MILPParetoSet)
{
    TempTotalParetoSet[ins].clear();
    // 合并所有Pareto集合
    TempTotalParetoSet[ins].insert(TempTotalParetoSet[ins].end(),MILPParetoSet[ins].begin(), MILPParetoSet[ins].end());
    TempTotalParetoSet[ins].insert(TempTotalParetoSet[ins].end(),QCCEAParetoSet[ins].begin(), QCCEAParetoSet[ins].end());

}

void ParetoOptimizer::deemphasizeDominatedSolutions(int ins,vector<vector<Individual>>& TempTotalParetoSet,
                                                    vector<vector<Individual>>& TotalParetoSet)
{
    // 弱化支配解
    for (int j = 0; j < TempTotalParetoSet[ins].size(); j++) {
        for (int i = 0; i < TempTotalParetoSet[ins].size(); i++) {
            if (TempTotalParetoSet[ins][i].flag == 0) {
                if (TempTotalParetoSet[ins][i].pareto_rel[j] == 1) {
                    TempTotalParetoSet[ins][i].flag = 999;
                }
            }
        }
    }

    // 将非支配解添加到 TotalParetoSet 中
    for (int i = 0; i < TempTotalParetoSet[ins].size(); i++) {
        if (TempTotalParetoSet[ins][i].flag == 0) {
            TotalParetoSet[ins].push_back(TempTotalParetoSet[ins][i]);
        }
    }
}

//Pareto_relation
void ParetoOptimizer::mainPareto_relation(vector<Individual>& Population)
{
    int pflag, pflag1;
    for (int i = 0; i < Population.size(); i++)
    {
        Population[i].flag = 0;
        for (int j = 0; j < Population.size(); j++)
        {
            pflag = 1;
            pflag1 = 0;
            if (Population[i].MS > Population[j].MS)
                pflag = 0;
            if (Population[i].MS < Population[j].MS)
                pflag1 = 1;
            if (Population[i].TEC > Population[j].TEC)
                pflag = 0;
            if (Population[i].TEC < Population[j].TEC)
                pflag1 = 1;


            if (pflag == 1 && pflag1 == 1)
                Population[j].pareto_rel[i] = 1;	//i占优j
            else
                Population[j].pareto_rel[i] = 0;
        }

    }
}

void ParetoOptimizer::removeDuplicateSolutions(int ins,vector<vector<Individual>>& TotalParetoSet,vector<vector<Individual>>& FinalTotalParetoSet)
{
    for (int i = 0; i < TotalParetoSet[ins].size(); i++) {
        bool fg = true;
        // 检查是否有相同的解（比较 MS 和 TEC）
        for (int j = 0; j < FinalTotalParetoSet[ins].size(); j++) {
            if ((FinalTotalParetoSet[ins][j].MS == TotalParetoSet[ins][i].MS) &&
                (FinalTotalParetoSet[ins][j].TEC == TotalParetoSet[ins][i].TEC)) {
                fg = false;
                break;
            }
        }
        // 如果没有重复解，添加到 FinalTotalParetoSet 中
        if (fg) {
            FinalTotalParetoSet[ins].push_back(TotalParetoSet[ins][i]);
        }
    }
}

void ParetoOptimizer::mainNormalize(vector<Individual>& totalPopulation,
                                 vector<Individual>& MOSAPopulation,
                                 vector<Individual>& CMOEAPopulation,
                                 vector<Individual>& GCCEAPopulation,
                                 vector<Individual>& TSCCPopulation,
                                  vector<Individual>& QCCEAPopulation) {
    // 初始化 min/max
    int minMS = totalPopulation[0].MS;
    float minTEC = totalPopulation[0].TEC;
    int maxMS = totalPopulation[0].MS;
    float maxTEC = totalPopulation[0].TEC;

    // 统一更新 min/max
    auto updateBounds = [&minMS, &maxMS, &minTEC, &maxTEC](const vector<Individual>& population) {
        for (const auto& ind : population) {
            if (ind.MS < minMS) minMS = ind.MS;
            if (ind.MS > maxMS) maxMS = ind.MS;
            if (ind.TEC < minTEC) minTEC = ind.TEC;
            if (ind.TEC > maxTEC) maxTEC = ind.TEC;
        }
    };

    updateBounds(totalPopulation);
    updateBounds(MOSAPopulation);
    updateBounds(CMOEAPopulation);
    updateBounds(GCCEAPopulation);
    updateBounds(TSCCPopulation);
    updateBounds(QCCEAPopulation);

    // 统一归一化
    auto normalizePopulation = [minMS, maxMS, minTEC, maxTEC](vector<Individual>& population) {
        for (auto& ind : population) {
            ind.normalMS = (ind.MS - static_cast<float>(minMS)) / (maxMS - minMS);
            ind.normalTEC = (ind.TEC - minTEC) / (maxTEC - minTEC);
        }
    };

    normalizePopulation(totalPopulation);
    normalizePopulation(MOSAPopulation);
    normalizePopulation(CMOEAPopulation);
    normalizePopulation(GCCEAPopulation);
    normalizePopulation(TSCCPopulation);
    normalizePopulation(QCCEAPopulation);
}


void ParetoOptimizer::mainNormalize(vector<Individual>& totalPopulation,
                                    vector<Individual>& QCCEAPopulation,
                                    vector<Individual>& QCCEANoQParetoSet,
                                    vector<Individual>& QCCEANoRapidParetoSet) {
    // 初始化 min/max
    int minMS = totalPopulation[0].MS;
    float minTEC = totalPopulation[0].TEC;
    int maxMS = totalPopulation[0].MS;
    float maxTEC = totalPopulation[0].TEC;

    // 统一更新 min/max
    auto updateBounds = [&minMS, &maxMS, &minTEC, &maxTEC](const vector<Individual>& population) {
        for (const auto& ind : population) {
            if (ind.MS < minMS) minMS = ind.MS;
            if (ind.MS > maxMS) maxMS = ind.MS;
            if (ind.TEC < minTEC) minTEC = ind.TEC;
            if (ind.TEC > maxTEC) maxTEC = ind.TEC;
        }
    };

    updateBounds(totalPopulation);
    updateBounds(QCCEAPopulation);
    updateBounds(QCCEANoQParetoSet);
    updateBounds(QCCEANoRapidParetoSet);


    // 统一归一化
    auto normalizePopulation = [minMS, maxMS, minTEC, maxTEC](vector<Individual>& population) {
        for (auto& ind : population) {
            ind.normalMS = (ind.MS - static_cast<float>(minMS)) / (maxMS - minMS);
            ind.normalTEC = (ind.TEC - minTEC) / (maxTEC - minTEC);
        }
    };

    normalizePopulation(totalPopulation);
    normalizePopulation(QCCEAPopulation);
    normalizePopulation(QCCEANoQParetoSet);
    normalizePopulation(QCCEANoRapidParetoSet);

}



void ParetoOptimizer::mainNormalize(vector<Individual>& totalPopulation,
                                    vector<Individual>& QCCEAPopulation,
                                    vector<Individual>& MILPParetoSet) {
    // 初始化 min/max
    int minMS = totalPopulation[0].MS;
    float minTEC = totalPopulation[0].TEC;
    int maxMS = totalPopulation[0].MS;
    float maxTEC = totalPopulation[0].TEC;

    // 统一更新 min/max
    auto updateBounds = [&minMS, &maxMS, &minTEC, &maxTEC](const vector<Individual>& population) {
        for (const auto& ind : population) {
            if (ind.MS < minMS) minMS = ind.MS;
            if (ind.MS > maxMS) maxMS = ind.MS;
            if (ind.TEC < minTEC) minTEC = ind.TEC;
            if (ind.TEC > maxTEC) maxTEC = ind.TEC;
        }
    };

    updateBounds(totalPopulation);
    updateBounds(QCCEAPopulation);
    updateBounds(MILPParetoSet);


    // 统一归一化
    auto normalizePopulation = [minMS, maxMS, minTEC, maxTEC](vector<Individual>& population) {
        for (auto& ind : population) {
            ind.normalMS = (ind.MS - static_cast<float>(minMS)) / (maxMS - minMS);
            ind.normalTEC = (ind.TEC - minTEC) / (maxTEC - minTEC);
        }
    };

    normalizePopulation(totalPopulation);
    normalizePopulation(QCCEAPopulation);
    normalizePopulation(MILPParetoSet);

}



void ParetoOptimizer::AdPrummainNormalize(vector<Individual>& TotalParetoPopulation,vector<Individual>& ParetoPopulation1,vector<Individual>& ParetoPopulation2,
        vector<Individual>& ParetoPopulation3,vector<Individual>& ParetoPopulation4,
        vector<Individual>& ParetoPopulation5,vector<Individual>& ParetoPopulation6,
        vector<Individual>& ParetoPopulation7,vector<Individual>& ParetoPopulation8,vector<Individual>& ParetoPopulation9) {
    // 1. 初始化 min/max
    int minMS = TotalParetoPopulation[0].MS;
    int maxMS = TotalParetoPopulation[0].MS;
    float minTEC = TotalParetoPopulation[0].TEC;
    float maxTEC = TotalParetoPopulation[0].TEC;

    // 2. 更新 min/max（所有前沿）
    auto updateMinMax = [&](const vector<Individual>& pop) {
        for (const auto& ind : pop) {
            if (ind.MS < minMS) minMS = ind.MS;
            if (ind.MS > maxMS) maxMS = ind.MS;
            if (ind.TEC < minTEC) minTEC = ind.TEC;
            if (ind.TEC > maxTEC) maxTEC = ind.TEC;
        }
    };

    updateMinMax(TotalParetoPopulation);
    updateMinMax(ParetoPopulation1);
    updateMinMax(ParetoPopulation2);
    updateMinMax(ParetoPopulation3);
    updateMinMax(ParetoPopulation4);
    updateMinMax(ParetoPopulation5);
    updateMinMax(ParetoPopulation6);
    updateMinMax(ParetoPopulation7);
    updateMinMax(ParetoPopulation8);
    updateMinMax(ParetoPopulation9);

    // 3. 归一化所有前沿
    auto normalizePopulation = [&](vector<Individual>& pop) {
        for (auto& ind : pop) {
            ind.normalMS = (ind.MS - static_cast<float>(minMS)) / (maxMS - minMS);
            ind.normalTEC = (ind.TEC - minTEC) / (maxTEC - minTEC);
        }
    };

    normalizePopulation(TotalParetoPopulation);
    normalizePopulation(ParetoPopulation1);
    normalizePopulation(ParetoPopulation2);
    normalizePopulation(ParetoPopulation3);
    normalizePopulation(ParetoPopulation4);
    normalizePopulation(ParetoPopulation5);
    normalizePopulation(ParetoPopulation6);
    normalizePopulation(ParetoPopulation7);
    normalizePopulation(ParetoPopulation8);
    normalizePopulation(ParetoPopulation9);
}

//void ParetoOptimizer::AdPrumComputeIGD(int ins, const vector<Individual>& TotalParetoPopulation, const vector<Individual>& ParetoPopulation1, const vector<Individual>& ParetoPopulation2,
//                      const vector<Individual>& ParetoPopulation3, const vector<Individual>& ParetoPopulation4, const vector<Individual>& ParetoPopulation5, const vector<Individual>& ParetoPopulation6,
//                      const vector<Individual>& ParetoPopulation7, const vector<Individual>& ParetoPopulation8, const vector<Individual>& ParetoPopulation9,
//                      double& IGD1, double& IGD2, double& IGD3, double& IGD4, double& IGD5, double& IGD6, double& IGD7,double& IGD8, double& IGD9)
//{
//
//
//    for (int j = 0; j < TotalParetoPopulation.size(); j++)
//    {
//        vector<double> tempValue;
//        tempValue.resize(ParetoPopulation1.size(), 0.0);
//        int index = 0;
//        double minValue = 1000000.0;
//        for (int i = 0; i < ParetoPopulation1.size(); i++)
//        {
//            tempValue[i] += (ParetoPopulation1[i].normalMS - TotalParetoPopulation[j].normalMS) * (ParetoPopulation1[i].normalMS - TotalParetoPopulation[j].normalMS);
//           // cout << tempValue[i] << endl;
//            tempValue[i] += (ParetoPopulation1[i].normalTEC - TotalParetoPopulation[j].normalTEC) * (ParetoPopulation1[i].normalTEC - TotalParetoPopulation[j].normalTEC);
//
//            double value = sqrt(tempValue[i]);
//            if (value < minValue)
//            {
//                minValue = value;
//                index = i;
//            }
//        }
//        IGD1 += minValue;
//    }
//    IGD1 = IGD1 / TotalParetoPopulation.size();
//  //  cout << "1反世代距离：" << IGD1 << endl;
//
//    for (int j = 0; j < TotalParetoPopulation.size(); j++)
//    {
//        vector<double> tempValue;
//        tempValue.resize(ParetoPopulation2.size(), 0.0);
//        int index = 0;
//        double minValue = 1000000.0;
//        for (int i = 0; i < ParetoPopulation2.size(); i++)
//        {
//            tempValue[i] += (ParetoPopulation2[i].normalMS - TotalParetoPopulation[j].normalMS) * (ParetoPopulation2[i].normalMS - TotalParetoPopulation[j].normalMS);
//
//            tempValue[i] += (ParetoPopulation2[i].normalTEC - TotalParetoPopulation[j].normalTEC) * (ParetoPopulation2[i].normalTEC - TotalParetoPopulation[j].normalTEC);
//
//            double value = sqrt(tempValue[i]);
//            if (value < minValue)
//            {
//                minValue = value;
//                index = i;
//            }
//        }
//        IGD2 += minValue;
//    }
//
//    IGD2 = IGD2 / TotalParetoPopulation.size();
//  //  cout << "2反世代距离：" << IGD2 << endl;
//
//    for (int j = 0; j < TotalParetoPopulation.size(); j++)
//    {
//        vector<double> tempValue;
//        tempValue.resize(ParetoPopulation3.size(), 0.0);
//        int index = 0;
//        double minValue = 1000000.0;
//        for (int i = 0; i < ParetoPopulation3.size(); i++)
//        {
//            tempValue[i] += (ParetoPopulation3[i].normalMS - TotalParetoPopulation[j].normalMS) * (ParetoPopulation3[i].normalMS - TotalParetoPopulation[j].normalMS);
//
//            tempValue[i] += (ParetoPopulation3[i].normalTEC - TotalParetoPopulation[j].normalTEC) * (ParetoPopulation3[i].normalTEC - TotalParetoPopulation[j].normalTEC);
//
//            double value = sqrt(tempValue[i]);
//            if (value < minValue)
//            {
//                minValue = value;
//                index = i;
//            }
//        }
//        IGD3 += minValue;
//    }
//
//    IGD3 = IGD3 / TotalParetoPopulation.size();
//  //  cout << "3反世代距离：" << IGD3 << endl;
//
//    for (int j = 0; j < TotalParetoPopulation.size(); j++)
//    {
//        vector<double> tempValue;
//        tempValue.resize(ParetoPopulation4.size(), 0.0);
//        int index = 0;
//        double minValue = 1000000.0;
//        for (int i = 0; i < ParetoPopulation4.size(); i++)
//        {
//            tempValue[i] += (ParetoPopulation4[i].normalMS - TotalParetoPopulation[j].normalMS) * (ParetoPopulation4[i].normalMS - TotalParetoPopulation[j].normalMS);
//
//            tempValue[i] += (ParetoPopulation4[i].normalTEC - TotalParetoPopulation[j].normalTEC) * (ParetoPopulation4[i].normalTEC - TotalParetoPopulation[j].normalTEC);
//
//            double value = sqrt(tempValue[i]);
//            if (value < minValue)
//            {
//                minValue = value;
//                index = i;
//            }
//        }
//        IGD4 += minValue;
//    }
//    IGD4 = IGD4 / TotalParetoPopulation.size();
//  //  cout << "4反世代距离：" << IGD4 << endl;
//
//
//    for (int j = 0; j < TotalParetoPopulation.size(); j++)
//    {
//        vector<double> tempValue;
//        tempValue.resize(ParetoPopulation5.size(), 0.0);
//        int index = 0;
//        double minValue = 1000000.0;
//        for (int i = 0; i < ParetoPopulation5.size(); i++)
//        {
//            tempValue[i] += (ParetoPopulation5[i].normalMS - TotalParetoPopulation[j].normalMS) * (ParetoPopulation5[i].normalMS - TotalParetoPopulation[j].normalMS);
//
//            tempValue[i] += (ParetoPopulation5[i].normalTEC - TotalParetoPopulation[j].normalTEC) * (ParetoPopulation5[i].normalTEC - TotalParetoPopulation[j].normalTEC);
//
//            double value = sqrt(tempValue[i]);
//            if (value < minValue)
//            {
//                minValue = value;
//                index = i;
//            }
//        }
//        IGD5 += minValue;
//    }
//    IGD5 = IGD5 / TotalParetoPopulation.size();
//  //  cout << "5反世代距离：" << IGD5 << endl;
//
//    for (int j = 0; j < TotalParetoPopulation.size(); j++)
//    {
//        vector<double> tempValue;
//        tempValue.resize(ParetoPopulation6.size(), 0.0);
//        int index = 0;
//        double minValue = 1000000.0;
//        for (int i = 0; i < ParetoPopulation6.size(); i++)
//        {
//            tempValue[i] += (ParetoPopulation6[i].normalMS - TotalParetoPopulation[j].normalMS) * (ParetoPopulation6[i].normalMS - TotalParetoPopulation[j].normalMS);
//
//            tempValue[i] += (ParetoPopulation6[i].normalTEC - TotalParetoPopulation[j].normalTEC) * (ParetoPopulation6[i].normalTEC - TotalParetoPopulation[j].normalTEC);
//
//            double value = sqrt(tempValue[i]);
//            if (value < minValue)
//            {
//                minValue = value;
//                index = i;
//            }
//        }
//        IGD6 += minValue;
//    }
//
//    IGD6 = IGD6 / TotalParetoPopulation.size();
//   // cout << "6反世代距离：" << IGD6 << endl;
//
//    for (int j = 0; j < TotalParetoPopulation.size(); j++)
//    {
//        vector<double> tempValue;
//        tempValue.resize(ParetoPopulation7.size(), 0.0);
//        int index = 0;
//        double minValue = 1000000.0;
//        for (int i = 0; i < ParetoPopulation7.size(); i++)
//        {
//            tempValue[i] += (ParetoPopulation7[i].normalMS - TotalParetoPopulation[j].normalMS) * (ParetoPopulation7[i].normalMS - TotalParetoPopulation[j].normalMS);
//
//            tempValue[i] += (ParetoPopulation7[i].normalTEC - TotalParetoPopulation[j].normalTEC) * (ParetoPopulation7[i].normalTEC - TotalParetoPopulation[j].normalTEC);
//
//            double value = sqrt(tempValue[i]);
//            if (value < minValue)
//            {
//                minValue = value;
//                index = i;
//            }
//        }
//        IGD7 += minValue;
//    }
//    IGD7 = IGD7 / TotalParetoPopulation.size();
//  //  cout << "7反世代距离：" << IGD7 << endl;
//
//    for (int j = 0; j < TotalParetoPopulation.size(); j++)
//    {
//        vector<double> tempValue;
//        tempValue.resize(ParetoPopulation8.size(), 0.0);
//        int index = 0;
//        double minValue = 1000000.0;
//        for (int i = 0; i < ParetoPopulation8.size(); i++)
//        {
//            tempValue[i] += (ParetoPopulation8[i].normalMS - TotalParetoPopulation[j].normalMS) * (ParetoPopulation8[i].normalMS - TotalParetoPopulation[j].normalMS);
//
//            tempValue[i] += (ParetoPopulation8[i].normalTEC - TotalParetoPopulation[j].normalTEC) * (ParetoPopulation8[i].normalTEC - TotalParetoPopulation[j].normalTEC);
//
//            double value = sqrt(tempValue[i]);
//            if (value < minValue)
//            {
//                minValue = value;
//                index = i;
//            }
//        }
//        IGD8 += minValue;
//    }
//    IGD8 = IGD8 / TotalParetoPopulation.size();
//  //  cout << "8反世代距离：" << IGD8 << endl;
//
//    for (int j = 0; j < TotalParetoPopulation.size(); j++)
//    {
//        vector<double> tempValue;
//        tempValue.resize(ParetoPopulation9.size(), 0.0);
//        int index = 0;
//        double minValue = 1000000.0;
//        for (int i = 0; i < ParetoPopulation9.size(); i++)
//        {
//            tempValue[i] += (ParetoPopulation9[i].normalMS - TotalParetoPopulation[j].normalMS) * (ParetoPopulation9[i].normalMS - TotalParetoPopulation[j].normalMS);
//
//            tempValue[i] += (ParetoPopulation9[i].normalTEC - TotalParetoPopulation[j].normalTEC) * (ParetoPopulation9[i].normalTEC - TotalParetoPopulation[j].normalTEC);
//
//            double value = sqrt(tempValue[i]);
//            if (value < minValue)
//            {
//                minValue = value;
//                index = i;
//            }
//        }
//        IGD9 += minValue;
//    }
//    IGD9 = IGD9 / TotalParetoPopulation.size();
//  //  cout << "9反世代距离：" << IGD9 << endl;
//
//}
// 简化版实现（供内部调用）
void ParetoOptimizer::AdPrumComputeIGD_Internal(
        const vector<Individual>& TotalParetoPopulation,
        const vector<vector<Individual>>& AllParetoPopulations,
        vector<double>& IGDs
) {
    IGDs.assign(AllParetoPopulations.size(), 0.0);

    for (const auto& refInd : TotalParetoPopulation) {
        for (size_t k = 0; k < AllParetoPopulations.size(); ++k) {
            double minDist = numeric_limits<double>::max();
            for (const auto& testInd : AllParetoPopulations[k]) {
                double dist = sqrt(
                        pow(testInd.normalMS - refInd.normalMS, 2) +
                        pow(testInd.normalTEC - refInd.normalTEC, 2)
                );
                minDist = min(minDist, dist);
            }
            IGDs[k] += minDist;
        }
    }

    for (auto& igd : IGDs) {
        igd /= TotalParetoPopulation.size();
    }
}

// 原接口（兼容现有调用）
void ParetoOptimizer::AdPrumComputeIGD(const vector<Individual>& TotalParetoPopulation,
        const vector<Individual>& ParetoPopulation1,const vector<Individual>& ParetoPopulation2,const vector<Individual>& ParetoPopulation3,
        const vector<Individual>& ParetoPopulation4,const vector<Individual>& ParetoPopulation5,
        const vector<Individual>& ParetoPopulation6,const vector<Individual>& ParetoPopulation7,
        const vector<Individual>& ParetoPopulation8,const vector<Individual>& ParetoPopulation9,double& IGD1, double& IGD2, double& IGD3, double& IGD4,
        double& IGD5, double& IGD6, double& IGD7, double& IGD8, double& IGD9
) {
    // 合并输入
    vector<vector<Individual>> AllPopulations = {
            ParetoPopulation1, ParetoPopulation2, ParetoPopulation3,
            ParetoPopulation4, ParetoPopulation5, ParetoPopulation6,
            ParetoPopulation7, ParetoPopulation8, ParetoPopulation9
    };

    // 调用简化版
    vector<double> IGDs(9);
    AdPrumComputeIGD_Internal(TotalParetoPopulation, AllPopulations, IGDs);

    // 分配结果
    IGD1 = IGDs[0]; IGD2 = IGDs[1]; IGD3 = IGDs[2];
    IGD4 = IGDs[3]; IGD5 = IGDs[4]; IGD6 = IGDs[5];
    IGD7 = IGDs[6]; IGD8 = IGDs[7]; IGD9 = IGDs[8];
}

void ParetoOptimizer::ComputeCoverage( const vector<Individual>& TotalParetoPopulation, const vector<Individual>& MOSAParetoPopulation, const vector<Individual>& CMOEAParetoPopulation,
                  const vector<Individual>& GCCEAParetoPopulation, const vector<Individual>& TSCCEAParetoPopulation, const std::vector<Individual>& QCCEAParetoPopulation,
                  double& MOSAC, double& CMOEAC, double& GCCEAC, double& TSCCEAC,double& QCCEAC)
{
    int MOSAcount = 0;
    for (int i = 0; i < MOSAParetoPopulation.size(); i++)
    {
        for (int j = 0; j < TotalParetoPopulation.size(); j++)
        {
            if ((MOSAParetoPopulation[i].MS == TotalParetoPopulation[j].MS) && (MOSAParetoPopulation[i].TEC == TotalParetoPopulation[j].TEC))
            {
                MOSAcount++;
            }
        }
    }
    MOSAC = static_cast<float>(MOSAcount) / TotalParetoPopulation.size();
  //  cout << "MOSA_Coverage：" << MOSAC << endl;

    int CMOEAcount = 0;
    for (int i = 0; i < CMOEAParetoPopulation.size(); i++)
    {
        for (int j = 0; j < TotalParetoPopulation.size(); j++)
        {
            if ((CMOEAParetoPopulation[i].MS == TotalParetoPopulation[j].MS) && (CMOEAParetoPopulation[i].TEC == TotalParetoPopulation[j].TEC))
            {
                CMOEAcount++;
            }
        }
    }
    CMOEAC = static_cast<float>(CMOEAcount) / TotalParetoPopulation.size();
   // cout << "CMOEA_Coverage：" << CMOEAC << endl;

    int GCCEAcount = 0;
    for (int i = 0; i < GCCEAParetoPopulation.size(); i++)
    {
        for (int j = 0; j < TotalParetoPopulation.size(); j++)
        {
            if ((GCCEAParetoPopulation[i].MS == TotalParetoPopulation[j].MS) && (GCCEAParetoPopulation[i].TEC == TotalParetoPopulation[j].TEC))
            {
                GCCEAcount++;
            }
        }
    }
    GCCEAC = static_cast<float>(GCCEAcount) / TotalParetoPopulation.size();
   // cout << "GCCEA_Coverage：" << GCCEAC << endl;

    int TSCCEAcount = 0;
    for (int i = 0; i < TSCCEAParetoPopulation.size(); i++)
    {
        for (int j = 0; j < TotalParetoPopulation.size(); j++)
        {
            if ((TSCCEAParetoPopulation[i].MS == TotalParetoPopulation[j].MS) && (TSCCEAParetoPopulation[i].TEC == TotalParetoPopulation[j].TEC))
            {
                TSCCEAcount++;
            }
        }
    }
    TSCCEAC = static_cast<float>(TSCCEAcount) / TotalParetoPopulation.size();
  //  cout << "TSCCEA_Coverage：" << TSCCEAC << endl;

    int QCCEAcount = 0;
    for (int i = 0; i <  QCCEAParetoPopulation.size(); i++)
    {
        for (int j = 0; j < TotalParetoPopulation.size(); j++)
        {
            if (( QCCEAParetoPopulation[i].MS == TotalParetoPopulation[j].MS) && ( QCCEAParetoPopulation[i].TEC == TotalParetoPopulation[j].TEC))
            {
                QCCEAcount++;
            }
        }
    }
    QCCEAC = static_cast<float>(QCCEAcount) / TotalParetoPopulation.size();
    //cout << " QCCEA_Coverage：" <<  QCCEAC << endl;
}


void ParetoOptimizer::ComputeCoverage_Str( const vector<Individual>& TotalParetoPopulation, const vector<Individual>& QCCEAParetoPopulation,
                                      const vector<Individual>& QCCEANoQParetoPopulation,
                                     const std::vector<Individual>& QCCEANoRapidParetoPopulation,
                                     double& QCCEAC ,double& QCCEANoQC, double& QCCEANoRapidC)
{
    int QCCEAcount = 0;
    for (int i = 0; i <  QCCEAParetoPopulation.size(); i++)
    {
        for (int j = 0; j < TotalParetoPopulation.size(); j++)
        {
            if (( QCCEAParetoPopulation[i].MS == TotalParetoPopulation[j].MS) && ( QCCEAParetoPopulation[i].TEC == TotalParetoPopulation[j].TEC))
            {
                QCCEAcount++;
            }
        }
    }
    QCCEAC = static_cast<float>(QCCEAcount) / TotalParetoPopulation.size();

    int QCCEANoQcount = 0;
    for (int i = 0; i <QCCEANoQParetoPopulation.size(); i++)
    {
        for (int j = 0; j < TotalParetoPopulation.size(); j++)
        {
            if ((QCCEANoQParetoPopulation[i].MS == TotalParetoPopulation[j].MS) && (QCCEANoQParetoPopulation[i].TEC == TotalParetoPopulation[j].TEC))
            {
                QCCEANoQcount++;
            }
        }
    }
   QCCEANoQC = static_cast<float>(QCCEANoQcount) / TotalParetoPopulation.size();


    int QCCEANoRapidcount = 0;
    for (int i = 0; i < QCCEANoRapidParetoPopulation.size(); i++)
    {
        for (int j = 0; j < TotalParetoPopulation.size(); j++)
        {
            if ((QCCEANoRapidParetoPopulation[i].MS == TotalParetoPopulation[j].MS) && (QCCEANoRapidParetoPopulation[i].TEC == TotalParetoPopulation[j].TEC))
            {
                QCCEANoRapidcount++;
            }
        }
    }
   QCCEANoRapidC = static_cast<float>(QCCEANoRapidcount) / TotalParetoPopulation.size();

}



void ParetoOptimizer::ComputeCoverage_MILP( const vector<Individual>& TotalParetoPopulation, const vector<Individual>& QCCEAParetoPopulation,
                                      const vector<Individual>& MILPParetoPopulation,
                                      double& QCCEAC ,double& MILPC)
{
    int QCCEAcount = 0;
    for (int i = 0; i <  QCCEAParetoPopulation.size(); i++)
    {
        for (int j = 0; j < TotalParetoPopulation.size(); j++)
        {
            if (( QCCEAParetoPopulation[i].MS == TotalParetoPopulation[j].MS) && ( QCCEAParetoPopulation[i].TEC == TotalParetoPopulation[j].TEC))
            {
                QCCEAcount++;
            }
        }
    }
    QCCEAC = static_cast<float>(QCCEAcount) / TotalParetoPopulation.size();

    int QCCEANoQcount = 0;
    for (int i = 0; i < MILPParetoPopulation.size(); i++)
    {
        for (int j = 0; j < TotalParetoPopulation.size(); j++)
        {
            if ((MILPParetoPopulation[i].MS == TotalParetoPopulation[j].MS) && (MILPParetoPopulation[i].TEC == TotalParetoPopulation[j].TEC))
            {
                QCCEANoQcount++;
            }
        }
    }
    MILPC = static_cast<float>(QCCEANoQcount) / TotalParetoPopulation.size();

}

// 比较函数，用于排序
bool  ParetoOptimizer::compareSolutions(const Individual& sol1, const Individual& sol2) {
    if (sol1.normalMS == sol2.normalMS)
        return sol1.normalTEC < sol2.normalTEC;
    return sol1.normalMS < sol2.normalMS;
}

void ParetoOptimizer::ComputeHV( const vector<Individual>& TotalParetoPopulation,
                                vector<Individual>& MOSAParetoPopulation,
                                vector<Individual>& CMOEAParetoPopulation,
                                vector<Individual>& GCCEAParetoPopulation,
                                vector<Individual>& TSCCEAParetoPopulation,
                                vector<Individual>& QCCEAParetoPopulation,
                                double& MOSAhv, double& CMOEAhv, double& GCCEAhv, double& TSCCEAhv, double& QCCEAhv)
{
    auto computeHV = [](vector<Individual>& paretoSet) -> double {
        if (paretoSet.empty()) return 0.0;

        // 假设 compareSolutions 是按 normalTEC 升序
        std::sort(paretoSet.begin(), paretoSet.end(), compareSolutions);

        double hv = 0.0;
        double preMS = 1.0;
        double preTEC = 1.0;

        for (const auto& ind : paretoSet) {
            double dMS = preMS - ind.normalMS;
            double dTEC = preTEC - ind.normalTEC;

            if (dMS > 0 && dTEC > 0) {
                hv += dMS * dTEC;
            }

            preTEC = ind.normalTEC;
        }

        return hv;
    };

    MOSAhv     = computeHV(MOSAParetoPopulation);
    CMOEAhv    = computeHV(CMOEAParetoPopulation);
    GCCEAhv    = computeHV(GCCEAParetoPopulation);
    TSCCEAhv   = computeHV(TSCCEAParetoPopulation);
    QCCEAhv    = computeHV(QCCEAParetoPopulation);
}



void ParetoOptimizer::ComputeHV_Str( const vector<Individual>& TotalParetoPopulation,
                                    vector<Individual>& QCCEAParetoPopulation,
                                    vector<Individual>& QCCEANoQParetoPopulation,
                                    vector<Individual>& QCCEANoRapidParetoPopulation,
                                    double& QCCEAHV, double& QCCEANoQHV, double& QCCEANoRapidHV)
{
    auto computeHypervolume = [](vector<Individual>& paretoPopulation) {
        // 按 normalTEC 升序排序（确保从右上往左下积累面积）
        std::sort(paretoPopulation.begin(), paretoPopulation.end(), [](const Individual& a, const Individual& b) {
            return a.normalTEC < b.normalTEC;
        });

        double preMS = 1.0;
        double preTEC = 1.0;
        double HV = 0.0;

        for (const auto& ind : paretoPopulation) {
            double dMS = preMS - ind.normalMS;
            double dTEC = preTEC - ind.normalTEC;

            // 防止负值干扰：确保数值合法（尤其是边界浮点误差）
            if (dMS > 0 && dTEC > 0) {
                HV += dMS * dTEC;
            }

            preTEC = ind.normalTEC;
        }

        return HV;
    };

    QCCEAHV = computeHypervolume(QCCEAParetoPopulation);
    QCCEANoQHV = computeHypervolume(QCCEANoQParetoPopulation);
    QCCEANoRapidHV = computeHypervolume(QCCEANoRapidParetoPopulation);
}


// 计算Hypervolume
void ParetoOptimizer::ComputeHV_MILP(const vector<Individual>& TotalParetoPopulation,
                                vector<Individual>& QCCEAParetoPopulation,
                                vector<Individual>& MILPParetoPopulation,
                                double& QCCEAHV, double& MILPHV)
{
    std::sort(QCCEAParetoPopulation.begin(), QCCEAParetoPopulation.end(), [](const Individual& a, const Individual& b) {
        return a.normalMS < b.normalMS;
    });
    std::sort(MILPParetoPopulation.begin(), MILPParetoPopulation.end(), [](const Individual& a, const Individual& b) {
        return a.normalMS < b.normalMS;
    });

    // 计算Hypervolume的辅助函数
    auto computeHypervolume = [](const vector<Individual>& paretoPopulation) {
        double preMS = 1.0;
        double preTEC = 1.0;
        double HV = 0.0;

        for (const auto& ind : paretoPopulation) {
            HV += (preMS - ind.normalMS) * (preTEC - ind.normalTEC);
            preTEC = ind.normalTEC;
        }

        return HV;
    };

    // 已排序的解集直接计算
    QCCEAHV = computeHypervolume(QCCEAParetoPopulation);
    MILPHV = computeHypervolume(MILPParetoPopulation);
}



//计算IGD
void ParetoOptimizer::ComputeIGD(const vector<Individual>& TotalParetoPopulation, const vector<Individual>& MOSAParetoPopulation, const vector<Individual>& CMOEAParetoPopulation,
                const vector<Individual>& GCCEAParetoPopulation, const vector<Individual>& TSCCEAParetoPopulation,const vector<Individual>& QCCEAParetoPopulation,
                double& MOSAIGD, double& CMOEAIGD, double& GCCEAIGD, double& TSCCEAIGD,double& QCCEAIGD)
{
    for (int j = 0; j < TotalParetoPopulation.size(); j++)
    {
        vector<double> tempValue;
        tempValue.resize(MOSAParetoPopulation.size(), 0.0);
        int index = 0;
        double minValue = 1000000.0;
        for (int i = 0; i < MOSAParetoPopulation.size(); i++)
        {
            tempValue[i] += (MOSAParetoPopulation[i].normalMS - TotalParetoPopulation[j].normalMS) * (MOSAParetoPopulation[i].normalMS - TotalParetoPopulation[j].normalMS);

            tempValue[i] += (MOSAParetoPopulation[i].normalTEC - TotalParetoPopulation[j].normalTEC) * (MOSAParetoPopulation[i].normalTEC - TotalParetoPopulation[j].normalTEC);

            double value = sqrt(tempValue[i]);
            if (value < minValue)
            {
                minValue = value;
                index = i;
            }
        }
        MOSAIGD += minValue;
    }
    MOSAIGD = MOSAIGD / TotalParetoPopulation.size();
  //  cout << "MOSA_IGD：" << MOSAIGD << endl;

    for (int j = 0; j < TotalParetoPopulation.size(); j++)
    {
        vector<double> tempValue;
        tempValue.resize(CMOEAParetoPopulation.size(), 0.0);
        int index = 0;
        double minValue = 1000000.0;
        for (int i = 0; i < CMOEAParetoPopulation.size(); i++)
        {
            tempValue[i] += (CMOEAParetoPopulation[i].normalMS - TotalParetoPopulation[j].normalMS) * (CMOEAParetoPopulation[i].normalMS - TotalParetoPopulation[j].normalMS);

            tempValue[i] += (CMOEAParetoPopulation[i].normalTEC - TotalParetoPopulation[j].normalTEC) * (CMOEAParetoPopulation[i].normalTEC - TotalParetoPopulation[j].normalTEC);

            double value = sqrt(tempValue[i]);
            if (value < minValue)
            {
                minValue = value;
                index = i;
            }
        }
        CMOEAIGD += minValue;
    }
    CMOEAIGD = CMOEAIGD / TotalParetoPopulation.size();
   // cout << "CMOEA_IGD：" << CMOEAIGD << endl;


    for (int j = 0; j < TotalParetoPopulation.size(); j++)
    {
        vector<double> tempValue;
        tempValue.resize(GCCEAParetoPopulation.size(), 0.0);
        int index = 0;
        double minValue = 1000000.0;
        for (int i = 0; i < GCCEAParetoPopulation.size(); i++)
        {
            tempValue[i] += (GCCEAParetoPopulation[i].normalMS - TotalParetoPopulation[j].normalMS) * (GCCEAParetoPopulation[i].normalMS - TotalParetoPopulation[j].normalMS);

            tempValue[i] += (GCCEAParetoPopulation[i].normalTEC - TotalParetoPopulation[j].normalTEC) * (GCCEAParetoPopulation[i].normalTEC - TotalParetoPopulation[j].normalTEC);

            double value = sqrt(tempValue[i]);
            if (value < minValue)
            {
                minValue = value;
                index = i;
            }
        }
        GCCEAIGD += minValue;
    }
    GCCEAIGD = GCCEAIGD / TotalParetoPopulation.size();
   // cout << "GCCEA_IGD：" << GCCEAIGD << endl;

    for (int j = 0; j < TotalParetoPopulation.size(); j++)
    {
        vector<double> tempValue;
        tempValue.resize(TSCCEAParetoPopulation.size(), 0.0);
        int index = 0;
        double minValue = 1000000.0;
        for (int i = 0; i < TSCCEAParetoPopulation.size(); i++)
        {
            tempValue[i] += (TSCCEAParetoPopulation[i].normalMS - TotalParetoPopulation[j].normalMS) * (TSCCEAParetoPopulation[i].normalMS - TotalParetoPopulation[j].normalMS);

            tempValue[i] += (TSCCEAParetoPopulation[i].normalTEC - TotalParetoPopulation[j].normalTEC) * (TSCCEAParetoPopulation[i].normalTEC - TotalParetoPopulation[j].normalTEC);

            double value = sqrt(tempValue[i]);
            if (value < minValue)
            {
                minValue = value;
                index = i;
            }
        }
        TSCCEAIGD += minValue;
    }
    TSCCEAIGD = TSCCEAIGD / TotalParetoPopulation.size();
   // cout << "TSCCEA_IGD：" << TSCCEAIGD << endl;


    for (int j = 0; j < TotalParetoPopulation.size(); j++)
    {
        vector<double> tempValue;
        tempValue.resize(QCCEAParetoPopulation.size(), 0.0);
        int index = 0;
        double minValue = 1000000.0;
        for (int i = 0; i < QCCEAParetoPopulation.size(); i++)
        {
            tempValue[i] += (QCCEAParetoPopulation[i].normalMS - TotalParetoPopulation[j].normalMS) * (QCCEAParetoPopulation[i].normalMS - TotalParetoPopulation[j].normalMS);

            tempValue[i] += (QCCEAParetoPopulation[i].normalTEC - TotalParetoPopulation[j].normalTEC) * (QCCEAParetoPopulation[i].normalTEC - TotalParetoPopulation[j].normalTEC);

            double value = sqrt(tempValue[i]);
            if (value < minValue)
            {
                minValue = value;
                index = i;
            }
        }
        QCCEAIGD += minValue;
    }
    QCCEAIGD = QCCEAIGD / TotalParetoPopulation.size();
   // cout << "QCCEA_IGD：" << QCCEAIGD << endl;

}


void ParetoOptimizer::ComputeIGD_Str(const vector<Individual>& TotalParetoPopulation,
                                 const vector<Individual>& QCCEAParetoPopulation,
                                 const vector<Individual>& QCCEANoQParetoPopulation,
                                 const vector<Individual>& QCCEANoRapidParetoPopulation,
                                 double& QCCEAIGD, double& QCCEANoQIGD, double& QCCEANoRapidIGD)
{
    // 计算 QCCEA IGD
    QCCEAIGD = 0.0;
    for (int j = 0; j < TotalParetoPopulation.size(); j++)
    {
        double minValue = std::numeric_limits<double>::max();  // 使用更适合的值
        for (int i = 0; i < QCCEAParetoPopulation.size(); i++)
        {
            double diffMS = QCCEAParetoPopulation[i].normalMS - TotalParetoPopulation[j].normalMS;
            double diffTEC = QCCEAParetoPopulation[i].normalTEC - TotalParetoPopulation[j].normalTEC;
            double value = sqrt(diffMS * diffMS + diffTEC * diffTEC);
            minValue = std::min(minValue, value);  // 更简洁的最小值计算
        }
        QCCEAIGD += minValue;
    }
    QCCEAIGD /= TotalParetoPopulation.size();  // 计算平均IGD

    // 计算 QCCEANoQ IGD
    QCCEANoQIGD = 0.0;
    for (int j = 0; j < TotalParetoPopulation.size(); j++)
    {
        double minValue = std::numeric_limits<double>::max();
        for (int i = 0; i < QCCEANoQParetoPopulation.size(); i++)
        {
            double diffMS = QCCEANoQParetoPopulation[i].normalMS - TotalParetoPopulation[j].normalMS;
            double diffTEC = QCCEANoQParetoPopulation[i].normalTEC - TotalParetoPopulation[j].normalTEC;
            double value = sqrt(diffMS * diffMS + diffTEC * diffTEC);
            minValue = std::min(minValue, value);
        }
        QCCEANoQIGD += minValue;
    }
    QCCEANoQIGD /= TotalParetoPopulation.size();

    // 计算 QCCEANoRapid IGD
    QCCEANoRapidIGD = 0.0;
    for (int j = 0; j < TotalParetoPopulation.size(); j++)
    {
        double minValue = std::numeric_limits<double>::max();
        for (int i = 0; i < QCCEANoRapidParetoPopulation.size(); i++)
        {
            double diffMS = QCCEANoRapidParetoPopulation[i].normalMS - TotalParetoPopulation[j].normalMS;
            double diffTEC = QCCEANoRapidParetoPopulation[i].normalTEC - TotalParetoPopulation[j].normalTEC;
            double value = sqrt(diffMS * diffMS + diffTEC * diffTEC);
            minValue = std::min(minValue, value);
        }
        QCCEANoRapidIGD += minValue;
    }
    QCCEANoRapidIGD /= TotalParetoPopulation.size();
}

void ParetoOptimizer::ComputeIGD_MILP( const vector<Individual>& TotalParetoPopulation,
                                 const vector<Individual>& QCCEAParetoPopulation,
                                 const vector<Individual>& MILPParetoPopulation,
                                 double& QCCEAIGD, double& MILPIGD)
{
    // 计算 QCCEA IGD
    QCCEAIGD = 0.0;
    for (int j = 0; j < TotalParetoPopulation.size(); j++)
    {
        double minValue = std::numeric_limits<double>::max();  // 使用更适合的值
        for (int i = 0; i < QCCEAParetoPopulation.size(); i++)
        {
            double diffMS = QCCEAParetoPopulation[i].normalMS - TotalParetoPopulation[j].normalMS;
            double diffTEC = QCCEAParetoPopulation[i].normalTEC - TotalParetoPopulation[j].normalTEC;
            double value = sqrt(diffMS * diffMS + diffTEC * diffTEC);
            minValue = std::min(minValue, value);  // 更简洁的最小值计算
        }
        QCCEAIGD += minValue;
    }
    QCCEAIGD /= TotalParetoPopulation.size();  // 计算平均IGD

    // 计算 QCCEANoQ IGD
    MILPIGD = 0.0;
    for (int j = 0; j < TotalParetoPopulation.size(); j++)
    {
        double minValue = std::numeric_limits<double>::max();
        for (int i = 0; i < MILPParetoPopulation.size(); i++)
        {
            double diffMS = MILPParetoPopulation[i].normalMS - TotalParetoPopulation[j].normalMS;
            double diffTEC = MILPParetoPopulation[i].normalTEC - TotalParetoPopulation[j].normalTEC;
            double value = sqrt(diffMS * diffMS + diffTEC * diffTEC);
            minValue = std::min(minValue, value);
        }
        MILPIGD += minValue;
    }
    MILPIGD /= TotalParetoPopulation.size();
}
