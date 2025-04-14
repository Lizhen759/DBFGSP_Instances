#include <iomanip>
#include <vector>
#include <algorithm>
#include <valarray>
#include <iostream>
#include "Individual.h"

Individual::Individual()
{

}

Individual::~Individual()
{

}

void Individual::Normalize(vector<Individual>& Pop, int& nadirpointMS, float& nadirpointTEC, int& idealpointMS, float& indealpointTEC)
{
    int minMS = INT_MAX;
    float minTEC = INT_MAX;

    int maxMS = INT_MIN;
    float maxTEC = INT_MIN;
    for (int i = 0; i < Pop.size(); i++)
    {
        if (Pop[i].MS < minMS)
            minMS = Pop[i].MS;

        if (Pop[i].MS > maxMS)
            maxMS = Pop[i].MS;

        if (Pop[i].TEC < minTEC)
            minTEC = Pop[i].TEC;

        if (Pop[i].TEC > maxTEC)
            maxTEC = Pop[i].TEC;

        //标志设为0
        Pop[i].flag = 0;

    }
    nadirpointMS = maxMS;
    nadirpointTEC = maxTEC;
    idealpointMS = minMS;
    indealpointTEC = minTEC;

    for (int i = 0; i < Pop.size(); i++)
    {
        Pop[i].normalMS = (static_cast<float>(Pop[i].MS - minMS) / static_cast<float>(maxMS - minMS));
        Pop[i].normalTEC = ((Pop[i].TEC - minTEC) / (maxTEC - minTEC));

    }

}

void Individual::calc_convergence_ind(vector<Individual>& Pop)
{
    float temp;
    for (int i = 0; i < Pop.size(); i++)
    {
        temp = 0;
        temp += (Pop[i].normalMS - 1.0) * (Pop[i].normalMS - 1.0);
        temp += (Pop[i].normalTEC - 1.0) * (Pop[i].normalTEC - 1.0);
        temp = sqrt(temp);
        temp = 1 / (temp + 1);
        Pop[i].convergence_ind = temp;
        //cout << "第" << i << "个个体的收敛指标：" << Pop[i].convergence_ind << endl;
    }

    sort(Pop.begin(), Pop.end(), [](const Individual& a, const Individual& b)
    {
        return a.convergence_ind < b.convergence_ind;
    });

}
//
//void Individual::calc_distribution_ind(vector<Individual>& Pop)
//{
//
//    float v1[maxpop];
//    float temp;
//    for (int i = 0; i < Pop.size(); i++)
//    {
//        v1[i] = 0;
//        v1[i] += (Pop[i].normalMS - 0.0) * (Pop[i].normalMS - 0.0);
//        v1[i] += (Pop[i].normalTEC - 0.0)* (Pop[i].normalTEC - 0.0);
//        v1[i] = sqrt(v1[i]);
//        //cout << "v[" << i << "]：" << v1[i] << endl;
//    }
//    for (int i = 0; i < Pop.size(); i++)
//    {
//        for (int j = i; j < Pop.size(); j++)
//        {
//            temp = 0;
//            temp += (Pop[i].normalMS - 0.0) * (Pop[j].normalMS - 0.0);
//            temp += (Pop[i].normalTEC - 0.0) * (Pop[j].normalTEC - 0.0);
//
//            if (v1[i] == 0 || v1[j] == 0)
//                temp = 1;
//            else
//                temp = temp / (v1[i] * v1[j]);
//
//            temp = 1 - temp;
//
//            Pop[i].distribution_ind[j] = temp;
//            Pop[j].distribution_ind[i] = temp;
//        }
//    }
//}

void Individual::calc_distribution_ind(vector<Individual>& Pop) {
    // 1. 检查 Pop 是否为空
    if (Pop.empty()) {
       cout<<"Pop is empty"<<endl;
    };

    // 手动初始化数组
    for (int i = 0; i < Pop.size(); i++) {
        for (int j = 0; j < Pop.size(); j++) {
            Pop[i].distribution_ind[j] = 0.0f; // 显式初始化
        }
    }
    // 3. 计算 v1[i] = sqrt(normalMS^2 + normalTEC^2)
    vector<float> v1(Pop.size(), 0.0f);
    for (int i = 0; i < Pop.size(); i++) {
        v1[i] = sqrt(Pop[i].normalMS * Pop[i].normalMS + Pop[i].normalTEC * Pop[i].normalTEC);
    }

    // 4. 计算余弦相似度并转换为距离
    for (int i = 0; i < Pop.size(); i++) {
        if (v1[i] == 0.0f) {
            // 该个体是理想点，与其他个体的距离设为 1.0
            for (int j = 0; j < Pop.size(); j++) {
                Pop[i].distribution_ind[j] = 1.0f;
                Pop[j].distribution_ind[i] = 1.0f;
            }
            continue;
        }

        for (int j = i; j < Pop.size(); j++) {
            if (v1[j] == 0.0f) {
                Pop[i].distribution_ind[j] = 1.0f;
                Pop[j].distribution_ind[i] = 1.0f;
                continue;
            }

            float cosine_sim = (Pop[i].normalMS * Pop[j].normalMS + Pop[i].normalTEC * Pop[j].normalTEC) / (v1[i] * v1[j]);
            cosine_sim = max(-1.0f, min(1.0f, cosine_sim)); // 确保在 [-1, 1] 范围内
            float dist = 1.0f - cosine_sim; // 距离 ∈ [0, 2]

            Pop[i].distribution_ind[j] = dist;
            Pop[j].distribution_ind[i] = dist;
        }
    }
}

//float Individual::calc_distribution_ind_QCCEA(vector<Individual>& CCEAPopulation, float normalMS, float normalTEC)
//{
//    float distribution_ind = 0.0;
//    int populationSize = CCEAPopulation.size();
//
//    if (populationSize > 1) {
//        // 重新计算归一化参数
//        float new_nadirMS = normalMS, new_idealMS = normalMS;
//        float new_nadirTEC = normalTEC, new_idealTEC = normalTEC;
//
//        for (const auto& ind : CCEAPopulation) {
//            if (ind.normalMS > new_nadirMS) new_nadirMS = ind.normalMS;
//            if (ind.normalMS < new_idealMS) new_idealMS = ind.normalMS;
//            if (ind.normalTEC > new_nadirTEC) new_nadirTEC = ind.normalTEC;
//            if (ind.normalTEC < new_idealTEC) new_idealTEC = ind.normalTEC;
//        }
//
//        const float epsilon = 1e-6;
//        if (fabs(new_nadirMS - new_idealMS) < epsilon || fabs(new_nadirTEC - new_idealTEC) < epsilon) {
//            return 0.0;
//        }
//
//        // 归一化所有个体
//        vector<pair<float, float>> normalized_population;
//        for (const auto& ind : CCEAPopulation) {
//            float normMS = (fabs(new_nadirMS - new_idealMS) < epsilon) ? 0.0 : (ind.normalMS - new_idealMS) / (new_nadirMS - new_idealMS);
//            float normTEC = (fabs(new_nadirTEC - new_idealTEC) < epsilon) ? 0.0 : (ind.normalTEC - new_idealTEC) / (new_nadirTEC - new_idealTEC);
//            normalized_population.emplace_back(normMS, normTEC);
//        }
//
//        // 归一化新解
//        float normMS = (fabs(new_nadirMS - new_idealMS) < epsilon) ? 0.0 : (normalMS - new_idealMS) / (new_nadirMS - new_idealMS);
//        float normTEC = (fabs(new_nadirTEC - new_idealTEC) < epsilon) ? 0.0 : (normalTEC - new_idealTEC) / (new_nadirTEC - new_idealTEC);
//
//        // 计算余弦相似度
//        float v1 = calculate_norm(normMS, normTEC);
//        float totalDistribution = 0.0;
//        int count = 0;
//
//        for (int i = 0; i < populationSize; i++) {
//            float popMS = normalized_population[i].first;
//            float popTEC = normalized_population[i].second;
//
//            if (fabs(popMS - normMS) < epsilon && fabs(popTEC - normTEC) < epsilon) {
//                continue;
//            }
//
//            float temp = (normMS * popMS) + (normTEC * popTEC);
//            float v2 = calculate_norm(popMS, popTEC);
//
//            temp = (v1 == 0 || v2 == 0) ? 0.0 : (1 - (temp / (v1 * v2)));
//
//            totalDistribution += temp;
//            count++;
//        }
//
//        distribution_ind = (count > 0) ? (totalDistribution / count) : 0.0;
//    }
//
//    return distribution_ind;
//}
//
//
//// 计算个体目标值的模长（L2范数）
//float Individual::calculate_norm(float normalMS, float normalTEC) {
//    return sqrt(normalMS * normalMS + normalTEC * normalTEC);
//}

