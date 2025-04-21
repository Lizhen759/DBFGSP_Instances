#ifndef PARETONORMALIZER_H
#define PARETONORMALIZER_H
#include <vector>
#include "Individual.h"
#include <unordered_set>
#include <tuple>
using namespace std;

// Pareto 归一化工具类
class ParetoOptimizer {
public:
    // 比较两个解的函数，用于排序
    static bool compareSolutions(const Individual& sol1, const Individual& sol2);


    void mergeParetoSets(int ins, vector<vector<Individual>> &TempTotalParetoSet,
                         const vector<vector<Individual>> &MOSAParetoSet,
                         const vector<vector<Individual>> &CMOEAParetoSet,
                         const vector<vector<Individual>> &CCEAParetoSet,
                         const vector<vector<Individual>> &TSCCEAParetoSet,
                         const vector<vector<Individual>> &QCCEAParetoSet);
    // 计算Pareto关系
    void mainPareto_relation(std::vector<Individual>& Population);

    //弱化支配
    void deemphasizeDominatedSolutions(int ins, vector<vector<Individual>> &TempTotalParetoSet,vector<vector<Individual>> &TotalParetoSet);

    //去掉重复
    void removeDuplicateSolutions(int ins, vector<vector<Individual>> &TotalParetoSet,vector<vector<Individual>> &FinalTotalParetoSet);

    // 归一化函数
    void mainNormalize(std::vector<Individual>& totalPopulation,
                   std::vector<Individual>& MOSAPopulation,
                   std::vector<Individual>& CMOEAPopulation,
                   std::vector<Individual>& GCCEAPopulation,
                   std::vector<Individual>& TSCCPopulation,
                   std::vector<Individual>& QCCEAPopulation);

    // 归一化函数（更复杂的版本）
    void AdPrummainNormalize(std::vector<Individual>& TotalParetoPopulation,
                             std::vector<Individual>& ParetoPopulation1, std::vector<Individual>& ParetoPopulation2,
                             std::vector<Individual>& ParetoPopulation3, std::vector<Individual>& ParetoPopulation4,
                             std::vector<Individual>& ParetoPopulation5, std::vector<Individual>& ParetoPopulation6,
                             std::vector<Individual>& ParetoPopulation7, std::vector<Individual>& ParetoPopulation8,
                             std::vector<Individual>& ParetoPopulation9);

    // 计算C值
    void ComputeCoverage(const std::vector<Individual>& TotalParetoPopulation,
                  const std::vector<Individual>& MOSAParetoPopulation, const std::vector<Individual>& CMOEAParetoPopulation,
                  const std::vector<Individual>& GCCEAParetoPopulation, const std::vector<Individual>& TSCCEAParetoPopulation, const std::vector<Individual>& QCCEAParetoPopulation,
                  double& MOSAC, double& CMOEAC, double& GCCEAC, double& TSCCEAC,double& QCCEAC);

    // 计算超体积（HV）
    void ComputeHV( const std::vector<Individual>& TotalParetoPopulation,
                   std::vector<Individual>& MOSAParetoPopulation, std::vector<Individual>& CMOEAParetoPopulation,
                   std::vector<Individual>& GCCEAParetoPopulation, std::vector<Individual>& TSCCEAParetoPopulation, std::vector<Individual>& QCCEAParetoPopulation,
                   double& MOSAhv, double& CMOEAhv, double& GCCEAhv, double& TSCCEAhv, double& QCCEAhv);

    // 计算IGD（指标值）
    void ComputeIGD( const std::vector<Individual>& TotalParetoPopulation, const std::vector<Individual>& MOSAParetoPopulation,
                    const std::vector<Individual>& CMOEAParetoPopulation, const std::vector<Individual>& GCCEAParetoPopulation,
                    const std::vector<Individual>& TSCCEAParetoPopulation,const std::vector<Individual>& QCCEAParetoPopulation,
                    double& MOSAIGD, double& CMOEAIGD, double& GCCEAIGD, double& TSCCEAIGD,double& QCCEAIGD);

    void
    AdPrumComputeIGD( const vector<Individual> &TotalParetoPopulation,
                     const vector<Individual> &ParetoPopulation1,
                     const vector<Individual> &ParetoPopulation2, const vector<Individual> &ParetoPopulation3,
                     const vector<Individual> &ParetoPopulation4, const vector<Individual> &ParetoPopulation5,
                     const vector<Individual> &ParetoPopulation6, const vector<Individual> &ParetoPopulation7,
                     const vector<Individual> &ParetoPopulation8, const vector<Individual> &ParetoPopulation9,
                     double &IGD1,
                     double &IGD2, double &IGD3, double &IGD4, double &IGD5, double &IGD6, double &IGD7, double &IGD8,
                     double &IGD9);



    void AdPrumComputeIGD_Internal(const vector<Individual> &TotalParetoPopulation,
                                   const vector<vector<Individual>> &AllParetoPopulations, vector<double> &IGDs);

    void mergeParetoSets(int ins, vector<vector<Individual>> &TempTotalParetoSet,
                         const vector<vector<Individual>> &QCCEAParetoSet,
                         const vector<vector<Individual>> &QCCEANoQParetoSet,
                         const vector<vector<Individual>> &QCCEANoRapidParetoSet);

    void mainNormalize(vector<Individual> &totalPopulation, vector<Individual> &QCCEAPopulation,
                       vector<Individual> &QCCEANoQParetoSet, vector<Individual> &QCCEANoRapidParetoSet);

    void ComputeCoverage_Str( const vector<Individual> &TotalParetoPopulation,
                         const vector<Individual> &QCCEAParetoPopulation,
                         const vector<Individual> &QCCEANoQParetoPopulation,
                         const vector<Individual> &QCCEANoRapidParetoPopulation, double &QCCEAC, double &QCCEANoQC,
                         double &QCCEANoRapidC);

    void
    ComputeIGD_Str( const vector<Individual> &TotalParetoPopulation,
               const vector<Individual> &QCCEAParetoPopulation,
               const vector<Individual> &QCCEANoQParetoPopulation,
               const vector<Individual> &QCCEANoRapidParetoPopulation,
               double &QCCEAIGD, double &QCCEANoQIGD, double &QCCEANoRapidIGD);

    void ComputeHV_Str( const vector<Individual> &TotalParetoPopulation, vector<Individual> &QCCEAParetoPopulation,
                   vector<Individual> &QCCEANoQParetoPopulation, vector<Individual> &QCCEANoRapidParetoPopulation,
                   double &QCCEAIGD, double &QCCEANoQIGD, double &QCCEANoRapidIGD);

    void mergeParetoSets(int ins, vector<vector<Individual>> &TempTotalParetoSet,
                         const vector<vector<Individual>> &QCCEAParetoSet,
                         const vector<vector<Individual>> &MILPParetoSet);

    void mainNormalize(vector<Individual> &totalPopulation, vector<Individual> &QCCEAPopulation,
                       vector<Individual> &MILPParetoSet);

    void ComputeCoverage_MILP( const vector<Individual> &TotalParetoPopulation,
                         const vector<Individual> &QCCEAParetoPopulation,
                         const vector<Individual> &MILPParetoPopulation,
                         double &QCCEAC, double &MILPC);

    void
    ComputeIGD_MILP( const vector<Individual> &TotalParetoPopulation,
               const vector<Individual> &QCCEAParetoPopulation,
               const vector<Individual> &MILPParetoPopulation, double &QCCEAIGD, double &MILPIGD);

    void ComputeHV_MILP(const vector<Individual> &TotalParetoPopulation, vector<Individual> &QCCEAParetoPopulation,
                   vector<Individual> &MILPParetoPopulation, double &QCCEAHV, double &MILPHV);

};


#endif // PARETONORMALIZER_H
