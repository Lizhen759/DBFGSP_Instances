#include <iostream>
#include <ostream>
#include <cassert>
#include "Problem.h"
#include "Individual.h"
#include "TSCCEA.h"
#include "GCCEA.h"
#include "ParetoReader.h"
#include "ParetoOptimizer.h"
#include "MOSA.h"
#include "CMOEA.h"
#include "QCCEA.h"
#include "QCCEANoQ.h"
#include "QCCEANoRapid.h"

void AdjustPerum();  // 先声明函数
void Compare();
void CompareStr();
void CompareMILP();

int main()
{
    // 初始化 TSCCEAFinalAfterRepParetoSet
    vector<vector<Individual>> TSCCEAFinalAfterRepParetoSet(405);
//    TSCCEA tsccea;
//    tsccea.RunEvolution(20, 10, TSCCEAFinalAfterRepParetoSet);

    // 初始化 CMOEAFinalAfterRepParetoSet
    vector<vector<Individual>> CMOEAFinalAfterRepParetoSet(405);
//    CMOEA cmoea;
//    cmoea.RunEvolution(20 , CMOEAFinalAfterRepParetoSet, 10, 4);

    // 初始化 MOSAFinalAfterRepParetoSet
    vector<vector<Individual>> MOSAFinalAfterRepParetoSet(405);
//    MOSA mosa;
//    mosa.RunEvolution(20, 10, MOSAFinalAfterRepParetoSet);

    // 初始化GCCEAFinalAfterRepParetoSet
    vector<vector<Individual>> GCCEAFinalAfterRepParetoSet(405);
//    GCCEA gccea;
//    gccea.RunEvolution(20, 10, GCCEAFinalAfterRepParetoSet);


    vector<vector<Individual>> QCCEAFinalAfterRepParetoSet(405);
    QCCEA qccea;
//    qccea.RunEvolution(20, QCCEAFinalAfterRepParetoSet, 10, 4);

    vector<vector<Individual>> QCCEANoQFinalAfterRepParetoSet(405);
//    QCCEANoQ qcceaNoQ;
//    qcceaNoQ.RunEvolution(20, QCCEANoQFinalAfterRepParetoSet, 10, 4);

    vector<vector<Individual>> QCCEANoRapidFinalAfterRepParetoSet(405);
//    QCCEANoRapid qcceaNoRapid;
//    qcceaNoRapid.RunEvolution(20, QCCEAANoRapidFinalAfterRepParetoSet, 10, 4);

//      Compare();
//    AdjustPerum();
//    CompareStr();
    CompareMILP();
    return 0;
}
//对比算法
void Compare(){
    //暂时
    vector<vector<Individual>> TempTotalParetoSet(405);

    //全部
    vector<vector<Individual>> TotalParetoSet(405);

    //最后
    vector<vector<Individual>> FinalTotalParetoSet(405);

    // 读取文件中的 TSCCEAPareto 解集
    std::vector<std::vector<Individual>> TSCCEAparetoSets(405);
    ParetoReader TSCCEAreader;
    TSCCEAreader.readParetoSetFromFile("../Result/TSCCEA_20_experiment.txt", TSCCEAparetoSets);

    // 读取文件中的 CMOEAPareto 解集
    std::vector<std::vector<Individual>> MOSAparetoSets(405);
    ParetoReader MOSAreader;
    MOSAreader.readParetoSetFromFile("../Result/MOSA_20_experiment.txt", MOSAparetoSets);

    // 读取文件中的 CMOEAPareto 解集
    std::vector<std::vector<Individual>> GCCEAparetoSets(405);
    ParetoReader GCCEAreader;
    GCCEAreader.readParetoSetFromFile("../Result/GCCEA_20_experiment.txt",  GCCEAparetoSets);

    // 读取文件中的 CMOEAPareto 解集
    std::vector<std::vector<Individual>> CMOEAparetoSets(405);
    ParetoReader CMOEAreader;
    CMOEAreader.readParetoSetFromFile("../Result/CMOEA_20_experiment.txt", CMOEAparetoSets);

    // 读取文件中的 CMOEAPareto 解集
    std::vector<std::vector<Individual>> QCCEAparetoSets(405);
    ParetoReader QCCEAreader;
    QCCEAreader.readParetoSetFromFile("../Result/QCCEA_20_experiment.txt", QCCEAparetoSets);

    ofstream ofileC;
    ofileC.open("../Result/Ind_Coverage.txt");  // 打开Ind_Coverage文件进行写入
    if (!ofileC.is_open())
    {
        cout << "../Result/Ind_Coverage.txt" << "\t is not open" << endl;
        exit(0);
    }

    ofstream ofileIGD;
    ofileIGD.open("../Result/Ind_IGD.txt");  // 打开Ind_IGD.txt文件进行写入
    if (!ofileIGD.is_open())
    {
        cout << "../Result/Ind_IGD.txt.txt" << "\t is not open" << endl;
        exit(0);
    }

    ofstream ofileHV;
    ofileHV.open("../Result/Ind_HV.txt");  // 打开Ind_HV.txt文件进行写入
    if (!ofileHV.is_open())
    {
        cout << "../Result/Ind_HV.txt" << "\t is not open" << endl;
        exit(0);
    }

    for (int ins = 0; ins < 405; ins++)
    {

        ParetoOptimizer paretoOptimizer;  // 创建 ParetoOptimizer 实例

        //合并
        TempTotalParetoSet[ins].clear(); // 清空当前 TempTotalParetoSet
        paretoOptimizer.mergeParetoSets(ins,TempTotalParetoSet,MOSAparetoSets,CMOEAparetoSets,
                                        GCCEAparetoSets,TSCCEAparetoSets,QCCEAparetoSets);

        // 寻找非支配解
        paretoOptimizer.mainPareto_relation(TempTotalParetoSet[ins]);
        // 弱化非支配解
        TotalParetoSet[ins].clear();
        paretoOptimizer.deemphasizeDominatedSolutions(ins,TempTotalParetoSet,TotalParetoSet);
        //去除重复
        FinalTotalParetoSet[ins].clear();
        paretoOptimizer.removeDuplicateSolutions(ins,TotalParetoSet,FinalTotalParetoSet);


        // 归一化
        paretoOptimizer.mainNormalize(FinalTotalParetoSet[ins], MOSAparetoSets[ins], CMOEAparetoSets[ins],GCCEAparetoSets[ins],
                                      TSCCEAparetoSets[ins], QCCEAparetoSets[ins]);
        // 计算 Coverage
        double MOSACoverage= 0.0;
        double CMOEACoverage= 0.0;
        double GCCEACoverage = 0.0;
        double TSCCEACoverage = 0.0;
        double QCCEACoverage = 0.0;
        paretoOptimizer.ComputeCoverage(FinalTotalParetoSet[ins], MOSAparetoSets[ins], CMOEAparetoSets[ins],
                                        GCCEAparetoSets[ins], TSCCEAparetoSets[ins],QCCEAparetoSets[ins], MOSACoverage, CMOEACoverage,
                                        GCCEACoverage, TSCCEACoverage,QCCEACoverage);
        ofileC << MOSACoverage << "," << CMOEACoverage << "," << GCCEACoverage << "," << TSCCEACoverage <<","<< QCCEACoverage <<endl;

        // 计算 IGD
        double MOSAIGD = 0.0;
        double CMOEAIGD = 0.0;
        double GCCEAIGD = 0.0;
        double TSCCEAIGD = 0.0;
        double QCCEAIGD = 0.0;
        paretoOptimizer.ComputeIGD(FinalTotalParetoSet[ins], MOSAparetoSets[ins], CMOEAparetoSets[ins],
                                   GCCEAparetoSets[ins], TSCCEAparetoSets[ins],QCCEAparetoSets[ins], MOSAIGD, CMOEAIGD, GCCEAIGD, TSCCEAIGD,QCCEAIGD);
        ofileIGD << MOSAIGD << "," << CMOEAIGD << "," << GCCEAIGD << "," << TSCCEAIGD << ","<< QCCEAIGD <<endl;

        // 计算 HV
        double MOSAHV = 0.0;
        double CMOEAHV = 0.0;
        double GCCEAHV = 0.0;
        double TSCCEAHV = 0.0;
        double QCCEAHV = 0.0;
        paretoOptimizer.ComputeHV(FinalTotalParetoSet[ins], MOSAparetoSets[ins], CMOEAparetoSets[ins],
                                  GCCEAparetoSets[ins], TSCCEAparetoSets[ins], QCCEAparetoSets[ins],MOSAHV, CMOEAHV, GCCEAHV, TSCCEAHV, QCCEAHV);
        ofileHV << MOSAHV << "," << CMOEAHV << "," << GCCEAHV << "," << TSCCEAHV << ","<< QCCEAHV <<endl;
    }
}

//策略验证
void CompareStr(){
    //暂时
    vector<vector<Individual>> TempTotalParetoSet(405);

    //全部
    vector<vector<Individual>> TotalParetoSet(405);

    //最后
    vector<vector<Individual>> FinalTotalParetoSet(405);

    // 读取文件中的QCCEAPareto 解集
    std::vector<std::vector<Individual>> QCCEAparetoSets(405);

    ParetoReader QCCEAreader;
    QCCEAreader.readParetoSetFromFile("../StrResult/QCCEA_20_experiment.txt", QCCEAparetoSets);

    // 读取文件中的 CMOEAPareto 解集
    std::vector<std::vector<Individual>> QCCEANoQparetoSets(405);
    ParetoReader QCCEANoQreader;
    QCCEANoQreader.readParetoSetFromFile("../StrResult/QCCEANoQ_20_experiment.txt", QCCEANoQparetoSets);

    // 读取文件中的 CMOEAPareto 解集
    std::vector<std::vector<Individual>> QCCEANoRapidparetoSets(405);
    ParetoReader QCCEANoRapidreader;
    QCCEANoRapidreader.readParetoSetFromFile("../StrResult/QCCEANoRapid_20_experiment.txt", QCCEANoRapidparetoSets);


    ofstream ofileC;
    ofileC.open("../StrResult/CompareStr_Ind_Coverage.txt");
    if (!ofileC.is_open())
    {
        cout << "../StrResult/CompareStr_Ind_Coverage.txt" << "\t is not open" << endl;
        exit(0);
    }

    ofstream ofileIGD;
    ofileIGD.open("../StrResult/CompareStr_Ind_IGD.txt");  // 打开Ind_IGD.txt文件进行写入
    if (!ofileIGD.is_open())
    {
        cout << "../StrResult/CompareStr_Ind_IGD.txt" << "\t is not open" << endl;
        exit(0);
    }

    ofstream ofileHV;
    ofileHV.open("../StrResult/CompareStr_Ind_HV.txt");  // 打开Ind_HV.txt文件进行写入
    if (!ofileHV.is_open())
    {
        cout << "../StrResult/CompareStr_Ind_HV.txt" << "\t is not open" << endl;
        exit(0);
    }

    for (int ins = 0; ins < 405; ins++)
    {

        ParetoOptimizer paretoOptimizer;  // 创建 ParetoOptimizer 实例

        //合并
        TempTotalParetoSet[ins].clear(); // 清空当前 TempTotalParetoSet
        paretoOptimizer.mergeParetoSets(ins,TempTotalParetoSet,QCCEAparetoSets,QCCEANoQparetoSets,QCCEANoRapidparetoSets);

        // 寻找非支配解
        paretoOptimizer.mainPareto_relation(TempTotalParetoSet[ins]);
        // 弱化非支配解
        TotalParetoSet[ins].clear();
        paretoOptimizer.deemphasizeDominatedSolutions(ins,TempTotalParetoSet,TotalParetoSet);
        //去除重复
        FinalTotalParetoSet[ins].clear();
        paretoOptimizer.removeDuplicateSolutions(ins,TotalParetoSet,FinalTotalParetoSet);

        // 归一化
        paretoOptimizer.mainNormalize(FinalTotalParetoSet[ins], QCCEAparetoSets[ins], QCCEANoQparetoSets[ins],QCCEANoRapidparetoSets[ins] );
        // 计算 Coverage
        double QCCEACoverage = 0.0;
        double QCCEANoQCoverage = 0.0;
        double QCCEANoRapidCoverage = 0.0;
        paretoOptimizer.ComputeCoverage_Str(FinalTotalParetoSet[ins], QCCEAparetoSets[ins],
                                        QCCEANoQparetoSets[ins],QCCEANoRapidparetoSets[ins],
                                        QCCEACoverage,QCCEANoQCoverage,QCCEANoRapidCoverage );
        ofileC << QCCEANoQCoverage << "," << QCCEANoRapidCoverage <<","<< QCCEACoverage <<endl;

        // 计算 IGD
        double QCCEAIGD = 0.0;
        double QCCEANoQIGD = 0.0;
        double QCCEANoRapidIGD = 0.0;
        paretoOptimizer.ComputeIGD_Str(FinalTotalParetoSet[ins], QCCEAparetoSets[ins],
                                   QCCEANoQparetoSets[ins], QCCEANoRapidparetoSets[ins],
                                   QCCEAIGD,QCCEANoQIGD,QCCEANoRapidIGD);
        ofileIGD << QCCEANoQIGD << "," << QCCEANoRapidIGD<< "," << QCCEAIGD <<endl;

        // 计算 HV
        double QCCEAHV = 0.0;
        double QCCEANoQHV = 0.0;
        double QCCEANoRapidHV = 0.0;
        paretoOptimizer.ComputeHV_Str( FinalTotalParetoSet[ins],  QCCEAparetoSets[ins],QCCEANoQparetoSets[ins], QCCEANoRapidparetoSets[ins],
                                   QCCEAHV,QCCEANoQHV,QCCEANoRapidHV);
        ofileHV <<QCCEANoQHV << "," << QCCEANoRapidHV << ","<< QCCEAHV <<endl;
    }
}

//调参
void AdjustPerum()
{
    vector<vector<vector<Individual>>> QCCEAResults(113);
    ParetoReader QCCEAreader;

    // 读取多个文件的 Pareto 集
    for (int i = 0; i < 9; i++) {
        string filePath = "../AdPrum_Result/" + to_string(i + 1) + ".txt";

        QCCEAreader.readParetoSetFromFile(filePath, QCCEAResults[i]);
    }


    vector<vector<Individual>> TempTotalParetoSet(113);
    vector<vector<Individual>> TotalParetoSet(113);
    vector<vector<Individual>> FinalTotalParetoSet(113);

    ofstream ofileIGD("../AdPrum_Result/AdPrum_Ind_IGD.txt");
    if (!ofileIGD.is_open()) {
        cout << "../AdPrum_Result/AdPrum_Ind_IGD.txt" << "\t is not open" << endl;
        exit(0);
    }

    ParetoOptimizer paretoOptimizer;

    for (int ins = 0; ins < 113; ins++) {



        TempTotalParetoSet[ins].clear();

        // 读取所有解并输出
        for (int fileIdx = 0; fileIdx < 9; fileIdx++) {
            for (const auto& indiv : QCCEAResults[fileIdx][ins]) {
                TempTotalParetoSet[ins].push_back(indiv);
             //   cout << indiv.MS << " " << indiv.TEC << endl;
            }
          //  cout << endl;
        }

        // 处理非支配解
        TotalParetoSet[ins].clear();
        paretoOptimizer.mainPareto_relation(TempTotalParetoSet[ins]);

        // 标记被支配的解
        for (auto& indiv : TempTotalParetoSet[ins]) {
            if (indiv.flag == 0) {
                for (int j = 0; j < TempTotalParetoSet[ins].size(); j++) {
                    if (indiv.pareto_rel[j] == 1) {
                        indiv.flag = 999;
                        break;
                    }
                }
            }
        }

        // 仅保留非支配解
        for (const auto& indiv : TempTotalParetoSet[ins]) {
            if (indiv.flag == 0) {
                TotalParetoSet[ins].push_back(indiv);
            }
        }

        // 去重
        FinalTotalParetoSet[ins].clear();
        unordered_set<string> seenSolutions;
        for (const auto& indiv : TotalParetoSet[ins]) {
            string key = to_string(indiv.MS) + "," + to_string(indiv.TEC);
            if (seenSolutions.insert(key).second) {
                FinalTotalParetoSet[ins].push_back(indiv);
            }
        }

        // 归一化
        paretoOptimizer.AdPrummainNormalize(FinalTotalParetoSet[ins],
                                            QCCEAResults[0][ins], QCCEAResults[1][ins],
                                            QCCEAResults[2][ins], QCCEAResults[3][ins],
                                            QCCEAResults[4][ins], QCCEAResults[5][ins],
                                            QCCEAResults[6][ins], QCCEAResults[7][ins],
                                            QCCEAResults[8][ins]);

        // 计算 IGD
        vector<double> IGDs(9, 0.0);
        paretoOptimizer.AdPrumComputeIGD( FinalTotalParetoSet[ins],
                                         QCCEAResults[0][ins], QCCEAResults[1][ins],
                                         QCCEAResults[2][ins], QCCEAResults[3][ins],
                                         QCCEAResults[4][ins], QCCEAResults[5][ins],
                                         QCCEAResults[6][ins], QCCEAResults[7][ins],
                                         QCCEAResults[8][ins],
                                         IGDs[0], IGDs[1], IGDs[2], IGDs[3],
                                         IGDs[4], IGDs[5], IGDs[6], IGDs[7], IGDs[8]);

        // 输出 IGD
        //ofileIGD << ins + 1;
        for (double igd : IGDs) {
            ofileIGD << "," << igd;
        }
        ofileIGD << endl;
    }
}


//策略验证
void CompareMILP(){
    //暂时
    vector<vector<Individual>> TempTotalParetoSet(405);

    //全部
    vector<vector<Individual>> TotalParetoSet(405);

    //最后
    vector<vector<Individual>> FinalTotalParetoSet(405);

    // 读取文件中的QCCEAPareto 解集
    std::vector<std::vector<Individual>> QCCEAparetoSets(405);
    ParetoReader QCCEAreader;
    QCCEAreader.readParetoSetFromFile("../MILPResult/QCCEA_20_experiment.txt", QCCEAparetoSets);

    // 读取文件中的 MILPPareto 解集
    std::vector<std::vector<Individual>> MILPparetoSets(405);
    ParetoReader MILPreader;
    MILPreader.readParetoSetFromFile("../MILPResult/MILP_experiment.txt", MILPparetoSets);


    ofstream ofileC;
    ofileC.open("../MILPResult/CompareMILP_Ind_Coverage.txt");
    if (!ofileC.is_open())
    {
        cout << "../MILPResult/CompareMILP_Ind_Coverage.txt" << "\t is not open" << endl;
        exit(0);
    }

    ofstream ofileIGD;
    ofileIGD.open("../MILPResult/CompareMILP_Ind_IGD.txt");  // 打开Ind_IGD.txt文件进行写入
    if (!ofileIGD.is_open())
    {
        cout << "../MILPResult/CompareMILP_Ind_IGD.txt" << "\t is not open" << endl;
        exit(0);
    }

    ofstream ofileHV;
    ofileHV.open("../MILPResult/CompareMILP_Ind_HV.txt");  // 打开Ind_HV.txt文件进行写入
    if (!ofileHV.is_open())
    {
        cout << "../MILPResult/CompareMILP_Ind_HV.txt" << "\t is not open" << endl;
        exit(0);
    }

    for (int ins = 0; ins < 405; ins++)
    {

        ParetoOptimizer paretoOptimizer;  // 创建 ParetoOptimizer 实例

        //合并
        TempTotalParetoSet[ins].clear(); // 清空当前 TempTotalParetoSet
        paretoOptimizer.mergeParetoSets(ins, TempTotalParetoSet, QCCEAparetoSets, MILPparetoSets);

        // 寻找非支配解
        paretoOptimizer.mainPareto_relation(TempTotalParetoSet[ins]);
        // 弱化非支配解
        TotalParetoSet[ins].clear();
        paretoOptimizer.deemphasizeDominatedSolutions(ins,TempTotalParetoSet,TotalParetoSet);
        //去除重复
        FinalTotalParetoSet[ins].clear();
        paretoOptimizer.removeDuplicateSolutions(ins,TotalParetoSet,FinalTotalParetoSet);

        // 归一化
        paretoOptimizer.mainNormalize(FinalTotalParetoSet[ins], QCCEAparetoSets[ins], MILPparetoSets[ins] );
        // 计算 Coverage
        double QCCEACoverage = 0.0;
        double MILPCoverage = 0.0;
        paretoOptimizer.ComputeCoverage_MILP( FinalTotalParetoSet[ins], QCCEAparetoSets[ins], MILPparetoSets[ins],
                                        QCCEACoverage, MILPCoverage );
        ofileC << MILPCoverage << "," << QCCEACoverage << endl;

        // 计算 IGD
        double QCCEAIGD = 0.0;
        double MILPIGD = 0.0;
        paretoOptimizer.ComputeIGD_MILP( FinalTotalParetoSet[ins], QCCEAparetoSets[ins],MILPparetoSets[ins], QCCEAIGD, MILPIGD);
        ofileIGD << MILPIGD << "," << QCCEAIGD << endl;

        // 计算 HV
        double QCCEAHV = 0.0;
        double MILPHV = 0.0;
        paretoOptimizer.ComputeHV_MILP(FinalTotalParetoSet[ins], QCCEAparetoSets[ins],  MILPparetoSets[ins],QCCEAHV, MILPHV);
        ofileHV  << MILPHV << "," << QCCEAHV << endl;
    }
}
