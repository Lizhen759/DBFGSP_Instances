#ifndef DBFGSP_PROBLEM_H
#define DBFGSP_PROBLEM_H

#include <string>
#include <vector>
using namespace std;
class Problem
{
public:
    void ReadInstanceFileNameList(string Dir);   //读取实例文件名
    void ReadInstance(int InsNo);   //读取实例
    Problem();

    virtual ~Problem();

protected:
    int m_Factories;
    int m_Machines;
    int m_Families;
    int m_Jobs;
    int m_SetupType;

    vector<string> m_InstanceFileNameList;
    vector<string> m_InstanceNameList;

    vector<vector<int>> m_JobsInEachFamily; //jobs in each family
    vector<vector<int>> m_JobSeqInFam; //Job sequence in each Family

    vector<vector<int>> m_JobOperPTime; //operation time on each machine for each job
    vector<int> m_JobTotalPTime; //all the operation time on all machines for a job;
    int m_AllJobTotalPTime;

    vector<vector<vector<int>>> m_SetupTime; //setup time between families
    vector<vector<int>> m_FamSumSetupTime; //准备时间的和 setup time between Families

    vector<vector<int>> m_TureJobOpertime;  //变速之后的处理时间
    vector<float> m_Speed;  //速度
    vector<vector<int>> m_SpeedMatrix;  //速度矩阵

    float UnitIdleEC;
    float UnitSetupEC;
    vector<vector<float>> UnitPEC;
    vector<float> TEC;

    vector<int> m_JobWeightTotalPTime;
    vector<double> m_FamTotalPTime; //all the operation time on all machines for the jobs in a family
    vector<double> m_FamWeightTotalPTime;
    vector<double> m_FamTotalPTimeOnFirstMachine; //the fam total processing time on the first machine
    vector<double> m_FamTotalPTimeOnLastMachine; //the fam total processing time on the last machine
    vector<vector<int>> m_FamMaxSetupTime; //
    vector<double> m_FamAvgMaxSetupTime; //average all the max setup times for Families
    vector<double> m_FamAvgSetupTime; //average all the setup times for families
    vector<double> m_FamTotalSkewness;

    vector<vector<int>> m_BestJobSeqInFam; //best solution
    vector<vector<int>> m_BestFacFamSeq;

    vector<vector<int>> m_JFDTime, m_TempJFDTime, m_BestJFDTime, m_tempJFDTime1;
    vector<vector<int>> m_JBDTime, m_TempJBDTime, m_BestJBDTime, m_tempJBDTime1;
    vector<vector<int>> m_Hierarchy, m_TempHierarchy, m_BestHierarchy,m_tempHierarchy1;


    vector<vector<vector<int>>> m_JFDTimePop, m_TempJFDTimePop, m_BestJFDTimePop, m_tempJFDTime1Pop;
    vector<vector<vector<int>>> m_JBDTimePop, m_TempJBDTimePop, m_BestJBDTimePop, m_tempJBDTime1Pop;
    vector<vector<vector<int>>> m_HierarchyPop, m_TempHierarchyPop, m_BestHierarchyPop,m_tempHierarchy1Pop;

    void GetJobTotalPTime();

    void GetJobWeightTotalPTime();

    void GetFamTotalPTime();

    void GetFamWeightTotalPTime();

    void GetFamTotalPTimeOnFirstMachine();

    void GetFamTotalPTimeOnLastMachine();

    void GetFamSumSetupTime();

    void GetFamMaxSetupTime();

    void GetFamAvgMaxSetupTime();

    void GetFamAvgSetupTime();

    void GetFamTotalSkewness();

    double CalculateSkewness(const vector<double> &data);

    void GetFamAvgSetupTime_QCCEA();

    void GetFamTotalPTimeOnLastMachine_QCCEA();

    void GetFamTotalPTimeOnFirstMachine_QCCEA();

    void GetFamTotalPTime_QCCEA();

    void  GetFamWeightTotalPTime_QCCEA();

};

#endif //DBFGSP_NEW_PROBLEM_H
