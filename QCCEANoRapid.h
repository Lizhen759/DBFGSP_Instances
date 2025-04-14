#ifndef DBFGSP_NEW_QCCEANoRapid_H
#define DBFGSP_NEW_QCCEANoRapid_H
#pragma once

#include <unordered_map>
#include "Individual.h"
#include "NOperator.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include "Base.h"
#include <valarray>
#include <sstream>
#include <numeric>
#include "rand.h"

class QCCEANoRapid :
        public NOperator, public Individual {
public:
    QCCEANoRapid();

    virtual ~ QCCEANoRapid();

    int RunEvolution(int CPUTime, vector<vector<Individual>> &QCCEANoRapidFinalAfterRepParetoSet, int AN, int mu);

private:

    long m_TimeLimit;

    double e = 0.85; //贪婪策略
    double discountFactor = 0.2; //折扣因子
    double rewardFactor = 0.75;//学习率

    //asdf
    int m_RefSize;
    int m_PopSize;

    double m_T;

    float thr_zeta;

    int m_NadirPointMS;    //最低点
    float m_NadirPointTEC;

    int m_IdealPointMS;      //理想点
    float m_IdealPointTEC;

    //档案集
    vector<vector<vector<int>>> m_ArchiveSpeedVector;
    vector<vector<vector<int>>> m_ArchiveFacFamSeqArray;
    vector<vector<vector<int>>> m_ArchiveJobSeqinFamArray;
    vector<vector<int>> m_ArchiveFacSpanArray;
    vector<vector<float>> m_ArchiveFacECArray;
    vector<int> m_ArchiveSpanArray;
    vector<float> m_ArchiveTECArray;

    vector<Individual> m_QCCEANoRapidPopulation;
    vector<Individual> m_QCCEANoRapidParetoSet;  //ParetoSet解集

    vector<vector<vector<int>>> m_FacFamSeqArray;
    vector<vector<vector<int>>> m_JobSeqinFamArray;
    vector<int> m_SpanArray1;
    vector<int> m_SpanArray2;
    vector<float> m_TECArray1;
    vector<float> m_TECArray2;

    int m_PS1;
    int m_PS2;
    vector<int> m_Map1;
    vector<int> m_Map2;

    void SetParameters(int AN,  long TimeLimit, double T);

    void InitialPop();

    void EvolutionProcess(int mu);

    /*********************组进化************************/
    void BasedInd_RandFamInFacTobestPos (vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                            int NadirPointMS, float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                                            vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS, float &ObjectTEC);

    void BasedInd_SwapFam(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                              int NadirPointMS, float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                              vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS, float &ObjectTEC);

    void BasedInd_ShiftFam(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,  int NadirPointMS, float NadirPointTEC, int IdealPointMS,
                               float IdealPointTEC,
                               vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC, int &ObjectMS,
                               float &ObjectTEC );

    void BasedIndLS_SetupSwap(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,int NadirPointMS,
                                  float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                                  vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS, float &ObjectTEC);

    void BasedInd_DeAndConFams(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,int NadirPointMS,
                                   float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                                   vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS, float &ObjectTEC);


    /**********************工件进化********************************/
    void BasedInd_CirJobInFamTobestPos_New(const vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                           vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                           int NadirPointMS, float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                                           vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS, float &ObjectTEC);
    void
    JobLocalSearch(const vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam, int NadirPointMS, float NadirPointTEC, int IdealPointMS,
                   float IdealPointTEC, vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC,
                   int &ObjectMS, float &ObjectTEC);

    void BasedInd_SwapJob(const vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                              int NadirPointMS, float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                              vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS, float &ObjectTEC);

    void BasedInd_JobsInFam(int Fac, const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,
                                vector<vector<int>> &JobSeqInFam,int Fam, int NadirPointMS, float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                                vector<Individual> &CCEAPopulation, int OrgMS, float OrgTEC, int &ObjectMS,
                                float &ObjectTEC);

    /*******************************参考集进化************************************/
    void BasedInd_RefDeconAndCon_New(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                     int NadirPointMS, float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                                     vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS,
                                     float &ObjectTEC, vector<double> &scores_combine);

    void  ApplyDeconAndCon_Ref (int DeAndConIndex ,vector<vector<int>> &Sol, vector<vector<int>> &JobSeqInFam,
                                vector<int> &FamsExtracted,unordered_map<int, vector<int>> &JobsExtracted, int NadirPointMS,float NadirPointTEC, int IdealPointMS,
                                float IdealPointTEC, vector<int>& FacSpan, vector<float>& FacEC, vector<Individual> &CCEAPopulation);


    void BasedInd_Destruction_FamsAndJobs_New_Opeator0(vector<vector<int>>& FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                       vector<int> &FamsExtracted,
                                                       unordered_map<int, vector<int>> &JobsExtracted);

    void BasedInd_Destruction_FamsAndJobs_New_Opeator1(vector<vector<int>>& FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                       vector<int> &FamsExtracted,
                                                       unordered_map<int, vector<int>> &JobsExtracted, vector<int>& FacSpan, vector<float>& FacEC);

    void BasedInd_Destruction_FamsAndJobs_New_Opeator2(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                       vector<int> &FamsExtracted,
                                                       unordered_map<int, vector<int>> &JobsExtracted,
                                                       vector<int> &FacSpan, vector<float> &FacEC);

    void BasedInd_Construction_FamsAndJobs_New_Opeator0(vector<vector<int>>& FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                        vector<int> &FamsExtracted,
                                                        unordered_map<int, vector<int>> &JobsExtracted, int NadirPointMS,
                                                        float NadirPointTEC, int IdealPointMS, float IdealPointTEC ,vector<Individual> &CCEAPopulation);

    void BasedInd_Construction_FamsAndJobs_New_Opeator1(vector<vector<int>>& FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                        vector<int> &FamsExtracted,
                                                        unordered_map<int, vector<int>> &JobsExtracted, int NadirPointMS,
                                                        float NadirPointTEC, int IdealPointMS, float IdealPointTEC ,vector<Individual> &CCEAPopulation);

    void BasedInd_Construction_FamsAndJobs_New_Opeator2(vector<vector<int>>& FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                        vector<int> &FamsExtracted,
                                                        unordered_map<int, vector<int>> &JobsExtracted, int NadirPointMS,
                                                        float NadirPointTEC, int IdealPointMS, float IdealPointTEC ,vector<Individual> &CCEAPopulation);

    static bool SA( int OrgMS, float OrgTEC, int AfterMS, float AfterTEC, double temperature);

    /*******************************档案集进化************************************/
    void UpdateArchiveGroupJobSet(int mu);

    void SaveTECandDeMS(vector<Individual> &CCEAPopulation);

    int DetermineAction(int state, const vector<vector<double>> &Q);

    vector<int> findMaxIndices(const vector<double> &vec);

    static int DetermineState( int ObjectMS, float ObjectTEC, const vector<Individual> &QCCEANoRapidPopulation);

    void UpdateQ(int oldstate, int newstate, int act, double reward, vector<vector<double>> &Q);

    static double CaculateReward(int oldMS, int newMS, float oldTEC, float newTEC, const vector<Individual> &population);

    void PerformAction_Fams(int act, vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqinFam,
                            vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC,
                            int &ObjectMS, float &ObjectTEC);

    static bool CompareByMakespan(const Individual &a, const Individual &b);

    static bool CompareByTEC(const Individual &a, const Individual &b);

};

#endif //DBFGSP_NEW_QCCEANoRapid_H
