//
// Created by 13954 on 2025/3/31.
//

#ifndef DBFGSP_NEW_QCCEANOQ_H
#define DBFGSP_NEW_QCCEANOQ_H

#pragma once

#include <unordered_map>
#include "Individual.h"
#include "NOperator_New.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include "Base.h"
#include <valarray>
#include <sstream>
#include <numeric>
#include "rand.h"

class QCCEANoQ :
        public NOperator_New, public Individual {
public:
    QCCEANoQ();

    virtual ~ QCCEANoQ();

    int RunEvolution(int CPUTime, vector<vector<Individual>> &QCCEANoQFinalAfterRepParetoSet, int AN, int mu);

private:

    long m_TimeLimit;
    

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

    vector<Individual> m_QCCEAPopulation;
    vector<Individual> m_QCCEAParetoSet;  //ParetoSet解集

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

    void InitialPop(vector<vector<vector<int>>> &JFDTimePop, vector<vector<vector<int>>> &HierarchyPop);

    void EvolutionProcess(int mu);

    /*********************组进化************************/
    void BasedInd_RandFamInFacTobestPos_New(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                            vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                            int NadirPointMS, float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                                            vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS, float &ObjectTEC);

    void BasedInd_SwapFam_New(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                              vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                              int NadirPointMS, float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                              vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS, float &ObjectTEC);

    void BasedInd_ShiftFam_New(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam, vector<vector<int>> &JFDTime,
                               vector<vector<int>> &Hierarchy, int NadirPointMS, float NadirPointTEC, int IdealPointMS,
                               float IdealPointTEC,
                               vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC, int &ObjectMS,
                               float &ObjectTEC );

    void BasedIndLS_SetupSwap_New(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                  vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy, int NadirPointMS,
                                  float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                                  vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS, float &ObjectTEC);

    void BasedInd_DeAndConFams_New(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                   vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy, int NadirPointMS,
                                   float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                                   vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS, float &ObjectTEC);


    /**********************工件进化********************************/
    void BasedInd_CirJobInFamTobestPos_New(const vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                           vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                           int NadirPointMS, float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                                           vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS, float &ObjectTEC);
    void
    JobLocalSearch(const vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam, vector<vector<int>> &JFDTime,
                   vector<vector<int>> &Hierarchy, int NadirPointMS, float NadirPointTEC, int IdealPointMS,
                   float IdealPointTEC, vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC,
                   int &ObjectMS, float &ObjectTEC);

    void BasedInd_SwapJob_New(const vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                              vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                              int NadirPointMS, float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                              vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS, float &ObjectTEC);

    void BasedInd_JobsInFam_New(int Fac, const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,
                                vector<vector<int>> &JobSeqInFam, vector<vector<int>> &JFDTime,
                                vector<vector<int>> &Hierarchy,
                                int Fam, int NadirPointMS, float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                                vector<Individual> &CCEAPopulation, int OrgMS, float OrgTEC, int &ObjectMS,
                                float &ObjectTEC);

    /*******************************参考集进化************************************/
    void BasedInd_RefDeconAndCon_New(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                     vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                     int NadirPointMS, float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                                     vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS,
                                     float &ObjectTEC, vector<double> &scores_combine);

    void  ApplyDeconAndCon_Ref (int DeAndConIndex ,vector<vector<int>> &Sol, vector<vector<int>> &JobSeqInFam,vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                vector<int> &FamsExtracted,unordered_map<int, vector<int>> &JobsExtracted, int NadirPointMS,float NadirPointTEC, int IdealPointMS,
                                float IdealPointTEC, vector<int>& FacSpan, vector<float>& FacEC, vector<Individual> &CCEAPopulation);


    void BasedInd_Destruction_FamsAndJobs_New_Opeator0(vector<vector<int>>& FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                       vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                                       vector<int> &FamsExtracted,
                                                       unordered_map<int, vector<int>> &JobsExtracted);

    void BasedInd_Destruction_FamsAndJobs_New_Opeator1(vector<vector<int>>& FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                       vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                                       vector<int> &FamsExtracted,
                                                       unordered_map<int, vector<int>> &JobsExtracted, vector<int>& FacSpan, vector<float>& FacEC);

    void BasedInd_Destruction_FamsAndJobs_New_Opeator2(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                       vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                                       vector<int> &FamsExtracted,
                                                       unordered_map<int, vector<int>> &JobsExtracted,
                                                       vector<int> &FacSpan, vector<float> &FacEC);

    void BasedInd_Construction_FamsAndJobs_New_Opeator0(vector<vector<int>>& FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                        vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                                        vector<int> &FamsExtracted,
                                                        unordered_map<int, vector<int>> &JobsExtracted, int NadirPointMS,
                                                        float NadirPointTEC, int IdealPointMS, float IdealPointTEC ,vector<Individual> &CCEAPopulation);

    void BasedInd_Construction_FamsAndJobs_New_Opeator1(vector<vector<int>>& FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                        vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                                        vector<int> &FamsExtracted,
                                                        unordered_map<int, vector<int>> &JobsExtracted, int NadirPointMS,
                                                        float NadirPointTEC, int IdealPointMS, float IdealPointTEC ,vector<Individual> &CCEAPopulation);

    void BasedInd_Construction_FamsAndJobs_New_Opeator2(vector<vector<int>>& FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                        vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                                        vector<int> &FamsExtracted,
                                                        unordered_map<int, vector<int>> &JobsExtracted, int NadirPointMS,
                                                        float NadirPointTEC, int IdealPointMS, float IdealPointTEC ,vector<Individual> &CCEAPopulation);

    static bool SA( int OrgMS, float OrgTEC, int AfterMS, float AfterTEC, double temperature);

    /*******************************档案集进化************************************/
    void UpdateArchiveGroupJobSet(int mu);

    void SaveTECandDeMS(vector<Individual> &CCEAPopulation);
    

};


#endif //DBFGSP_NEW_QCCEANoQNOQ_H
