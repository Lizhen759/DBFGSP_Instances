#pragma once
#include "NOperator.h"
#include <unordered_map>
#include "Individual.h"
#include "Heu.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include "Base.h"
#include <valarray>
#include <sstream>
#include <numeric>

class CMOEA :
        public Heu, public Individual
{
public:
    CMOEA();
    virtual ~CMOEA();
    int RunEvolution(int CPUTime, vector<vector<Individual>>& CMOEAFinalAfterRepParetoSet, int AN, int mu);

private:
    long m_InitTime;
    long m_TimeLimit;
    int m_FamD;
    int m_BlockLengthChoice;
    int m_NoImproveNumber;
    int m_Temperature;

    //asdf
    int m_RefSize;
    int m_Popsize;

    float thr_zeta;

    int m_nadirpointMS;    //最低点
    float m_nadirpointTEC;

    int m_idealpointMS;	  //理想点
    float m_idealpointTEC;

    vector<vector<vector<int>>> m_RefSpeedVector;
    vector<vector<vector<int>>> m_RefFacFamSeqArray;
    vector<vector<vector<int>>> m_RefJobSeqinFamArray;
    vector<vector<int>> m_RefFacSpanArray;
    vector<vector<float>> m_RefFacECArray;
    vector<int> m_RefSpanArray;
    vector<float> m_RefTECArray;
    vector <int> m_nRefCriFacArray;

    vector<Individual> CMOEAPopulation;
    vector<Individual> CMOEAParetoSet;  //ParetoSet解集

    vector<vector<vector<int>>> m_FacFamSeqArray;
    vector<vector<vector<int>>> m_JobSeqinFamArray;
    vector<int> m_SpanArray1;
    vector<int> m_SpanArray2;
    vector<float> m_TECArray1;
    vector<float> m_TECArray2;
    vector <int> m_nCriFacArray1;
    vector <int> m_nCriFacArray2;

    int m_PS1;
    int m_PS2;
    vector<int> m_Map1;
    vector<int> m_Map2;
    vector<bool> m_bFlag1;// 为参考解组序列更新 for the reference familiy sequence updated
    vector<bool> m_bFlag2;// 为参考解工件序列更新 for the reference job sequence updated
    vector<int> m_Age1;
    vector<int> m_Age2;
    vector<int> m_RecordSpan;

    void SetParameters(int BlockLengthChoice, int NoImproveNumber, int NumberOfExtractedFams, long TimeLimit, int AllJobTotalPTime);

    void InitialPop(vector<vector<vector<int>>>& JFDTimePop, vector<vector<vector<int>>>& JBDTimePop);

    void EvolutionProcess(int mu);

    void BasedindRandFamInFacTobestPos(vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqinFam, vector<vector<int>>& JFDTime, vector<vector<int>>& JBDTime, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual>& CMOEAPopulation, int& ObjectMS, float& ObjectTEC);

    void BasedindSwapFam(vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, vector<vector<int>>& JFDTime, vector<vector<int>>& JBDTime,
                         int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual>& CMOEAPopulation, int& ObjectMS, float& ObjectTEC);

    void BasedindCirJobInFamTobestPos(const vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, vector<vector<int>>& JFDTime, vector<vector<int>>& JBDTime,
                                      int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual>& CMOEAPopulation, int& ObjectMS, float& ObjectTEC);

    void BasedindSwapJob(const vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, vector<vector<int>>& JFDTime, vector<vector<int>>& JBDTime,
                         int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual>& CMOEAPopulation, int& ObjectMS, float& ObjectTEC);

    void BasedindRefInsert(vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, vector<vector<int>>& JFDTime, vector<vector<int>>& JBDTime,
                           int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual>& CMOEAPopulation, int& ObjectMS, float& ObjectTEC);

    void BasedindDestruction_FamsAndJobs(vector<vector<int>>& Sol, vector<vector<int>>& JobSeqInFam, vector<vector<int>>& JFDTime, vector<vector<int>>& JBDTime,
                                         vector<int>& FamsExtracted, unordered_map<int, vector<int>>& JobsExtracted, int FamsExtractedMethods);

    void Basedind_Construction_FamsAndJobs(vector<vector<int>>& Sol, vector<vector<int>>& JobSeqInFam, vector<vector<int>>& JFDTime, vector<vector<int>>& JBDTime,
                                           vector<int>& FamsExtracted, unordered_map<int, vector<int>>& JobsExtracted, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC);

    void UpdateArchiveGroupJobSet(int mu);

    void SaveTECandDeMS(vector<Individual>& CMOEAPopulation);
};
