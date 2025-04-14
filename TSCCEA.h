#pragma once
#include "Heu.h"
#include "Individual.h"
#include "Base.h"
#include "Problem.h"
#include "NOperator.h"
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <sstream>
#include <algorithm>
class TSCCEA :
        public Heu, public Individual
{
public:
    TSCCEA();
    virtual ~TSCCEA();
    void RunEvolution(int CPUFactor, int Rep, vector<vector<Individual>>& TSCCEAFinalAfterRepParetoSet);

private:
    //xin
    vector<Individual> TSCCEAPopulation;
    vector<Individual> TSCCEAParetoSet;  //ParetoSet解集
    int m_nadirpointMS;    //最低点
    float m_nadirpointTEC;

    int m_idealpointMS;	  //理想点
    float m_idealpointTEC;

    int m_RefSize;
    vector<int> m_RefSpanArray; //档案集
    vector<float> m_RefTECArray;
    vector<vector<int>> m_RefFacSpanArray;
    vector<vector<vector<int>>> m_RefSpeedVector; //参考集速度矩阵
    vector<vector<vector<int>>> m_RefFacFamSeqArray;
    vector<vector<int>> m_RefFamSeqArray;
    vector<vector<vector<int>>> m_RefJobSeqinFamArray;

    vector<vector<vector<int>>> m_FacFamSeqArray;
    vector<vector<int>> m_FamSeqArray;
    vector<vector<vector<int>>> m_JobSeqinFamArray;

    vector<int> m_SpanArray1;
    vector<int> m_SpanArray2;
    vector<float> m_TECArray1;
    vector<float> m_TECArray2;

    int m_DesLen;
    int m_PS1;
    int m_PS2;
    vector<int> m_Map1;
    vector<int> m_Map2;
    vector<bool> m_bFlag1;// for the reference familiy sequence updated
    vector<bool> m_bFlag2;// for the reference job sequence updated
    vector<int> m_Age1;
    vector<int> m_Age2;
    int m_AgeLimit;
    long m_InitTime;
    long m_TimeLimit;



private:
    void InitPop();

    void ReInitPop();

    void SetParameters(int RefSize, int PS1, int PS2, int AgeLimit, int DesLen, long TimeLimit);

    int Destruction_Construction(int Len, vector<vector<int>>& FacFamSeq, vector<vector<int>> JobSeqinFam, vector<int>& FacSpan);

    void Evolution();
};

