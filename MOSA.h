#pragma once
#include "Heu.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include "Individual.h"
#include "Base.h"
#include <valarray>
#include <sstream>
class MOSA :
        public Heu, public Individual
{
public:
    MOSA();
    virtual ~MOSA();
    void RunEvolution(int CPUFactor, int Rep, vector<vector<Individual>>& MOSAFinalAfterRepParetoSet);

private:

    int iteration_number;
    double T0;
    double TF;
    double cool_ratio;
    //xin
    Individual MOSASol;
    vector<Individual> MOSAPopulation;
    vector<Individual> MOSAParetoSet;  //ParetoSet解集
    int m_nadirpointMS;    //最低点
    float m_nadirpointTEC;

    int m_idealpointMS;	  //理想点
    float m_idealpointTEC;

    vector<float> m_RefTECArray;
    vector<vector<vector<int>>> m_RefSpeedVector; //参考集速度矩阵
    vector<vector<int>> m_RefFacSpanArray;
    vector<vector<vector<int>>> m_RefFacFamSeqArray;
    vector<vector<int>> m_RefFamSeqArray;
    vector<vector<vector<int>>> m_RefJobSeqinFamArray;
    vector<int> m_RefSpanArray;
    int m_RefSize;

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
    vector<int> m_RecordSpan;



private:
    void InitPop();

    void SetParameters(int RefSize, int PS1, int PS2, int AgeLimit, int DesLen, long TimeLimit);
    //int Destruction_Construction(int Len, vector<int>& FamSeq, vector<vector<int>> JobSeqinFam, int& Span);
    //重载
    //void MOSAReInitPop();
    int Destruction_Construction(int Len, vector<vector<int>>& FacFamSeq, vector<vector<int>> JobSeqinFam, vector<int>& FacSpan);
    void MOSAEvolution();
    void NGP1(vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqinFam, vector<vector<int>>& m_SpeedMatrix);
    void NGP2(vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqinFam, vector<vector<int>>& m_SpeedMatrix);
    void NGP3(vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqinFam, vector<vector<int>>& m_SpeedMatrix);
    void LS_strategy(vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqinFam, vector<vector<int>>& m_SpeedMatrix);

};

