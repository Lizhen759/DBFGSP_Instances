#include "NOperator.h"
#include <algorithm>
#include <iostream>
#include "Base.h"
#pragma once
using namespace Base;
class Heu :
        virtual public NOperator
{
public:
    Heu();

    virtual ~Heu();

protected:
    void SortJobsinFam(int Factor, int SortMethod, vector<vector<int>> &JobSeqinFam); //Method=0:LPT; Method=1:SPT

    void SortJobsinOtherFam(int Factor,int SortMethod,vector<vector<int>> &JobSeqinFam);

    void SortFam(int Factor, int SortMethod, vector<int> &FamPrmu);

    void SortFam(int SortMethod, vector<int> &FamPermu);

    void SortOtherFam(int Factor, int SortMethod, vector<int> &FamSeq);

    int NEHFam(vector<int> FamPrmuInFac, vector<vector<int>> JobSeqinFam, vector<int> &FamSeqInFac, int &SpanInFac);

    int NEHFam(vector<int> FamPrmu, vector<vector<int>> JobSeqinFam, vector<vector<int>> &FacFamSeq, vector<int> &FacSpan);

    void JPA_TS(vector<vector<int>> &FacFamSeq);

    void JPA_G(vector<vector<int>> &FacFamSeq, vector<vector<int>> &InitialJobinFamSeq);

    void setDifference(vector<int> Originseq, vector<int> Removeseq, vector<int> &Restseq);

    void LRX(vector<vector<int>> FamPop, vector<vector<int>> &JobSeqinFamPop);

    void schedulejob(vector<int> &jobseq, vector<int> &MachReadyTime);

    void SortJobsInFam(int SortMethod, vector<vector<int>> &JobSeqInFam);


};

