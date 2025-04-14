
#ifndef DBFGSP_NEW_NOPERATOR_H
#define DBFGSP_NEW_NOPERATOR_H

#pragma once

#include "Problem.h"
#include "Individual.h"
#include "Base.h"

class NOperator :
        virtual public Problem
{
public:

    NOperator();

    virtual ~NOperator();

    vector<int> MachReadyTime;

    void CheckSol(vector<int> FamSeq, vector<vector<int>> JobSeqinFam, int Span);

    void CheckSol(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam, int Span);

    void CheckSolTEC(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam, float TEC);

    int GetSpan(vector<int> FamSeq, vector<vector<int>> JobSeqinFam);

    int GetSpan(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam);

    int GetSpan(vector<vector<int>> FamSeq, vector<vector<int>> JobSeqinFam, vector<int> &FacSpan);

    int GetSpan_Forward(const vector<vector<int>> &FamSeq, const vector<vector<int>> &JobSeqInFam,vector<int> &FacSpan);

    float GetTEC(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam);

    float GetPerFacEC(vector<int> FamSeqInFac, vector<vector<int>> JobSeqInFam, vector<vector<int>> &JFDTime);

    int GetSpanPerFacByJFD(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,const vector<vector<int>> &JFDTime);

    float GetTECForAllFacByJFD(const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam, const vector<vector<int>> &JFDTime);

    float GetTECForPerFacByJFD(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam, const vector<vector<int>> &JFDTime);

    int GetSpan_Forward_InFactory(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam);

    int GetJFDTime_Forward(const vector<vector<int>> &FamSeq, const vector<vector<int>> &JobSeqInFam,vector<vector<int>> &JFDTime, vector<int> &FacSpan);

    int GetJFDTime_Forward_InFactory(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,vector<vector<int>> &JFDTime);

    int GetJBDTime_Backward(const vector<vector<int>> &FamSeq, const vector<vector<int>> &JobSeqinFam,vector<vector<int>> &JBDTime, vector<int> &FacSpan);

    int GetJBDTime_Backward_InFactory(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqinFam, vector<vector<int>> &JBDTime);

    void GetMSandTECForPerandToalFac(const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,vector<int> &FacSpan, vector<float> &FacEC, int &Makespan, float&TotalEC);

    void GetMSandTECForPerandToalFacByJFD(const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,const vector<vector<int>> &JFDTime, vector<int> &FacSpan,vector<float> &FacEC,int &Makespan, float &TotalEC);

    int GetSol_Include(vector<int> FamPrmu, vector<vector<int>> JobSeqInFam, vector<vector<int>> &FacFamSeq,vector<int> &FacSpan);

    int GetDelayTime_Forward(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam, vector<vector<int>> &JFDTime,vector<int> &FacSpan, vector<vector<int>> &DelayTime);

    int GetDelayTime_Forward_InFactory(vector<int> FamSeqInFac, vector<vector<int>> JobSeqInFam,vector<vector<int>> &JFDTime,vector<vector<int>> &DelayTime);
    // Ind
    float
    GetIndForPerFacAfterInsertFam_DR(int Fac, const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,
                                     const vector<vector<int>> &JobSeqInFam, const vector<vector<int>> &JFDTime,
                                     const vector<vector<int>> &JBDTime, int InsFam, int Pos, int nadirpointMS,
                                     float nadirpointTEC, int idealpointMS, float idealpointTEC);
    // Ind
    float
    GetIndForPerFacAfterInsertJob_DR(int fac, const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,
                                     const vector<vector<int>> &JobSeqInFam, const vector<vector<int>> &JFDTime,
                                     const vector<vector<int>> &JBDTime, int InsFam, int FamPos, int InsJob,
                                     int JobPos, int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                     float idealpointTEC);
    // Ind
    float GetIndForPerFacAfterInsertFam(int Fac, const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,const vector<vector<int>> &JobSeqInFam,
                                        const vector<vector<int>> &JFDTime,const vector<vector<int>> &JBDTime,
                                        int InsFam, int Pos, int nadirpointMS,float nadirpointTEC, int idealpointMS,
                                        float idealpointTEC, int OrgMS,int OrgTEC,vector<Individual> &CMOEAPopulation);

    int FindBestPosToInsertFam_NoAC(const vector<vector<int>> &FamSeq, const vector<vector<int>> &JobSeqinFam,int InsFam, int &bestFac, int &bestPos);

    int FindBestPosToInsertFam_NoAC_InFactory(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqinFam,int InsFam, int &bestPos);

    int FindBestPosToInsertJob_NoAC(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqinFam, int Fam,int InsJob, int &bestPos);

    int FindBestPosToInsertFam_NoAC_Span(vector<vector<int>> &FamSeq, const vector<vector<int>> &JobSeqinFam,int InsFam, int &bestFac,int &bestPos);

    int FindBestPosToInsertFam_NoAC_InFactory_Span(vector<vector<int>> &FamSeq,const vector<vector<int>> &JobSeqInFam, int InsFam, int &bestPos, int factory);

    int FindBestPosToInsertJob_NoAC_Span(const vector<vector<int>> &FamSeq, vector<vector<int>> &JobSeqinFam,vector<int> &FacTT_AllScenarios, int Fam, int InsJob,int &bestPos);

    int FindBestPosToInsertFam_InFactory(vector <int> FamSeqInFac, vector <vector <int>> JobSeqinFam,  vector<vector<int>> JFDTime,vector<vector<int>> JBDTime, int InsFam, int& bestPos);


    int GetSpanForPerFacAfterInsertFam(const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam,const vector<vector<int>>& JFDTime, const vector<vector<int>>& JBDTime, int InsFam, int Pos);

    float GetECForPerFacAfterInsertFam(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,const vector<vector<int>> &JFDTime, const vector<vector<int>> &JBDTime, int InsFam,int Pos);

    float FindBestPosToInsertFamForAllFac_TEC(const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,const vector<vector<int>> &JFDTime, const vector<vector<int>> &JBDTime,int InsFam, int &BestFac, int &BestPos);

    int FindBestPosToInsertFamForAllFac_Makespan(const vector<vector<int>> &FacFamSeq,const vector<vector<int>> &JobSeqInFam,const vector<vector<int>> &JFDTime, const vector<vector<int>> &JBDTime,int InsFam, int &BestFac, int &BestPos);

    float FindBestPosToInsertFamForPerFac_TEC(const vector<int> &NewFamSeqInFac, const vector<vector<int>> &JobSeqInFam,const vector<vector<int>> &JFDTime, const vector<vector<int>> &JBDTime,int InsFam, int &BestPos);

    int FindBestPosToInsertFamForPerFac_Makespan(const vector<int> &NewFamSeqInFac, const vector<vector<int>> &JobSeqInFam,const vector<vector<int>> &JFDTime, const vector<vector<int>> &JBDTime,int InsFam, int &BestPos);

    int FindBestPosToInsertFam(vector<vector<int>> FamSeq, vector<vector<int>> JobSeqinFam,vector<vector<int>> JFDTime, vector<vector<int>> JBDTime, int InsFam, int &bestFac, int &bestPos);


    // Ind
    float FindBestPosToInsertFamForAllFac_Ind(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& JCTime, const vector<vector<int>>& JSTime,
                                              int InsFam, int& BestFac, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, int OrgMS, int OrgTEC, vector<Individual>& CMOEAPopulation);
    // Ind

    float FindBestPosToInsertFamForPerFac_Ind(int Fac, const vector<vector<int>>& FacFamSeq, const vector<int>& NewFamSeqInFac, const vector<vector<int>>& JobSeqInFam,
                                              const vector<vector<int>>& JCTime, const vector<vector<int>>& JSTime, int InsFam, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, int OrgMS, int OrgTEC, vector<Individual>& CMOEAPopulation);

    // Ind
    float
    FindBestPosToInsertFamForAllFac_Ind_DR(const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                           const vector<vector<int>> &JFDTime, const vector<vector<int>> &JBDTime,
                                           int InsFam,
                                           int &BestFac, int &BestPos, int nadirpointMS, float nadirpointTEC,
                                           int idealpointMS, float idealpointTEC);
    // Ind
    float
    FindBestPosToInsertFamForPerFac_Ind_DR(int Fac, const vector<vector<int>> &FacFamSeq,
                                           const vector<int> &NewFamSeqInFac,
                                           const vector<vector<int>> &JobSeqInFam, const vector<vector<int>> &JFDTime,
                                           const vector<vector<int>> &JBDTime, int InsFam, int &BestPos,
                                           int nadirpointMS,
                                           float nadirpointTEC, int idealpointMS, float idealpointTEC);

    // Ind
    float GetIndForPerFacAfterInsertJob(int fac, const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,
                                        const vector<vector<int>> &JobSeqInFam, const vector<vector<int>> &JFDTime,
                                        const vector<vector<int>> &JBDTime, int InsFam, int FamPos, int InsJob,
                                        int JobPos,
                                        int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                        int OrgMS, float OrgTEC, vector<Individual> &CMOEAPopulation);
    // Ind
    void Basedind_JobsInFam(int fac, const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,
                            vector<vector<int>> &JobSeqInFam, vector<vector<int>> &JFDTime, vector<vector<int>> &JBDTime,
                            int Fam, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                            vector<Individual> &CMOEAPopulation, int OrgMS, float OrgTEC, int &ObjectMS,
                            float &ObjectTEC);

    void Pareto_relation(vector<Individual> &Population);

    void Speed_mutation(vector<Individual> &Population, vector<Individual> &Temp);

    void Speed_mutation_new(vector<Individual> &Population, vector<Individual> &Temp);

    void Speed_mutation(vector<Individual> &CMOEAPopulation, vector<Individual> &Temp,
                        vector<Individual> &tureCMOEAPopulation);

    void SwapFam(vector<int> &FamSeq);

    void SwapFam(vector<vector<int>> &FamSeq);

    void SwapFaminFac(vector<int> &FamSeqInFac);

    void SwapFamBetweenFacs(vector<vector<int>> &FamSeq);

    void SwapJob(vector<vector<int>> &JobSeqinFam);

    void JobSwap(vector<int> FamSeq, vector<vector<int>> &JobSeqinFam, int &Span);

    void JobInsert(vector<int> FamSeq, vector<vector<int>> &JobSeqinFam, int &Span);

    void FamInsert(vector<int> &FamSeq, vector<vector<int>> JobSeqinFam, int &Span);

    void IG_DR(vector<vector<int>> FamSeqinFac, vector<vector<int>> &JobSeqinFam, vector<int> FacSpan);

    void InsertFamBetweenFac(vector<vector<int>> &FamSeqinFac, vector<vector<int>> JobSeqinFam);


    void RefreshJFDJBDTime_InFactory(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqinFam,vector<vector<int>> &JFDTime, vector<vector<int>> &JBDTime, int ForwardFamPos,int BackwardFamPos, int FowardJobPos, int BackwardJobPos);

    void CopyJFDJBDTime(const vector<vector<int>>& srcJFDTime, const vector<vector<int>>& srcJBDTime,
                        vector<vector<int>>& dstJFDTime, vector<vector<int>>& dstJBDTime);

    vector<int> MergeSequences(const vector<int> &SortedSeq1, const vector<int> &SortedSeq2);

    void CombineSortedSequences(int SortMethod1, int SortMethod2, vector<int> &FamPermu);

    void SortJobsInFam(int SortMethod, vector<vector<int>> &JobSeqInFam);

    void SortFam(int SortMethod, vector<int> &FamPermu);

    float
    FindBestPosToInsertFamForAllFac_TEC(const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                        int InsFam, int &BestFac, int &BestPos);

    float FindBestPosToInsertFamForPerFac_TEC(const vector<int> &NewFamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                              int InsFam, int &BestPos);

    float GetTECForPerFac(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam);

    double
    FindBestPosToInsertFamForAllFacs_MC(const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                        int InsFam, int &BestFac, int &BestPos);

    double GetMCForPerFacAfterInsertFam(int Fac, const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,
                                        const vector<vector<int>> &JobSeqInFam, int InsFam, int Pos);

    double
    FindBestPosToInsertFam_InFactory_MC(int Fac, const vector<vector<int>> &FacFamSeq,
                                        const vector<int> &NewFamSeqInFac,
                                        const vector<vector<int>> &JobSeqInFam, int InsFam, int &BestPos);

    float
    FindBestPosToInsertFamForAllFac_Ind_NoAC(const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                             int InsFam, int &BestFac, int &BestPos, int nadirpointMS, float nadirpointTEC,
                                             int idealpointMS, float idealpointTEC, int OrgMS, int OrgTEC,
                                             vector<Individual> &CMOEAPopulation);

    float
    FindBestPosToInsertFamForPerFac_Ind_NoAC(int Fac, const vector<vector<int>> &FacFamSeq,
                                             const vector<int> &NewFamSeqInFac,
                                             const vector<vector<int>> &JobSeqInFam, int InsFam, int &BestPos,
                                             int nadirpointMS,
                                             float nadirpointTEC, int idealpointMS, float idealpointTEC, int OrgMS,
                                             int OrgTEC,
                                             vector<Individual> &CMOEAPopulation);

    float GetIndForPerFacAfterInsertFam_NoAC(int Fac, const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,
                                             const vector<vector<int>> &JobSeqInFam, int InsFam, int Pos, int nadirpointMS,
                                             float nadirpointTEC, int idealpointMS, float idealpointTEC, int OrgMS,
                                             int OrgTEC,
                                             vector<Individual> &CMOEAPopulation);

    float GetTECForAllFac(const vector<vector<int>> &FamSeq, const vector<vector<int>> &JobSeqInFam);

    // Ind
    float
    FindBestPosToInsertFamForAllFac_Ind_DR_NoAC(const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                                int InsFam,
                                                int &BestFac, int &BestPos, int nadirpointMS, float nadirpointTEC,
                                                int idealpointMS, float idealpointTEC,vector<Individual> &CCEAPopulation);
    // Ind
    float
    FindBestPosToInsertFamForPerFac_Ind_DR_NoAC(int Fac, const vector<vector<int>> &FacFamSeq,
                                                const vector<int> &NewFamSeqInFac,
                                                const vector<vector<int>> &JobSeqInFam,int InsFam, int &BestPos,
                                                int nadirpointMS,
                                                float nadirpointTEC, int idealpointMS, float idealpointTEC,vector<Individual> &CCEAPopulation);

    float
    GetIndForPerFacAfterInsertFam_DR_NoAC(int Fac, const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,
                                          const vector<vector<int>> &JobSeqInFam, int InsFam, int Pos, int nadirpointMS,
                                          float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual> &CCEAPopulation);

    int GetDelayTime_Forward_InFactory_NoAC(vector<int> FamSeqInFac, vector<vector<int>> JobSeqInFam,
                                            vector<vector<int>> &DelayTime);

    int GetDelayTime_Forward_NoAC(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam, vector<int> &FacSpan,
                                  vector<vector<int>> &DelayTime);

    float
    FindBestPosToInsertJobForPerFac_Ind_NoAC(int fac, const vector<vector<int>> &FacFamSeq,
                                             const vector<int> &FamSeqInFac,
                                             const vector<vector<int>> JobSeqInFam, int Fam, int FamPos, int InsJob,
                                             int &BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                             float idealpointTEC, int OrgMS, float OrgTEC,
                                             vector<Individual> &CCEAPopulation);

    float
    GetIndForPerFacAfterInsertJob_NoAC(int fac, const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,
                                       const vector<vector<int>> &JobSeqInFam, int InsFam, int FamPos, int InsJob,
                                       int JobPos, int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                       float idealpointTEC, int OrgMS, float OrgTEC,
                                       vector<Individual> &CCEAPopulation);

    void Speed_mutation_NoAC(vector<Individual> &QCCEAPopulation, vector<Individual> &Temp,
                             vector<Individual> &tureQCCEAPopulation);

    float
    GetIndForPerFacAfterInsertJob_DR_NoAC(int fac, const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,
                                          const vector<vector<int>> &JobSeqInFam, int InsFam, int FamPos, int InsJob,
                                          int JobPos, int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                          float idealpointTEC, vector<Individual> &CCEAPopulation);
};

#endif //DBFGSP_NEW_NOPERATOR_H