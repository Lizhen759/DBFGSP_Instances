#ifndef DBFGSP_NEW_NOPERATOR_NEW_H
#define DBFGSP_NEW_NOPERATOR_NEW_H

#include "Problem.h"
#include "Individual.h"
#include <iostream>
#include <numeric>
#include <valarray>
#include <algorithm>
#include <iomanip>
#include "Base.h"

using namespace Base;

class NOperator_New : public Problem {

public:

    NOperator_New();

    virtual ~NOperator_New();

    vector<int> MachReadyTime;

    void Pareto_relation(vector<Individual> &Population);

    //Heu
    void SortJobsInFam(int SortMethod, vector<vector<int>> &JobSeqInFam);

    void CombineSortedSequences(int SortMethod1, int SortMethod2, vector<int> &FamPermu);

    vector<int> MergeSequences(const vector<int> &SortedSeq1, const vector<int> &SortedSeq2);

    void SortFam(int SortMethod, vector<int> &FamPermu);

    //NOperator
    void CheckSol(vector<int> FamSeq, vector<vector<int>> JobSeqinFam, int Span);

    void CheckSol(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam, int Span);


    void CheckSolTEC(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam, float TEC);

    float GetTEC(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam);

    float GetPerFacEC(vector<int> FamSeqInFac, vector<vector<int>> JobSeqInFam, vector<vector<int>> &JFDTime);

    int GetSpan(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam);

    int GetSpan(vector<int> FamSeq, vector<vector<int>> JobSeqinFam);

    int GetSpanPerFacByJFD(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                           const vector<vector<int>> &JFDTime);

    int GetSpanForAllFacByJFD(const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                              const vector<vector<int>> &JFDTime);
    int
    GetDelayTime_Forward_InFactory(vector<int> FamSeqInFac, vector<vector<int>> JobSeqinFam,
                                   vector<vector<int>> &JFDTime,vector<vector<int>> &DelayTime);

    int
    GetDelayTime_Forward(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam, vector<vector<int>> &JFDTime,
                         vector<int> &FacSpan, vector<vector<int>> &DelayTime);


    int GetJFDTime_Forward(const vector<vector<int>> &FamSeq, const vector<vector<int>> &JobSeqInFam,
                           vector<vector<int>> &JFDTime, vector<int> &FacSpan);

    int GetJFDTime_Forward_InFactory(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                     vector<vector<int>> &JFDTime);

    float GetTECForPerFacByJFD(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                               const vector<vector<int>> &JFDTime);

    float GetTECForAllFacByJFD(const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                               const vector<vector<int>> &JFDTime);


    void GetMSandTECForPerandToalFacByJFD(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam,
                                          const vector<vector<int>>& JFDTime, vector<int>& FacSpan, vector<float>& FacEC, int& Makespan, float& TotalEC);

    double GetMCForPerFacAfterInsertFam_New(int Fac , const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                            const vector<vector<int>> &JFDTime, const vector<vector<int>> &Hierarchy,
                                            int InsFam, int Pos);


    void GetJFDHierarchy_Forward_New(const vector<vector<int>> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                    vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy);

    void GetJFDHierarchy_Forward_InFactory_New(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                               vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy);

    int GetSpanForPerFacAfterInsertFam_New(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                           const vector<vector<int>> &JFDTime, const vector<vector<int>> &hierarchy,
                                           int InsFam, int Pos);

    int GetSpanForPerFacAfterInsertJob_New(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                           const vector<vector<int>> &JFDTime,const vector<vector<int>> &Hierarchy,int InsFam, int FamPos, int InsJob,
                                           int JobPos);

    float GetTECForPerFacAfterInsertFam_New(const vector<int> &FamSeqInFac, const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                            const vector<vector<int>> &JFDTime, const vector<vector<int>> &Hierarchy,
                                            int InsFam, int Pos);

    float GetTECForPerFacAfterInsertJob_New(const vector<int> &FamSeqInFac, const vector<vector<int>> &FacFamSeq,
                                             vector<vector<int>> &JobSeqInFam, const vector<vector<int>> &JFDTime,
                                            const vector<vector<int>> &hierarchy, int InsFam, int FamPos, int InsJob,
                                            int JobPos);
    float
    GetIndForPerFacAfterInsertFam_New(int Fac, const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,
                                      const vector<vector<int>> &JobSeqInFam, const vector<vector<int>> &JFDTime,
                                      const vector<vector<int>> &Hierarchy, int InsFam, int Pos, int nadirpointMS,
                                      float nadirpointTEC, int idealpointMS, float idealpointTEC,vector<int>& FacSpan, vector<float>& FacEC, int OrgMS,
                                      float OrgTEC, vector<Individual> &CCEAPopulation);

    float
    GetIndForPerFacAfterInsertJob_New(int fac, const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,
                                      const vector<vector<int>> &JobSeqInFam, const vector<vector<int>> &JFDTime,
                                      const vector<vector<int>> &Hierarchy, int InsFam, int FamPos, int InsJob,
                                      int JobPos, int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                      float idealpointTEC, int OrgMS, float OrgTEC,
                                      vector<Individual> &CCEAPopulation);

    float
    GetIndForPerFacAfterInsertJob_DR_New(int fac, const vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,
                                         const vector<vector<int>> &JobSeqInFam, const vector<vector<int>> &JFDTime,
                                         const vector<vector<int>> &Hierarchy, int InsFam, int FamPos, int InsJob,
                                         int JobPos, int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                         float idealpointTEC    ,vector<Individual> &CCEAPopulation);


    float
    GetIndForPerFacAfterInsertFam_DR_New(int Fac,   const  vector<vector<int>> &FacFamSeq, const vector<int> &FamSeqInFac,
                                         const vector<vector<int>> &JobSeqInFam, const vector<vector<int>> &JFDTime,
                                         const vector<vector<int>> &Hierarchy, int InsFam, int Pos, int nadirpointMS,
                                         float nadirpointTEC, int idealpointMS, float idealpointTEC ,vector<Individual> &CCEAPopulation);

    double
    FindBestPosToInsertFamForAllFacs_MC_New(const vector<vector<int>> &FacFamSeq,
                                            const vector<vector<int>> &JobSeqInFam,
                                            const vector<vector<int>> &JFDTime, const vector<vector<int>> &Hierarchy,
                                            int InsFam, int &BestFac, int &BestPos);

    double
    FindBestPosToInsertFam_InFactory_MC_New(int Fac , const vector<vector<int>> &FacFamSeq,const vector<int> &NewFamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                            const vector<vector<int>> &JFDTime, const vector<vector<int>> &Hierarchy,
                                            int InsFam, int &BestPos);


    float
    FindBestPosToInsertFam_InFactory_TEC_New(const vector<int> &NewFamSeqInFac, const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                             const vector<vector<int>> &JFDTime, const vector<vector<int>> &Hierarchy,
                                             int InsFam, int &BestPos);


    float
    FindBestPosInsertFamAllFacs_TEC_New(const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                        const vector<vector<int>> &JFDTime, const vector<vector<int>> &Hierarchy,
                                        int InsFam, int &BestFac, int &BestPos);


    int
    FindBestPosToInsertFam_InFactory_MakeSpan_New(const vector<int> &NewFamSeqInFac,
                                                  const vector<vector<int>> &JobSeqInFam,
                                                  const vector<vector<int>> &JFDTime,
                                                  const vector<vector<int>> &Hierarchy,
                                                  int InsFam, int &BestPos);

    int FindBestPosToInsertFamForAllFacs_MakeSpan_New(const vector<vector<int>> &FacFamSeq,
                                                      const vector<vector<int>> &JobSeqInFam,
                                                      const vector<vector<int>> &JFDTime,
                                                      const vector<vector<int>> &Hierarchy, int InsFam, int &BestFac,
                                                      int &BestPos);

    float FindBestPosToInsertJobForPerFac_Ind_New(int fac ,const vector<vector<int>> &FacFamSeq,const vector<int>& FamSeqInFac, const vector<vector<int>> JobSeqInFam,
                                                   const vector<vector<int>>& JFDTime, const vector<vector<int>>& Hierarchy, int Fam, int FamPos, int InsJob,
                                                   int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                                   float idealpointTEC, int OrgMS, float OrgTEC,
                                                   vector<Individual> &CCEAPopulation);
    float FindBestPosToInsertFamForPerFac_Ind_New(int Fac, const vector<vector<int>> &FacFamSeq,
                                                  const vector<int> &NewFamSeqInFac,
                                                  const vector<vector<int>> &JobSeqInFam,
                                                  const vector<vector<int>> &JFDTime,
                                                  const vector<vector<int>> &Hierarchy,
                                                  int InsFam, int &BestPos, int nadirpointMS, float nadirpointTEC,
                                                  int idealpointMS, float idealpointTEC, vector<int> &FacSpan,
                                                  vector<float> &FacEC, int OrgMS, float OrgTEC,
                                                  vector<Individual> &CCEAPopulation);
    float
    FindBestPosToInsertFamForAllFac_Ind_New(const vector<vector<int>> &FacFamSeq,
                                            const vector<vector<int>> &JobSeqInFam,
                                            const vector<vector<int>> &JFDTime, const vector<vector<int>> &Hierarchy,
                                            int InsFam, int &BestFac, int &BestPos, int nadirpointMS,
                                            float nadirpointTEC,
                                            int idealpointMS, float idealpointTEC, vector<int> &FacSpan,
                                            vector<float> &FacEC, int OrgMS,float OrgTEC,
                                            vector<Individual> &CCEAPopulation);

    float
    FindBestPosToInsertFamForAllFac_Ind_DR_New(const vector<vector<int>> &FacFamSeq,
                                               const vector<vector<int>> &JobSeqInFam,
                                               const vector<vector<int>> &JFDTime, const vector<vector<int>> &Hierarchy,
                                               int InsFam, int &BestFac, int &BestPos, int nadirpointMS,
                                               float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual> &CCEAPopulation);

    float FindBestPosToInsertFamForPerFac_Ind_DR_New(int Fac, const vector<vector<int>> &FacFamSeq,
                                                     const vector<int> &NewFamSeqInFac,
                                                     const vector<vector<int>> &JobSeqInFam,
                                                     const vector<vector<int>> &JFDTime,
                                                     const vector<vector<int>> &Hierarchy, int InsFam, int &BestPos,
                                                     int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                                     float idealpointTEC, vector<Individual> &CCEAPopulation);

    void scheduling_one_job(int job, vector<int> &machine_ready);

    void scheduling_one_job(int job, vector<int> &machine_ready, vector<int> &one_job_hierarchy);


    bool scheduling_one_job(int job, vector<int> &machine_ready, const vector<int> &one_job_hierarchy);

    bool scheduling_one_job_with_Update_Return(int job, vector<int> &machine_ready, vector<int> &one_job_hierarchy);

    void
    RefreshJFDTimeHierarchy_InFactory_Insert(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                             vector<vector<int>> &JFDTime, vector<vector<int>> &hierarchy,
                                             int StartFamPos, int StartJobPos, int EndJobPos);

    void RefreshJFDTimeHierarchy_InFactory_Erase(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                                 vector<vector<int>> &JFDTime, vector<vector<int>> &hierarchy,
                                                 int StartFamPos, int StartJobPos);

    void CopyJFDHierarchy(const vector<vector<int>> &srcJFDTime, const vector<vector<int>> &srcHierarchy,
                          vector<vector<int>> &dstJFDTime, vector<vector<int>> &dstHierarchy);

    void
    Speed_mutation(vector<Individual> &CCEAPopulation, vector<Individual> &Temp,
                   vector<Individual> &tureCCEAPopulation);


};


#endif //DBFGSP_NEW_NOPERATOR_NEW_H
