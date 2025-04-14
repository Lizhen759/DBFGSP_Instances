#include "QCCEANoRapid.h"

#include <random>

QCCEANoRapid::QCCEANoRapid() {

}

QCCEANoRapid::~QCCEANoRapid() {
}

void
QCCEANoRapid::SetParameters(int AN, long TimeLimit, double T) {

    this->m_TimeLimit = TimeLimit;

    this->m_PopSize = 1000;
    this->m_RefSize = AN;
    this->m_PS1 = AN;
    this->m_PS2 = AN;
    this->thr_zeta = 1.0;

    this->m_T = T;

    m_QCCEANoRapidPopulation.clear();
    m_QCCEANoRapidPopulation.resize(m_RefSize);

    m_NadirPointMS = -1;
    m_NadirPointTEC = -1;

    m_IdealPointMS = -1;
    m_IdealPointTEC = -1;

    //档案集 Reference set
    m_ArchiveSpanArray.clear();//档案集完工时间
    m_ArchiveSpanArray.resize(m_RefSize);

    m_ArchiveTECArray.clear();//档案集总能耗
    m_ArchiveTECArray.resize(m_RefSize);

    m_ArchiveSpeedVector.clear();
    m_ArchiveSpeedVector.resize(m_RefSize);

    m_ArchiveFacFamSeqArray.clear();//档案集工厂组
    m_ArchiveFacFamSeqArray.resize(m_RefSize);

    m_ArchiveJobSeqinFamArray.clear();//档案集组工件
    m_ArchiveJobSeqinFamArray.resize(m_RefSize);

    m_ArchiveFacSpanArray.clear();//工厂完工时间
    m_ArchiveFacSpanArray.resize(m_RefSize);

    m_ArchiveFacECArray.clear();//工厂能耗
    m_ArchiveFacECArray.resize(m_RefSize);

    //组种群
    m_SpanArray1.clear();//完工时间
    m_SpanArray1.resize(m_PS1);

    m_TECArray1.clear(); //能耗
    m_TECArray1.resize(m_PS1);

    m_Map1.clear();//协作者在档案集的下标
    m_Map1.resize(m_PS1);

    m_FacFamSeqArray.clear();//组序列
    m_FacFamSeqArray.resize(m_PS1);

    //工件种群
    m_SpanArray2.clear();//完工时间
    m_SpanArray2.resize(m_PS2);

    m_TECArray2.clear(); //能耗
    m_TECArray2.resize(m_PS2);


    m_Map2.clear();//协作者在档案集的下标
    m_Map2.resize(m_PS2);


    m_JobSeqinFamArray.clear();//工件序列
    m_JobSeqinFamArray.resize(m_PS2);
}


void QCCEANoRapid::BasedInd_RandFamInFacTobestPos(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                               int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                               vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS,
                                               float &ObjectTEC) {


    int CriFac = 0;
    int max = INT_MIN;
    //关键工厂
    for (int i = 0; i < FacFamSeq.size(); i++) {
        if (max < FacSpan[i]) {
            max = FacSpan[i];
            CriFac = i;
        }
    }
    int BestFac = -1, BestPos = -1;
    float Ind = 0.0;

    vector<vector<int>> TempFacFamSeq;
    TempFacFamSeq.clear();
    TempFacFamSeq = FacFamSeq;

    for (int fam = 0; fam < TempFacFamSeq[CriFac].size(); fam++) {
        if (FacFamSeq[CriFac].size() < 2)
            continue;

        int CurFam = TempFacFamSeq[CriFac][fam];

        auto it = find(FacFamSeq[CriFac].begin(), FacFamSeq[CriFac].end(), CurFam); //定义迭代器
        int j = it - FacFamSeq[CriFac].begin();
        FacFamSeq[CriFac].erase(it);

        BestFac = -1;
        BestPos = -1;
        Ind = this->FindBestPosToInsertFamForAllFac_Ind_NoAC(FacFamSeq, JobSeqInFam,CurFam, BestFac, BestPos, idealpointMS, nadirpointTEC,
                                                          idealpointMS, idealpointTEC, ObjectMS, ObjectTEC ,CCEAPopulation);
        FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);


        vector<int> AfterInsertSpanFac(FacFamSeq.size(), 0);
        vector<float> AfterInsertECFac(FacFamSeq.size(), 0);
        int AfterInsertMS = -1;
        float AfterInsertTEC = -1;
        GetMSandTECForPerandToalFac(FacFamSeq, JobSeqInFam, AfterInsertSpanFac, AfterInsertECFac, AfterInsertMS, AfterInsertTEC);
        ObjectMS = AfterInsertMS;
        ObjectTEC = AfterInsertTEC;
    }
}

void QCCEANoRapid::BasedInd_SwapFam(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,int nadirpointMS,
                                 float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                 vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC,
                                 int &ObjectMS, float &ObjectTEC) {

    int OrgMS = ObjectMS;
    float OrgTEC = ObjectTEC;

    //判断是否替换最低点和理想点
    if (OrgMS > nadirpointMS)
        nadirpointMS = OrgMS;
    if (OrgMS < idealpointMS)
        idealpointMS = OrgMS;

    //判断是否替换最低点和理想点
    if (OrgTEC > nadirpointTEC)
        nadirpointTEC = OrgTEC;
    if (OrgTEC < idealpointTEC)
        idealpointTEC = OrgTEC;

    //归一化
    float OrgnormalMS = -1;
    float OrgnormalTEC = -1;

    OrgnormalMS = (static_cast<float>(OrgMS - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    OrgnormalTEC = ((OrgTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));

    //计算交换之前指标Ind
    float OrgConvergence_ind = 0;
    OrgConvergence_ind += (OrgnormalMS - 1.0f) * (OrgnormalMS - 1.0f);
    OrgConvergence_ind += (OrgnormalMS - 1.0f) * (OrgnormalTEC - 1.0f);
    OrgConvergence_ind = sqrt(OrgConvergence_ind);
    OrgConvergence_ind = 1 / (OrgConvergence_ind + 1);


    float  OrgDistribution_ind = 0;
    Individual::calc_distribution_ind(CCEAPopulation);
    for (int i = 0; i < CCEAPopulation.size(); i++) {
        for (int j = i + 1; j < CCEAPopulation.size(); j++) {
            OrgDistribution_ind += CCEAPopulation[i].distribution_ind[j];
        }
    }
    OrgDistribution_ind = OrgDistribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

    float OrgCombined_ind= 0.5f*OrgConvergence_ind+0.5f*OrgDistribution_ind;

    for (int i = 0; i < m_Families / m_Factories; i++) {

        int CriFac = -1;
        int max = -1;
        int OptFac = -1;
        int min = INT_MAX;
        //关键工厂和makespan最小的工厂
        for (int i = 0; i < FacFamSeq.size(); i++) {
            if (max < FacSpan[i]) {
                max = FacSpan[i];
                CriFac = i;
            }
            if (min > FacSpan[i]) {
                min = FacSpan[i];
                OptFac = i;
            }
        }
        // 确保关键工厂和最优工厂不相同
        if (CriFac == OptFac) {
            break;
        }

        if (FacFamSeq[CriFac].empty()) {
            break;
        }



        int Old_Cri_Span = FacSpan[CriFac];
        float Old_Cri_EC = FacEC[CriFac];
        int Old_Opt_Span = FacSpan[OptFac];
        float Old_Opt_EC = FacEC[OptFac];

        int Pos1 = wyt_rand(FacFamSeq[CriFac].size());
        int Fam1 = FacFamSeq[CriFac][Pos1];
        FacFamSeq[CriFac].erase(FacFamSeq[CriFac].begin() + Pos1);
        if (FacFamSeq[OptFac].empty()) {
            break;
        }
        int Pos2 = rand() % FacFamSeq[OptFac].size();
        int Fam2 = FacFamSeq[OptFac][Pos2];
        FacFamSeq[OptFac].erase(FacFamSeq[OptFac].begin() + Pos2);

        FacFamSeq[CriFac].insert(FacFamSeq[CriFac].begin() + Pos1, Fam2);

        FacSpan[CriFac] = GetSpan_Forward_InFactory(FacFamSeq[CriFac], JobSeqInFam);
        FacEC[CriFac] = GetTECForPerFac(FacFamSeq[CriFac], JobSeqInFam);
        FacFamSeq[OptFac].insert(FacFamSeq[OptFac].begin() + Pos2, Fam1);
        FacSpan[OptFac] =GetSpan_Forward_InFactory(FacFamSeq[OptFac], JobSeqInFam);
        FacEC[OptFac] = GetTECForPerFac(FacFamSeq[OptFac], JobSeqInFam);

        // 计算交换后的 Makespan/TEC
        int AfterMakespan = *max_element(FacSpan.begin(), FacSpan.end());
        float AfterTEC = accumulate(FacEC.begin(), FacEC.end(), 0.0f);

        //判断是否替换最低点和理想点
        if (AfterMakespan > nadirpointMS)
            nadirpointMS = AfterMakespan;
        if (AfterMakespan < idealpointMS)
            idealpointMS = AfterMakespan;

        //判断是否替换最低点和理想点
        if (AfterTEC > nadirpointTEC)
            nadirpointTEC = AfterTEC;
        if (AfterTEC < idealpointTEC)
            idealpointTEC = AfterTEC;

        //判断是否比操作前改进，若改进则加入参考集
        bool flag = true;
        if ((AfterMakespan < OrgMS && AfterTEC <= OrgTEC) || (AfterMakespan <= OrgMS && AfterTEC < OrgTEC)) {
            for (int i = 0; i < CCEAPopulation.size(); i++) {
                if ((AfterMakespan == CCEAPopulation[i].MS) && (AfterTEC == CCEAPopulation[i].TEC)) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                Individual tempIndi;
                tempIndi.m_FacFamSeqArray = FacFamSeq;
                tempIndi.m_JobSeqInFamArray = JobSeqInFam;
                tempIndi.MS = AfterMakespan;
                tempIndi.TEC = AfterTEC;
                tempIndi.m_SpeedVector = m_SpeedMatrix;  //速度矩阵
                CCEAPopulation.push_back(tempIndi);
            }
        }

        //归一化
        float AfternormalMS = -1;
        float AfternormalTEC = -1;


        AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) /static_cast<float>(nadirpointMS - idealpointMS));
        AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));

        //计算指标Ind
        float Afterconvergence_ind = 0;
        Afterconvergence_ind += (AfternormalMS - 1.0f) * (AfternormalMS - 1.0f);
        Afterconvergence_ind += (AfternormalTEC - 1.0f) * (AfternormalTEC - 1.0f);
        Afterconvergence_ind = sqrt(Afterconvergence_ind);
        Afterconvergence_ind = 1 / (Afterconvergence_ind + 1);

        float  AfterDistribution_ind = 0;
        Individual::calc_distribution_ind(CCEAPopulation);
        for (int i = 0; i < CCEAPopulation.size(); i++) {
            for (int j = i + 1; j < CCEAPopulation.size(); j++) {
                AfterDistribution_ind += CCEAPopulation[i].distribution_ind[j];
            }
        }
        AfterDistribution_ind = AfterDistribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

        float AfterCombined_ind = 0.5f * AfterDistribution_ind + 0.5f * Afterconvergence_ind;

        if (AfterCombined_ind < OrgCombined_ind) {
            ObjectMS = AfterMakespan;
            ObjectTEC = AfterTEC;
        } else {
            FacFamSeq[OptFac].erase(FacFamSeq[OptFac].begin() + Pos2);
            FacFamSeq[CriFac].erase(FacFamSeq[CriFac].begin() + Pos1);
            FacFamSeq[OptFac].insert(FacFamSeq[OptFac].begin() + Pos2, Fam2);
            FacFamSeq[CriFac].insert(FacFamSeq[CriFac].begin() + Pos1, Fam1);
            FacSpan[CriFac]= Old_Cri_Span ;
            FacEC[CriFac] = Old_Cri_EC ;
            FacSpan[OptFac] = Old_Opt_Span;
            FacEC[OptFac] = Old_Opt_EC ;
        }
    }

}

void QCCEANoRapid::BasedIndLS_SetupSwap(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                     int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                     vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC,
                                     int &ObjectMS, float &ObjectTEC) {


    int OrgMS = ObjectMS;
    float OrgTEC = ObjectTEC;

    //判断是否替换最低点和理想点
    if (OrgMS > nadirpointMS)
        nadirpointMS = OrgMS;
    if (OrgMS < idealpointMS)
        idealpointMS = OrgMS;

    //判断是否替换最低点和理想点
    if (OrgTEC > nadirpointTEC)
        nadirpointTEC = OrgTEC;
    if (OrgTEC < idealpointTEC)
        idealpointTEC = OrgTEC;

    //归一化
    float OrgnormalMS = -1;
    float OrgnormalTEC = -1;

    OrgnormalMS = (static_cast<float>(OrgMS - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    OrgnormalTEC = ((OrgTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));

    //计算交换之前指标Ind
    float OrgConvergence_ind = 0;
    OrgConvergence_ind += (OrgnormalMS - 1.0f) * (OrgnormalMS - 1.0f);
    OrgConvergence_ind += (OrgnormalTEC - 1.0f) * (OrgnormalTEC - 1.0f);
    OrgConvergence_ind = sqrt(OrgConvergence_ind);
    OrgConvergence_ind = 1 / (OrgConvergence_ind + 1);
    //cout << "个体的收敛指标：" << OrgConvergence_ind << endl;

    float  OrgDistribution_ind = 0;
    Individual::calc_distribution_ind(CCEAPopulation);
    for (int i = 0; i < CCEAPopulation.size(); i++) {
        for (int j = i + 1; j < CCEAPopulation.size(); j++) {
            OrgDistribution_ind += CCEAPopulation[i].distribution_ind[j];
        }
    }
    OrgDistribution_ind = OrgDistribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

    float OrgCombined_ind= 0.5f*OrgConvergence_ind+0.5f*OrgDistribution_ind;


    for (int i = 0; i < m_Families / m_Factories; i++) {
        int CriMS = 0;
        int OptMS = INT_MAX;
        int CriFac = 0, OptFac = 0;
        //关键工厂
        for (int i = 0; i < FacFamSeq.size(); i++) {
            if (CriMS < FacSpan[i]) {
                CriFac = i;  // 关键工厂
                CriMS = FacSpan[i];  // 更新最大值
            }
            if (OptMS >= FacSpan[i])  // 排除关键工厂
            {
                OptFac = i;  // 最优工厂
                OptMS = FacSpan[i];  // 更新最小值
            }
        }

        int CurFam;
        int PreFam;
        int NowFam1;
        int NowFam2;
        int Max_Setup = 0;
        int Pos1 = -1;
        int Pos2 = -1;
        int Min_Setup = INT_MAX;

        PreFam = FacFamSeq[CriFac][0];
        for (int Fam = 0; Fam < FacFamSeq[CriFac].size(); Fam++) {
            CurFam = FacFamSeq[CriFac][Fam];
            if (this->m_SetupTime[0][PreFam][CurFam] > Max_Setup) {
                Max_Setup = this->m_SetupTime[0][PreFam][CurFam];
                Pos1 = Fam;
                NowFam1 = CurFam;
            }
            PreFam = CurFam;
        }
        // 如果该工厂的组序列为空，直接跳过
        if (FacFamSeq[OptFac].empty())
        {
            continue;  // 或者 return; 具体看逻辑是否需要继续后续操作
        }
        PreFam = FacFamSeq[OptFac][0];
        for (int Fam = 0; Fam < FacFamSeq[OptFac].size(); Fam++) {
            CurFam = FacFamSeq[OptFac][Fam];
            if (this->m_SetupTime[0][PreFam][CurFam] < Min_Setup) {
                Min_Setup = this->m_SetupTime[0][PreFam][CurFam];
                Pos2 = Fam;
                NowFam2 = CurFam;
            }
            PreFam = CurFam;
        }


        int Old_Cri_Span = FacSpan[CriFac];
        float Old_Cri_EC = FacEC[CriFac];
        int Old_Opt_Span = FacSpan[OptFac];
        float Old_Opt_EC = FacEC[OptFac];

        FacFamSeq[CriFac].erase(FacFamSeq[CriFac].begin() + Pos1);
        FacFamSeq[OptFac].erase(FacFamSeq[OptFac].begin() + Pos2);


        //交换，在crifac工厂插入组
        FacFamSeq[CriFac].insert(FacFamSeq[CriFac].begin() + Pos1, NowFam2);
        FacSpan[CriFac] = GetSpan_Forward_InFactory(FacFamSeq[CriFac], JobSeqInFam);
        FacEC[CriFac] = GetTECForPerFac(FacFamSeq[CriFac], JobSeqInFam);

        //交换，在optfac工厂插入组
        FacFamSeq[OptFac].insert(FacFamSeq[OptFac].begin() + Pos2, NowFam1);

        FacSpan[OptFac] =  GetSpan_Forward_InFactory(FacFamSeq[OptFac], JobSeqInFam);
        FacEC[OptFac] =  GetTECForPerFac(FacFamSeq[OptFac], JobSeqInFam);

        // 计算交换后的 Makespan/TEC
        int AfterMakespan = *max_element(FacSpan.begin(), FacSpan.end());
        float AfterTEC = accumulate(FacEC.begin(), FacEC.end(), 0.0f);
        //判断是否替换最低点和理想点
        if (AfterMakespan > nadirpointMS)
            nadirpointMS = AfterMakespan;
        if (AfterMakespan < idealpointMS)
            idealpointMS = AfterMakespan;
        if (AfterTEC > nadirpointTEC)
            nadirpointTEC = AfterTEC;
        if (AfterTEC < idealpointTEC)
            idealpointTEC = AfterTEC;

        //判断是否比操作前改进，若改进则加入参考集
        bool flag = true;
        if (((AfterMakespan < OrgMS) && (AfterTEC < OrgTEC)) || ((AfterMakespan < OrgMS) && (AfterTEC == OrgTEC)) ||
            ((AfterMakespan == OrgMS) && (AfterTEC < OrgTEC)) || (AfterMakespan < idealpointMS) ||
            (AfterTEC < idealpointTEC)) {
            for (int i = 0; i < CCEAPopulation.size(); i++) {
                if ((AfterMakespan == CCEAPopulation[i].MS) && (AfterTEC == CCEAPopulation[i].TEC)) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                //cout << "*********************" << endl;
                Individual tempIndi;
                tempIndi.m_FacFamSeqArray = FacFamSeq;
                tempIndi.m_JobSeqInFamArray = JobSeqInFam;
                tempIndi.MS = AfterMakespan;
                tempIndi.TEC = AfterTEC;
                tempIndi.m_SpeedVector = m_SpeedMatrix;  //速度矩阵
                CCEAPopulation.push_back(tempIndi);
            }
        }

        //归一化
        float AfternormalMS = -1;
        float AfternormalTEC = -1;
        //cout << endl << "normalize" << endl;

        AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) /
                         static_cast<float>(nadirpointMS - idealpointMS));
        AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
        //cout << "个体归一化后的MS：" << AfternormalMS << "\tTEC：" << AfternormalTEC << endl;

        //计算指标Ind
        float Afterconvergence_ind = 0;
        Afterconvergence_ind += (AfternormalMS - 1.0f) * (AfternormalMS - 1.0f);
        Afterconvergence_ind += (AfternormalTEC - 1.0f) * (AfternormalTEC - 1.0f);
        Afterconvergence_ind = sqrt(Afterconvergence_ind);
        Afterconvergence_ind = 1 / (Afterconvergence_ind + 1);

        float  AfterDistribution_ind = 0;
        Individual::calc_distribution_ind(CCEAPopulation);
        for (int i = 0; i < CCEAPopulation.size(); i++) {
            for (int j = i + 1; j < CCEAPopulation.size(); j++) {
                AfterDistribution_ind += CCEAPopulation[i].distribution_ind[j];
            }
        }
        AfterDistribution_ind = AfterDistribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

        float AfterCombined_ind = 0.5f * AfterDistribution_ind + 0.5f * Afterconvergence_ind;


        if (AfterCombined_ind < OrgCombined_ind) {
            ObjectMS = AfterMakespan;
            ObjectTEC = AfterTEC;
        } else {
            FacFamSeq[CriFac].erase(FacFamSeq[CriFac].begin() + Pos1);
            FacFamSeq[OptFac].erase(FacFamSeq[OptFac].begin() + Pos2);
            FacFamSeq[CriFac].insert(FacFamSeq[CriFac].begin() + Pos1, NowFam1);
            FacFamSeq[OptFac].insert(FacFamSeq[OptFac].begin() + Pos2, NowFam2);
            FacSpan[CriFac]= Old_Cri_Span ;
            FacEC[CriFac] = Old_Cri_EC ;
            FacSpan[OptFac] = Old_Opt_Span;
            FacEC[OptFac] = Old_Opt_EC ;
        }
    }

}

void QCCEANoRapid::BasedInd_ShiftFam(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                  int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                  vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC,
                                  int &ObjectMS, float &ObjectTEC) {
    int OrgMS = ObjectMS;
    float OrgTEC = ObjectTEC;

    //判断是否替换最低点和理想点
    if (OrgMS > nadirpointMS)
        nadirpointMS = OrgMS;
    if (OrgMS < idealpointMS)
        idealpointMS = OrgMS;

    //判断是否替换最低点和理想点
    if (OrgTEC > nadirpointTEC)
        nadirpointTEC = OrgTEC;
    if (OrgTEC < idealpointTEC)
        idealpointTEC = OrgTEC;

    //归一化
    float OrgnormalMS = -1;
    float OrgnormalTEC = -1;

    OrgnormalMS = (static_cast<float>(OrgMS - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    OrgnormalTEC = ((OrgTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));

    //计算交换之前指标Ind
    float OrgConvergence_ind = 0;
    OrgConvergence_ind += (OrgnormalMS - 1.0f) * (OrgnormalMS - 1.0f);
    OrgConvergence_ind += (OrgnormalTEC - 1.0f) * (OrgnormalTEC - 1.0f);
    OrgConvergence_ind = sqrt(OrgConvergence_ind);
    OrgConvergence_ind = 1 / (OrgConvergence_ind + 1);
    //cout << "个体的收敛指标：" << OrgConvergence_ind << endl;

    float  OrgDistribution_ind = 0;
    Individual::calc_distribution_ind(CCEAPopulation);
    for (int i = 0; i < CCEAPopulation.size(); i++) {
        for (int j = i + 1; j < CCEAPopulation.size(); j++) {
            OrgDistribution_ind += CCEAPopulation[i].distribution_ind[j];
        }
    }
    OrgDistribution_ind = OrgDistribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

    float OrgCombined_ind= 0.5f*OrgConvergence_ind+0.5f*OrgDistribution_ind;
    int facIndex1 = -1;
    int facIndex2 = -1;

    if (facIndex1 == -1 && facIndex2 == -1) {
        do {
            facIndex1 = wyt_rand(FacFamSeq.size());
            facIndex2 = wyt_rand(FacFamSeq.size());
        } while (facIndex1 == facIndex2 || FacFamSeq[facIndex1].size() < 2 || FacFamSeq[facIndex2].empty());
    }


    int pos = wyt_rand(FacFamSeq[facIndex1].size());


    int Old_Span1 =FacSpan[facIndex1];
    float Old_EC1 = FacEC[facIndex1] ;
    int Old_Span2=FacSpan[facIndex2] ;
    float Old_EC2 =FacEC[facIndex2]  ;

    int CurFam = FacFamSeq[facIndex1][pos];
    FacFamSeq[facIndex1].erase(FacFamSeq[facIndex1].begin() + pos);

    FacSpan[facIndex1] = GetSpan_Forward_InFactory(FacFamSeq[facIndex1], JobSeqInFam);
    FacEC[facIndex1] =  GetTECForPerFac(FacFamSeq[facIndex1], JobSeqInFam);


    int bestPos = -1;
    float Ind = FindBestPosToInsertFamForPerFac_Ind_NoAC(facIndex2, FacFamSeq, FacFamSeq[facIndex2], JobSeqInFam, CurFam, bestPos,
                                                        nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, OrgMS,
                                                        OrgTEC, CCEAPopulation);

    FacFamSeq[facIndex2].insert(FacFamSeq[facIndex2].begin() + bestPos, CurFam);
    FacSpan[facIndex2] = GetSpan_Forward_InFactory(FacFamSeq[facIndex2], JobSeqInFam);
    FacEC[facIndex2] = GetTECForPerFac(FacFamSeq[facIndex2], JobSeqInFam);

    int AfterMakespan = *max_element(FacSpan.begin(), FacSpan.end());
    float AfterTEC = accumulate(FacEC.begin(), FacEC.end(), 0.0f);

    if (AfterMakespan > nadirpointMS)
        nadirpointMS = AfterMakespan;
    if (AfterMakespan < idealpointMS)
        idealpointMS = AfterMakespan;
    if (AfterTEC > nadirpointTEC)
        nadirpointTEC = AfterTEC;
    if (AfterTEC < idealpointTEC)
        idealpointTEC = AfterTEC;

    //归一化
    float AfternormalMS = -1;
    float AfternormalTEC = -1;

    AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) /static_cast<float>(nadirpointMS - idealpointMS));
    AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));

    //计算指标Ind
    float Afterconvergence_ind = 0;
    Afterconvergence_ind += (AfternormalMS - 1.0f) * (AfternormalMS - 1.0f);
    Afterconvergence_ind += (AfternormalTEC - 1.0f) * (AfternormalTEC - 1.0f);
    Afterconvergence_ind = sqrt(Afterconvergence_ind);
    Afterconvergence_ind = 1 / (Afterconvergence_ind + 1);

    float  AfterDistribution_ind = 0;
    Individual::calc_distribution_ind(CCEAPopulation);
    for (int i = 0; i < CCEAPopulation.size(); i++) {
        for (int j = i + 1; j < CCEAPopulation.size(); j++) {
            AfterDistribution_ind += CCEAPopulation[i].distribution_ind[j];
        }
    }
    AfterDistribution_ind = AfterDistribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

    float AfterCombined_ind = 0.5f * AfterDistribution_ind + 0.5f * Afterconvergence_ind;

    if (AfterCombined_ind < OrgCombined_ind) {
        ObjectMS = AfterMakespan;
        ObjectTEC = AfterTEC;
    } else {
        // 结果没有改善，需要回滚
        FacFamSeq[facIndex2].erase(FacFamSeq[facIndex2].begin() + bestPos);
        FacFamSeq[facIndex1].insert(FacFamSeq[facIndex1].begin() + pos, CurFam);
        FacSpan[facIndex1]= Old_Span1 ;
        FacEC[facIndex1] = Old_EC1 ;
        FacSpan[facIndex2]= Old_Span2 ;
        FacEC[facIndex2] = Old_EC2 ;
    }

}


void QCCEANoRapid::BasedInd_DeAndConFams(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                      int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                      vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC,
                                      int &ObjectMS, float &ObjectTEC) {

    int CriFac;
    float max = FLT_MIN;

    for (int i = 0; i < FacFamSeq.size(); i++) {
        if (max < FacEC[i]) {
            max = FacEC[i];
            CriFac = i;
        }
    }

    int Len = FacFamSeq[CriFac].size();
    vector<int> ExtractFamSeq;
    bool groupSelected = false;
    do {
        int Fac;
        if (ExtractFamSeq.size() < Len / 2) //从关键工厂挑选d/2个组
            Fac = CriFac;
        else  //随机工厂挑选剩下的d/2个组
            Fac = rand() % m_Factories;
        if (FacFamSeq[Fac].size() > 1) {
            int Pos = rand() % FacFamSeq[Fac].size();
            ExtractFamSeq.push_back(FacFamSeq[Fac][Pos]);  //将挑选出的组添加到ExtractFamSeq
            FacFamSeq[Fac].erase(FacFamSeq[Fac].begin() + Pos);  //从原序列中删除pos位置的组
            groupSelected = true;
        }

    } while (ExtractFamSeq.size() < Len && groupSelected);

    // 将挑选出来的d个组插入到最好的位置 Insert the extracted Families into the best Positions
    for (int Fam = 0; Fam < ExtractFamSeq.size(); Fam++) {

        int CurFam = ExtractFamSeq[Fam];  //当前组
        int bestFac = -1, bestPos = -1;  //最好的工厂，位置
        if (Fam == ExtractFamSeq.size() - 1) {
            float Ind = this->FindBestPosToInsertFamForAllFac_Ind_DR_NoAC(FacFamSeq, JobSeqInFam,
                                                                         CurFam, bestFac, bestPos, nadirpointMS,
                                                                         nadirpointTEC, idealpointMS, idealpointTEC,CCEAPopulation);
        } else {
            float randValue = wyt_rand_include_right(0.0f,1.0f);
            if (randValue < 0.5) {

                FindBestPosToInsertFam_NoAC_Span(FacFamSeq, JobSeqInFam, CurFam,
                                                              bestFac, bestPos);

            } else {
                FindBestPosToInsertFamForAllFac_TEC(FacFamSeq, JobSeqInFam,  CurFam, bestFac,
                                                    bestPos);

            }
        }
        // 在最好的工厂的最好位置中插入当前组 Insert CurFam to bestPos at bestFac
        FacFamSeq[bestFac].insert(FacFamSeq[bestFac].begin() + bestPos, CurFam);
    }

}


void QCCEANoRapid::BasedInd_SwapJob(const vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                 int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                 vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC,
                                 int &ObjectMS, float &ObjectTEC) {

    int OrgMS = ObjectMS;
    float OrgTEC = ObjectTEC;

    //判断是否替换最低点和理想点
    if (OrgMS > nadirpointMS)
        nadirpointMS = OrgMS;
    if (OrgMS < idealpointMS)
        idealpointMS = OrgMS;

    //判断是否替换最低点和理想点
    if (OrgTEC > nadirpointTEC)
        nadirpointTEC = OrgTEC;
    if (OrgTEC < idealpointTEC)
        idealpointTEC = OrgTEC;

    //归一化
    float OrgnormalMS = -1;
    float OrgnormalTEC = -1;
    //cout << endl << "normalize" << endl;

    OrgnormalMS = (static_cast<float>(OrgMS - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    OrgnormalTEC = ((OrgTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
    //cout << "原始个体归一化后的MS：" << OrgnormalMS << "\tTEC：" << OrgnormalTEC << endl;

    //计算交换之前指标Ind
    float OrgConvergence_ind = 0;
    OrgConvergence_ind += (OrgnormalMS - 1.0f) * (OrgnormalMS - 1.0f);
    OrgConvergence_ind += (OrgnormalTEC - 1.0f) * (OrgnormalTEC - 1.0f);
    OrgConvergence_ind = sqrt(OrgConvergence_ind);
    OrgConvergence_ind = 1 / (OrgConvergence_ind + 1);
    //cout << "个体的收敛指标：" << OrgConvergence_ind << endl;

    float  OrgDistribution_ind = 0;
    Individual::calc_distribution_ind(CCEAPopulation);
    for (int i = 0; i < CCEAPopulation.size(); i++) {
        for (int j = i + 1; j < CCEAPopulation.size(); j++) {
            OrgDistribution_ind += CCEAPopulation[i].distribution_ind[j];
        }
    }
    OrgDistribution_ind = OrgDistribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

    float OrgCombined_ind= 0.5f*OrgConvergence_ind+0.5f*OrgDistribution_ind;


    int fac = -1;
    int min = -1;

    //关键工厂
    for (int i = 0; i < FacFamSeq.size(); i++) {
        if (min < FacSpan[i]) {
            min = FacSpan[i];
            fac = i;
        }

    }

    for (int Fam = 0; Fam < FacFamSeq[fac].size(); Fam++) {
        int CurFam = FacFamSeq[fac][Fam]; //当前组
        if (JobSeqInFam[CurFam].size() > 1) //组工件数大于1，进行局部搜索
        {
            int pt1 = 0;
            int pt2 = 0;

            int originalFacSpan = FacSpan[fac];

            if (JobSeqInFam[CurFam].size() > 2) {

                pt1 = rand() % JobSeqInFam[CurFam].size();
                do {
                    pt2 = rand() % JobSeqInFam[CurFam].size();
                } while (pt1 == pt2);

                if (pt1 > pt2) {
                    int t = pt1;
                    pt1 = pt2;
                    pt2 = t;
                }
            } else if (JobSeqInFam[CurFam].size() == 2) {
                pt1 = 0;
                pt2 = 1;
            }

            vector<vector<int>> TempJobSeqInFam = JobSeqInFam;

            auto it = find(begin(FacFamSeq[fac]), end(FacFamSeq[fac]), CurFam);
            if (it != end(FacFamSeq[fac])) {
                int famPos = it - begin(FacFamSeq[fac]);
                int Job = JobSeqInFam[CurFam][pt1];
                JobSeqInFam[CurFam][pt1] = JobSeqInFam[CurFam][pt2];
                JobSeqInFam[CurFam][pt2] = Job;


            }


            FacSpan[fac] = GetSpan_Forward_InFactory(FacFamSeq[fac], JobSeqInFam);

            int AfterMakespan = *max_element(FacSpan.begin(), FacSpan.end());

            //判断是否替换最低点和理想点
            if (AfterMakespan > nadirpointMS)
                nadirpointMS = AfterMakespan;
            if (AfterMakespan < idealpointMS)
                idealpointMS = AfterMakespan;


            float AfterTEC = 0;
            FacEC[fac] = GetTECForPerFac(FacFamSeq[fac], JobSeqInFam);
            //cout << "关键工厂" << fac << "的能耗：" << FacEC[fac] << endl;

            AfterTEC = accumulate(FacEC.begin(), FacEC.end(), 0.0f);
            //cout << "总能耗：" << AfterTEC << endl;

            //判断是否替换最低点和理想点
            if (AfterTEC > nadirpointTEC)
                nadirpointTEC = AfterTEC;
            if (AfterTEC < idealpointTEC)
                idealpointTEC = AfterTEC;

            //判断是否比操作前改进，若改进则加入参考集
            bool flag = true;
            if (((AfterMakespan < OrgMS) && (AfterTEC < OrgTEC)) || ((AfterMakespan < OrgMS) && (AfterTEC == OrgTEC)) ||
                ((AfterMakespan == OrgMS) && (AfterTEC < OrgTEC)) || (AfterMakespan < idealpointMS) ||
                (AfterTEC < idealpointTEC)) {
                for (int i = 0; i < CCEAPopulation.size(); i++) {
                    if ((AfterMakespan == CCEAPopulation[i].MS) && (AfterTEC == CCEAPopulation[i].TEC)) {
                        flag = false;
                        break;
                    }
                }
                if (flag) {
                    //cout << "*********************" << endl;
                    Individual tempIndi;
                    tempIndi.m_FacFamSeqArray = FacFamSeq;
                    tempIndi.m_JobSeqInFamArray = JobSeqInFam;
                    tempIndi.MS = AfterMakespan;
                    tempIndi.TEC = AfterTEC;
                    tempIndi.m_SpeedVector = m_SpeedMatrix;  //速度矩阵
                    CCEAPopulation.push_back(tempIndi);
                }
            }

            //归一化
            float AfternormalMS = -1;
            float AfternormalTEC = -1;
            //cout << endl << "normalize" << endl;

            AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) /
                             static_cast<float>(nadirpointMS - idealpointMS));
            AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
            //cout << "个体归一化后的MS：" << AfternormalMS << "\tTEC：" << AfternormalTEC << endl;

            //计算指标Ind
            float Afterconvergence_ind = 0;
            Afterconvergence_ind += (AfternormalMS - 1.0f) * (AfternormalMS - 1.0f);
            Afterconvergence_ind += (AfternormalTEC - 1.0f) * (AfternormalTEC - 1.0f);
            Afterconvergence_ind = sqrt(Afterconvergence_ind);
            Afterconvergence_ind = 1 / (Afterconvergence_ind + 1);
            //cout << "交换后的个体的收敛指标：" << Afterconvergence_ind << endl;

            float  AfterDistribution_ind = 0;
            Individual::calc_distribution_ind(CCEAPopulation);
            for (int i = 0; i < CCEAPopulation.size(); i++) {
                for (int j = i + 1; j < CCEAPopulation.size(); j++) {
                    AfterDistribution_ind += CCEAPopulation[i].distribution_ind[j];
                }
            }
            AfterDistribution_ind = AfterDistribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

            float AfterCombined_ind = 0.5f * AfterDistribution_ind + 0.5f * Afterconvergence_ind;


            if (AfterCombined_ind < OrgCombined_ind ) {
                ObjectMS = AfterMakespan;
                ObjectTEC = AfterTEC;
            } else {
                JobSeqInFam = TempJobSeqInFam;
                FacSpan[fac] = originalFacSpan;
            }
        }
        //检查函数
        //CheckSol(FacFamSeq, JobSeqInFam, ObjectMS);
        //CheckSolTEC(FacFamSeq, JobSeqInFam, ObjectTEC);
    }
}



void QCCEANoRapid::JobLocalSearch(const vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                           int NadirPointMS, float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                           vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS, float &ObjectTEC)
{

    for (int FacIndex = 0; FacIndex < m_Factories; ++FacIndex)
    {
        for (int FamIndex = 0; FamIndex < FacFamSeq[FacIndex].size(); ++FamIndex)
        {
            int CurFam = FacFamSeq[FacIndex][FamIndex];
            //如果组中只有一个工件，则继续
            if (JobSeqInFam[CurFam].size() == 1)
            {
                continue;
            }

            vector<int> SeqForExtracting = JobSeqInFam[CurFam];//工件序列
            shuffle(SeqForExtracting.begin(), SeqForExtracting.end(), rand_generator());//随机打乱序列
            int Counter = 0;
            int Pos = 0;
            while (Counter < ceil(SeqForExtracting.size() / 2))
            {

                vector<int> tempFacSpan(this->m_Factories, 0);
                vector<int> tempFacEC(this->m_Factories, 0);
                int CurJob = SeqForExtracting[Pos];
                auto it = find(JobSeqInFam[CurFam].begin(), JobSeqInFam[CurFam].end(), CurJob);
                int OrgPos = distance(JobSeqInFam[CurFam].begin(), it);

                JobSeqInFam[CurFam].erase(it);

                int bestPos = -1;
                auto re = FindBestPosToInsertJobForPerFac_Ind_NoAC(FacIndex,FacFamSeq,FacFamSeq[FacIndex], JobSeqInFam,CurFam,FamIndex, CurJob, bestPos, NadirPointMS,
                                                                  NadirPointTEC, IdealPointMS,  IdealPointTEC,ObjectMS, ObjectTEC, CCEAPopulation );

                JobSeqInFam[CurFam].insert(JobSeqInFam[CurFam].begin() + bestPos, CurJob);

                FacSpan[FacIndex] = GetSpan_Forward_InFactory(FacFamSeq[FacIndex], JobSeqInFam);
                FacEC[FacIndex] = GetTECForPerFac(FacFamSeq[FacIndex], JobSeqInFam);

                Pos = (Pos + 1) % (int) ceil(JobSeqInFam[CurFam].size() / 2);
                Counter += 1;
            }
        }
    }
}

void QCCEANoRapid::BasedInd_RefDeconAndCon_New(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                        int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                        vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC,
                                        int &ObjectMS, float &ObjectTEC, vector<double> &scores_combined) {

    int OrgMS = ObjectMS;
    float OrgTEC = ObjectTEC;

    //判断是否替换最低点和理想点
    if (OrgMS > nadirpointMS)
        nadirpointMS = OrgMS;
    if (OrgMS < idealpointMS)
        idealpointMS = OrgMS;

    //判断是否替换最低点和理想点
    if (OrgTEC > nadirpointTEC)
        nadirpointTEC = OrgTEC;
    if (OrgTEC < idealpointTEC)
        idealpointTEC = OrgTEC;

    // 归一化
    float OrgnormalMS = (static_cast<float>(OrgMS - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    float OrgnormalTEC = ((OrgTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));

    float OrgConvergence_ind = 0;
    OrgConvergence_ind += (OrgnormalMS - 1.0f) * (OrgnormalMS - 1.0f);
    OrgConvergence_ind += (OrgnormalTEC - 1.0f) * (OrgnormalTEC - 1.0f);
    OrgConvergence_ind = sqrt(OrgConvergence_ind);
    OrgConvergence_ind = 1 / (OrgConvergence_ind + 1);

    float  OrgDistribution_ind = 0;
    Individual::calc_distribution_ind(CCEAPopulation);
    for (int i = 0; i < CCEAPopulation.size(); i++) {
        for (int j = i + 1; j < CCEAPopulation.size(); j++) {
            OrgDistribution_ind += CCEAPopulation[i].distribution_ind[j];
        }
    }
    OrgDistribution_ind = OrgDistribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

    float OrgCombined_ind= 0.5f*OrgConvergence_ind+0.5f*OrgDistribution_ind;

    vector<vector<int>> TempFacFamSeq = FacFamSeq;
    vector<vector<int>> TempJobSeqInFam = JobSeqInFam;

    vector<int> FamsExtracted;
    unordered_map<int, vector<int>> JobsExtracted;

    // 轮盘赌选择破坏+重构组合
    auto select_operator = [](const vector<double> &weights) -> int {
        float total_weight = accumulate(weights.begin(), weights.end(), 0.0f);
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, total_weight);
        float r = dis(gen);
        float cumulative = 0.0;
        for (size_t i = 0; i < weights.size(); ++i) {
            cumulative += weights[i];
            if (r <= cumulative) {
                return i;
            }
        }
        return 0;
    };

    // 选择破坏重构组合索引
    int DeconConIndex = select_operator(scores_combined);

    // 执行对应的破坏重构操作
    ApplyDeconAndCon_Ref(DeconConIndex, FacFamSeq, JobSeqInFam, FamsExtracted,
                         JobsExtracted,nadirpointMS, nadirpointTEC,
                         idealpointMS, idealpointTEC, FacSpan, FacEC,CCEAPopulation);

    // 计算新目标值
    int AfterMakespan = GetSpan(FacFamSeq, JobSeqInFam);
    float AfterTEC = GetTECForAllFac(FacFamSeq, JobSeqInFam);

    //判断是否替换最低点和理想点
    if (AfterTEC > nadirpointTEC)
        nadirpointTEC = AfterTEC;
    if (AfterTEC < idealpointTEC)
        idealpointTEC = AfterTEC;

    //判断新解是否优于当前解，若改进则加入参考集
    bool flag = true;
    if (((AfterMakespan < OrgMS) && (AfterTEC < OrgTEC)) ||
        ((AfterMakespan < OrgMS) && (AfterTEC == OrgTEC)) ||
        ((AfterMakespan == OrgMS) && (AfterTEC < OrgTEC)) ||
        (AfterMakespan < idealpointMS) || (AfterTEC < idealpointTEC)) {

        // 检查新解是否已存在于参考集
        bool flag = true;
        for (const auto &individual: CCEAPopulation) {
            if (AfterMakespan == individual.MS && AfterTEC == individual.TEC) {
                flag = false;
                break;
            }
        }

        if (flag) {
            // 新解未在参考集中，添加到参考集
            Individual tempIndi;
            tempIndi.m_FacFamSeqArray = FacFamSeq;
            tempIndi.m_JobSeqInFamArray = JobSeqInFam;
            tempIndi.MS = AfterMakespan;
            tempIndi.TEC = AfterTEC;
            tempIndi.m_SpeedVector = m_SpeedMatrix;  // 速度矩阵
            CCEAPopulation.push_back(tempIndi);
            // 阶梯式加分
            scores_combined[DeconConIndex] += 1.5; // 破坏算子加分
            scores_combined[DeconConIndex] += 1.5;  // 重构算子加分
        } else {
            scores_combined[DeconConIndex] += 1.0; // 破坏算子加分
            scores_combined[DeconConIndex] += 1.0;  // 重构算子加分
        }
    } else {
        // 新解未改进，执行其他操作（如模拟退火的接受准则）
        double Temperature = this->m_T * accumulate(this->m_JobTotalPTime.begin(), this->m_JobTotalPTime.end(), 0) /
                             (this->m_Factories * this->m_Jobs * this->m_Machines);
        bool accept = SA(OrgMS, OrgTEC, AfterMakespan, AfterTEC, Temperature);
        if (accept) {
            // 接受较差解，给予较低的加分
            scores_combined[DeconConIndex] += 0.5; // 破坏算子加分
            scores_combined[DeconConIndex] += 0.5;  // 重构算子加分
        } else {
            // 不接受较差解，不加分
            scores_combined[DeconConIndex] += 0.0; // 破坏算子加分
            scores_combined[DeconConIndex] += 0.0;  // 重构算子加分
        }
    }

    //归一化
    float AfternormalMS = -1;
    float AfternormalTEC = -1;

    AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) /static_cast<float>(nadirpointMS - idealpointMS));
    AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));


    //计算指标Ind
    float Afterconvergence_ind = 0;
    Afterconvergence_ind += (AfternormalMS - 1.0f) * (AfternormalMS - 1.0f);
    Afterconvergence_ind += (AfternormalTEC - 1.0f) * (AfternormalTEC - 1.0f);
    Afterconvergence_ind = sqrt(Afterconvergence_ind);
    Afterconvergence_ind = 1 / (Afterconvergence_ind + 1);

    float  AfterDistribution_ind = 0;
    Individual::calc_distribution_ind(CCEAPopulation);
    for (int i = 0; i < CCEAPopulation.size(); i++) {
        for (int j = i + 1; j < CCEAPopulation.size(); j++) {
            AfterDistribution_ind += CCEAPopulation[i].distribution_ind[j];
        }
    }
    AfterDistribution_ind = AfterDistribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

    float AfterCombined_ind = 0.5f * AfterDistribution_ind + 0.5f * Afterconvergence_ind;

    if (AfterCombined_ind < OrgCombined_ind) {
        ObjectMS = AfterMakespan;
        ObjectTEC = AfterTEC;
    } else {
        FacFamSeq = TempFacFamSeq;
        JobSeqInFam = TempJobSeqInFam;
    }

}

// temperature：当前温度
// 返回true表示接受新解，返回false表示不接受
bool QCCEANoRapid::SA(int OrgMS, float OrgTEC, int AfterMS, float AfterTEC, double Temperature) {

    // Metropolis 准则
    int deltaMS = OrgMS - AfterMS;
    float deltaTEC = OrgTEC - AfterTEC;
    float delta = static_cast<float>(deltaMS) + deltaTEC;
    double r1 = wyt_rand_include_right(0.1f, 1.0f);
    if (r1 < exp(delta) / Temperature) {
        return true;
    }
    return false;
}

void QCCEANoRapid::ApplyDeconAndCon_Ref(int DeAndConIndex, vector<vector<int>> &FacFamSeq,
                                 vector<vector<int>> &JobSeqInFam,
                                 vector<int> &FamsExtracted, unordered_map<int, vector<int>> &JobsExtracted,
                                 int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                 vector<int> &FacSpan, vector<float> &FacEC,vector<Individual> &CCEAPopulation) {
    // 根据选择的算子执行破坏和重构操作
    if (DeAndConIndex == 0) {
//        cout<<"00"<<endl;
        BasedInd_Destruction_FamsAndJobs_New_Opeator0(FacFamSeq, JobSeqInFam,  FamsExtracted,
                                                      JobsExtracted);
        BasedInd_Construction_FamsAndJobs_New_Opeator0(FacFamSeq, JobSeqInFam, FamsExtracted,
                                                       JobsExtracted,
                                                       nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC,CCEAPopulation);
//        cout<<"00 end "<<endl;
    } else if (DeAndConIndex== 1) {
//        cout<<"11"<<endl;
        BasedInd_Destruction_FamsAndJobs_New_Opeator1(FacFamSeq, JobSeqInFam,  FamsExtracted,
                                                      JobsExtracted, FacSpan, FacEC);
        BasedInd_Construction_FamsAndJobs_New_Opeator1(FacFamSeq, JobSeqInFam,  FamsExtracted,
                                                       JobsExtracted,
                                                       nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC,CCEAPopulation);
//        cout<<"11 end "<<endl;
    } else {
//        cout<<"22"<<endl;
        BasedInd_Destruction_FamsAndJobs_New_Opeator2(FacFamSeq, JobSeqInFam,  FamsExtracted,
                                                      JobsExtracted, FacSpan, FacEC);
        BasedInd_Construction_FamsAndJobs_New_Opeator2(FacFamSeq, JobSeqInFam, FamsExtracted,
                                                       JobsExtracted,
                                                       nadirpointMS, nadirpointTEC, idealpointMS,
                                                       idealpointTEC,CCEAPopulation);
//        cout<<"22 end"<<endl;
    }
}

void
QCCEANoRapid::BasedInd_Destruction_FamsAndJobs_New_Opeator0(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                     vector<int> &FamsExtracted,
                                                     unordered_map<int, vector<int>> &JobsExtracted) {
    unordered_map<int, pair<int, int>> FamPosErasedFromFac;
    auto UpdateFamErasedFromFac = [&FamPosErasedFromFac](int Pos, int fac) {
        if (FamPosErasedFromFac.find(fac) == end(FamPosErasedFromFac)) {
            FamPosErasedFromFac[fac] = {Pos, Pos - 1};
        } else {
            if (Pos < FamPosErasedFromFac[fac].first) {
                FamPosErasedFromFac[fac].first = Pos;
            }
        }
    };

    FamsExtracted.clear();
    JobsExtracted.clear();
    int D_base = 2;
    if (m_Families == 40) {
        D_base = 3;
    } else if (m_Families == 60) {
        D_base = 4;
    }
    int D_variable = wyt_rand(1, 8);
    int D_extract = wyt_rand(D_base, D_base + D_variable);

    // 先检查是否有至少 2 组的工厂
    std::vector<int> availableFactories;
    for (int i = 0; i < this->m_Factories; ++i) {
        if (FacFamSeq[i].size() > 1) {
            availableFactories.push_back(i);  // 只存储至少 2 组的工厂
        }
    }

    // 如果没有符合条件的工厂，直接返回
    if (availableFactories.empty()) return;

    // 随机选择一个合适的工厂
    int Fac = availableFactories[rand() % availableFactories.size()];

    // 限制 D_extract 不超过该工厂的组数
    D_extract = std::min(D_extract, (int) FacFamSeq[Fac].size());

    while (FamsExtracted.size() < D_extract) {
        int Pos = rand() % FacFamSeq[Fac].size();
        int FamExt = FacFamSeq[Fac][Pos];

        FacFamSeq[Fac].erase(FacFamSeq[Fac].begin() + Pos);
        FamsExtracted.push_back(FamExt);
        UpdateFamErasedFromFac(Pos, Fac);
    }


    //从工件组里面提取工件
    for (int i = 0; i < FamsExtracted.size(); i++) //工件组里面的工件提取一半
    {
        int Fam = FamsExtracted[i];
        if (JobSeqInFam[Fam].size() >= 3) {
            for (int j = 0; j <= JobSeqInFam[Fam].size() / 2; j++) {
                int JobPos = rand() % JobSeqInFam[Fam].size();
                int Job = JobSeqInFam[Fam][JobPos];
                JobSeqInFam[Fam].erase(JobSeqInFam[Fam].begin() + JobPos);
                JobsExtracted[Fam].push_back(Job);
            }
        }
    }
}


void
QCCEANoRapid::BasedInd_Destruction_FamsAndJobs_New_Opeator1(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                     vector<int> &FamsExtracted,
                                                     unordered_map<int, vector<int>> &JobsExtracted,
                                                     vector<int> &FacSpan, vector<float> &FacEC) {
    unordered_map<int, pair<int, int>> FamPosErasedFromFac;
    auto UpdateFamErasedFromFac = [&FamPosErasedFromFac](int Pos, int fac) {
        if (FamPosErasedFromFac.find(fac) == end(FamPosErasedFromFac)) {
            FamPosErasedFromFac[fac] = {Pos, Pos - 1};
        } else {
            if (Pos < FamPosErasedFromFac[fac].first) {
                FamPosErasedFromFac[fac].first = Pos;
            }
        }
    };
    FamsExtracted.clear();
    JobsExtracted.clear();
    int CriFac = 0;
    float max = FLT_MIN;
    for (int i = 0; i < FacFamSeq.size(); i++) {
        float p = wyt_rand_include_right(0.0f, 1.0f);
        if (p < 0.5) {
            if (max < FacSpan[i]) {
                max = FacSpan[i];
                CriFac = i;
            }
        } else {
            if (max < FacEC[i]) {
                max = FacEC[i];
                CriFac = i;
            }
        }
    }

    int FamD1 = wyt_rand(2, 7);
    FamD1 = std::min(FamD1, (int) FacFamSeq[CriFac].size());

    while (FamsExtracted.size() < FamD1) {
        int Pos = rand() % FacFamSeq[CriFac].size();
        int FamExt = FacFamSeq[CriFac][Pos];
        FacFamSeq[CriFac].erase(FacFamSeq[CriFac].begin() + Pos);
        FamsExtracted.push_back(FamExt);
        UpdateFamErasedFromFac(Pos, CriFac);
    }
}


void
QCCEANoRapid::BasedInd_Destruction_FamsAndJobs_New_Opeator2(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                     vector<int> &FamsExtracted,
                                                     unordered_map<int, vector<int>> &JobsExtracted,
                                                     vector<int> &FacSpan, vector<float> &FacEC) {
    unordered_map<int, pair<int, int>> FamPosErasedFromFac;
    auto UpdateFamErasedFromFac = [&FamPosErasedFromFac](int Pos, int fac) {
        if (FamPosErasedFromFac.find(fac) == end(FamPosErasedFromFac)) {
            FamPosErasedFromFac[fac] = {Pos, Pos - 1};
        } else {
            if (Pos < FamPosErasedFromFac[fac].first) {
                FamPosErasedFromFac[fac].first = Pos;
            }
        }
    };

    FamsExtracted.clear();
    JobsExtracted.clear();

    int CriFac = 0; // 默认取第一个工厂，避免未初始化警告
    float max = FLT_MIN;
    for (int i = 0; i < FacFamSeq.size(); i++) {
        float p = wyt_rand_include_right(0.0f, 1.0f);
        if (p < 0.5) {
            if (max < FacSpan[i]) {
                max = FacSpan[i];
                CriFac = i;
            }
        } else {
            if (max < FacEC[i]) {
                max = FacEC[i];
                CriFac = i;
            }
        }
    }

    int FamD1 = wyt_rand(2, 7);
    int CriFacFamD = FamD1 / 2;   // 关键工厂抹去一半
    CriFacFamD = std::min(CriFacFamD, (int) FacFamSeq[CriFac].size()); // 避免超出范围

    // 1. **从关键工厂 CriFac 删除族**
    for (int i = 0; i < CriFacFamD; i++) {
        int Pos = rand() % FacFamSeq[CriFac].size();
        int FamExt = FacFamSeq[CriFac][Pos];
        FacFamSeq[CriFac].erase(FacFamSeq[CriFac].begin() + Pos);
        FamsExtracted.push_back(FamExt);
        UpdateFamErasedFromFac(Pos, CriFac);
    }

    // 2. **从其他工厂删除族**
    vector<int> OtherFacs;
    for (int i = 0; i < FacFamSeq.size(); i++) {
        if (i != CriFac && !FacFamSeq[i].empty()) {
            OtherFacs.push_back(i);
        }
    }

    // 确保 OtherFacs 不是空的
    if (!OtherFacs.empty()) {
        int OtherFacFamD = FamD1 - CriFacFamD; // 其他工厂抹去剩下的族
        int deletedCount = 0;

        while (deletedCount < OtherFacFamD && !OtherFacs.empty()) {
            // 重新随机选择一个工厂
            int RandomFacIndex = rand() % OtherFacs.size();
            int SelectedFac = OtherFacs[RandomFacIndex];

            // 计算当前工厂最多能删除的族数
            int extractNum = std::min(OtherFacFamD - deletedCount, (int) FacFamSeq[SelectedFac].size());

            for (int i = 0; i < extractNum; ++i) {
                int Pos = rand() % FacFamSeq[SelectedFac].size();
                int FamExt = FacFamSeq[SelectedFac][Pos];
                FacFamSeq[SelectedFac].erase(FacFamSeq[SelectedFac].begin() + Pos);
                FamsExtracted.push_back(FamExt);
                UpdateFamErasedFromFac(Pos, SelectedFac);
            }

            deletedCount += extractNum;
            // 如果某个工厂已经没有族了，就从 OtherFacs 移除
            if (FacFamSeq[SelectedFac].empty()) {
                OtherFacs.erase(OtherFacs.begin() + RandomFacIndex);
            }
        }
    }
    // 从工件组里面提取工件
    for (int i = 0; i < FamsExtracted.size(); i++) // 工件组里面的工件提取一半
    {
        int Fam = FamsExtracted[i];
        if (JobSeqInFam[Fam].size() >= 3) {
            //            cout << "从族 " << Fam << " 删除的工件: ";
            for (int j = 0; j <= JobSeqInFam[Fam].size() / 2; j++) {
                int JobPos = rand() % JobSeqInFam[Fam].size();
                int Job = JobSeqInFam[Fam][JobPos];
                //                cout << Job << " ";  // 输出被删除的工件
                JobSeqInFam[Fam].erase(JobSeqInFam[Fam].begin() + JobPos);
                JobsExtracted[Fam].push_back(Job);
            }
        }
    }
}


void
QCCEANoRapid::BasedInd_Construction_FamsAndJobs_New_Opeator0(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                      vector<int> &FamsExtracted,
                                                      unordered_map<int, vector<int>> &JobsExtracted, int nadirpointMS,
                                                      float nadirpointTEC, int idealpointMS, float idealpointTEC,vector<Individual> &CCEAPopulation) {

    while (FamsExtracted.size() > 0) {
        int BestFac = -1;
        int BestPos = -1;
        int Pos = rand() % FamsExtracted.size();
        int CurFam = FamsExtracted[Pos];
        FamsExtracted.erase(FamsExtracted.begin() + Pos);
        float Ind = this->FindBestPosToInsertFamForAllFac_Ind_DR_NoAC(FacFamSeq, JobSeqInFam, CurFam,
                                                                     BestFac, BestPos, nadirpointMS, nadirpointTEC,
                                                                     idealpointMS, idealpointTEC,CCEAPopulation);
        FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);

        int FamPos = find(begin(FacFamSeq[BestFac]), end(FacFamSeq[BestFac]), CurFam) - begin(FacFamSeq[BestFac]);

        for (int i = 0; i < JobsExtracted[CurFam].size(); i++) {
            int BestJobPos = -1;
            float minFacInd = INT_MAX;
            int CurJob = JobsExtracted[CurFam][i];
            for (int Pos = 0; Pos <= JobSeqInFam[CurFam].size(); Pos++) {
                float FacInd = GetIndForPerFacAfterInsertJob_DR_NoAC(BestFac, FacFamSeq, FacFamSeq[BestFac], JobSeqInFam,
                                                                   CurFam, FamPos, CurJob,
                                                                    Pos, nadirpointMS, nadirpointTEC, idealpointMS,
                                                                    idealpointTEC, CCEAPopulation);

                if (FacInd < minFacInd) {
                    minFacInd = FacInd;
                    BestPos = Pos;
                }

            }
            JobSeqInFam[CurFam].insert(JobSeqInFam[CurFam].begin() + BestPos, JobsExtracted[CurFam][i]);

        }
    }
}


void
QCCEANoRapid::BasedInd_Construction_FamsAndJobs_New_Opeator1(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                      vector<int> &FamsExtracted,
                                                      unordered_map<int, vector<int>> &JobsExtracted, int nadirpointMS,
                                                      float nadirpointTEC, int idealpointMS, float idealpointTEC,vector<Individual> &CCEAPopulation) {
    // 遍历所有被移除的工件族
    for (auto &CurFam: FamsExtracted) {
        int BestFac = -1;
        int BestPos = -1;
        float minFacInd = INT_MAX;
        // 在所有工厂中找到最优的插入位置
        for (int Fac = 0; Fac < this->m_Factories; ++Fac) {
            int maxPos = FacFamSeq[Fac].size(); // 最大的可能插入位置
            // 遍历所有可能的位置
            for (int Pos = 0; Pos <= maxPos; ++Pos) {
                // 计算当前插入位置的指标
                float tempFacInd = this->GetIndForPerFacAfterInsertFam_DR_NoAC(Fac, FacFamSeq, FacFamSeq[Fac],
                                                                              JobSeqInFam,  CurFam,
                                                                              Pos, nadirpointMS, nadirpointTEC,
                                                                              idealpointMS, idealpointTEC, CCEAPopulation);
                // 更新最优位置
                if (tempFacInd < minFacInd) {
                    minFacInd = tempFacInd;
                    BestPos = Pos;
                    BestFac = Fac;
                    // 尝试相邻位置（局部搜索）
                    for (int offset: {-1, 1}) {
                        int neighborPos = Pos + offset;
                        if (neighborPos >= 0 && neighborPos <= maxPos) {
                            float neighborFacInd = this->GetIndForPerFacAfterInsertFam_DR_NoAC(Fac, FacFamSeq,
                                                                                              FacFamSeq[Fac],
                                                                                              JobSeqInFam,CurFam,
                                                                                              neighborPos, nadirpointMS,
                                                                                              nadirpointTEC,
                                                                                              idealpointMS,
                                                                                              idealpointTEC, CCEAPopulation);
                            // 更新最优位置
                            if (neighborFacInd < minFacInd) {
                                minFacInd = neighborFacInd;
                                BestPos = neighborPos;
                                BestFac = Fac;
                            }
                        }
                    }
                } else {
                    Pos += 2;
                }
            }
        }

        // 插入工件族到最优位置
        FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
    }

}


void
QCCEANoRapid::BasedInd_Construction_FamsAndJobs_New_Opeator2(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                      vector<int> &FamsExtracted,
                                                      unordered_map<int, vector<int>> &JobsExtracted, int nadirpointMS,
                                                      float nadirpointTEC, int idealpointMS, float idealpointTEC,vector<Individual> &CCEAPopulation) {

    // 将挑选出来的d个组插入到最好的位置 Insert the extracted Families into the best Positions
    for (int Fam = 0; Fam < FamsExtracted.size(); Fam++) {

        int CurFam = FamsExtracted[Fam];  //当前组
        int BestFac = -1;
        int BestPos = -1;  //最好的工厂，位置

        if (Fam == FamsExtracted.size() - 1) {
            float Ind = this->FindBestPosToInsertFamForAllFac_Ind_DR_NoAC(FacFamSeq, JobSeqInFam,
                                                                         CurFam, BestFac, BestPos, nadirpointMS,
                                                                         nadirpointTEC, idealpointMS, idealpointTEC, CCEAPopulation);
        } else {
            float randValue = wyt_rand_include_right(0.0f, 1.0f);

            if (randValue < 0.5) {
                FindBestPosToInsertFam_NoAC_Span(FacFamSeq, JobSeqInFam,  CurFam, BestFac, BestPos);
            } else {
                FindBestPosToInsertFamForAllFac_TEC(FacFamSeq, JobSeqInFam, CurFam, BestFac,BestPos);
            }
        }
        // 在最好的工厂的最好位置中插入当前组 Insert CurFam to bestPos at bestFac
        FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);

        int FamPos = find(begin(FacFamSeq[BestFac]), end(FacFamSeq[BestFac]), CurFam) - begin(FacFamSeq[BestFac]);
        for (int i = 0; i < JobsExtracted[CurFam].size(); i++) {
            int BestJobPos = -1;
            float minFacInd = INT_MAX;
            int CurJob = JobsExtracted[CurFam][i];
            for (int Pos = 0; Pos <= JobSeqInFam[CurFam].size(); Pos++) {

                float FacInd = GetIndForPerFacAfterInsertJob_DR_NoAC(BestFac, FacFamSeq, FacFamSeq[BestFac], JobSeqInFam,
                                                                  CurFam, FamPos, CurJob,
                                                                    Pos, nadirpointMS, nadirpointTEC, idealpointMS,
                                                                    idealpointTEC, CCEAPopulation);

                if (FacInd < minFacInd) {
                    minFacInd = FacInd;
                    BestPos = Pos;
                }

            }
            JobSeqInFam[CurFam].insert(JobSeqInFam[CurFam].begin() + BestPos, JobsExtracted[CurFam][i]);
        }
    }
}


void QCCEANoRapid::UpdateArchiveGroupJobSet(int Muu)
{

    vector<Individual> Temp;
    Temp.clear();

    //占优关系
    Pareto_relation(m_QCCEANoRapidPopulation);

    //弱化支配解 de-emphasize dominated solutions
    for (int j = 0; j < m_QCCEANoRapidPopulation.size(); j++)
    {
        for (int i = 0; i < m_QCCEANoRapidPopulation.size(); i++)
        {
            if (m_QCCEANoRapidPopulation[i].flag == 0)
            {
                if (m_QCCEANoRapidPopulation[i].pareto_rel[j] == 1)
                {
                    m_QCCEANoRapidPopulation[i].flag = 999;
                }
            }
        }
    }


    for (int i = 0; i < m_QCCEANoRapidPopulation.size(); i++)
    {
        if (m_QCCEANoRapidPopulation[i].flag == 0)
            Temp.push_back(m_QCCEANoRapidPopulation[i]);
    }

    //参数：mu(Archive种群倍数)
    int mu = Muu;
    if (Temp.size() < mu * m_RefSize)
    {
        while (true)
        {
            for (int i = 0; i < m_QCCEANoRapidPopulation.size(); i++)
            {
                if (m_QCCEANoRapidPopulation[i].flag == 999)
                {
                    Temp.push_back(m_QCCEANoRapidPopulation[i]);
                    if (Temp.size() == 10 * m_RefSize)
                        break;
                }
            }
            break;
        }
    }

    m_QCCEANoRapidPopulation.clear();
    for (int i = 0; i < Temp.size(); i++)
        m_QCCEANoRapidPopulation.push_back(Temp[i]);


    //归一化
    Normalize(m_QCCEANoRapidPopulation, m_NadirPointMS, m_NadirPointTEC, m_IdealPointMS, m_IdealPointTEC);

    //计算收敛指标
    calc_convergence_ind(m_QCCEANoRapidPopulation);

    //计算分散指标
    calc_distribution_ind(m_QCCEANoRapidPopulation);

    int n, n1, n2, n3, nrank;
    int the_one;
    n = 0;   //Qs.size  预选择的个体集合
    n1 = 0;  //Q.size   参考集的个体数量
    n2 = 0;  //Qth.size 由distribution threshold弱化的集合
    n3 = 0;  //Qd.size  被支配个体的集合

    vector<Individual> TempQCCEANoRapidPopulation;
    TempQCCEANoRapidPopulation.clear();

    //占优关系
    Pareto_relation(m_QCCEANoRapidPopulation);

    n1 = m_QCCEANoRapidPopulation.size();

    while (n1 > 0)
    {
        for (int i = 0; i < m_QCCEANoRapidPopulation.size(); i++)
        {
            if (m_QCCEANoRapidPopulation[i].flag == 0)
            {
                the_one = i;
                break;
            }
        }

        m_QCCEANoRapidPopulation[the_one].flag = 1;
        TempQCCEANoRapidPopulation.push_back(m_QCCEANoRapidPopulation[the_one]);
        n++;


        //弱化附近解 de-emphasize neighbors
        for (int i = 0; i < m_QCCEANoRapidPopulation.size(); i++)
        {
            if (m_QCCEANoRapidPopulation[i].flag == 0)
            {
                if (m_QCCEANoRapidPopulation[i].distribution_ind[the_one] < thr_zeta)
                {
                    m_QCCEANoRapidPopulation[i].flag = 999;
                    n2++;
                }
            }
        }

        //弱化支配解 de-emphasize dominated solutions
        for (int i = 0; i < m_QCCEANoRapidPopulation.size(); i++)
        {
            if (m_QCCEANoRapidPopulation[i].flag == 0)
            {
                if (m_QCCEANoRapidPopulation[i].pareto_rel[the_one] == 1)
                {
                    m_QCCEANoRapidPopulation[i].flag = 999;
                    n3++;
                }
            }
        }

        //剩余数量 number of the rest
        n1 = 0;
        for (int i = 0; i < m_QCCEANoRapidPopulation.size(); i++)
        {
            if (m_QCCEANoRapidPopulation[i].flag == 0)
            {
                n1++;
            }
        }
    }

    //参数 threshold
    float oldthr_zeta = thr_zeta;
    float ratio = n * 1.0f / (m_RefSize * 1.0f);
    if (n3 < (mu - 1) * m_RefSize)
        thr_zeta = thr_zeta * exp((ratio - 1.0f) / (2 * 1.0f));
    else
        thr_zeta = oldthr_zeta;


    while (n < m_RefSize)
    {
        //剩余数量 number of the rest
        n1 = 0;
        for (int i = 0; i < m_QCCEANoRapidPopulation.size(); i++)
        {
            if (m_QCCEANoRapidPopulation[i].flag == 0)
            {
                n1++;
            }
        }

        if (n1 == 0)
        {
            for (int i = 0; i < m_QCCEANoRapidPopulation.size(); i++)
            {
                if (m_QCCEANoRapidPopulation[i].flag == 999)
                {
                    m_QCCEANoRapidPopulation[i].flag = 0;
                }
            }
        }

        //选择一个
        for (int i = 0; i < m_QCCEANoRapidPopulation.size(); i++)
        {
            if (m_QCCEANoRapidPopulation[i].flag == 0)
            {
                the_one = i;
                break;
            }
        }
        m_QCCEANoRapidPopulation[the_one].flag = 1;
        TempQCCEANoRapidPopulation.push_back(m_QCCEANoRapidPopulation[the_one]);
        n++;

        //弱化附近解 de-emphasize neighbors
        for (int i = 0; i < m_QCCEANoRapidPopulation.size(); i++)
        {
            if (m_QCCEANoRapidPopulation[i].flag == 0)
            {
                if (m_QCCEANoRapidPopulation[i].distribution_ind[the_one] < oldthr_zeta)
                {
                    m_QCCEANoRapidPopulation[i].flag = 999;
                    n2++;
                }
            }
        }

        //弱化支配解 de-emphasize dominated solutions
        for (int i = 0; i < m_QCCEANoRapidPopulation.size(); i++)
        {
            if (m_QCCEANoRapidPopulation[i].flag == 0)
            {
                if (m_QCCEANoRapidPopulation[i].pareto_rel[the_one] == 1)
                {
                    m_QCCEANoRapidPopulation[i].flag = 999;
                    n3++;
                }
            }
        }
    }

    vector<Individual> Temp1;
    Temp1.clear();
    // m_RefSize
    for (int i = 0; i < TempQCCEANoRapidPopulation.size(); i++)
    {
        Temp1.push_back(TempQCCEANoRapidPopulation[i]);
    }

    //变速
    Speed_mutation_NoAC(TempQCCEANoRapidPopulation, Temp1, m_QCCEANoRapidPopulation);

    for (int PS = 0; PS < m_RefSize; PS++)
    {
        // 更新档案集（赋值）
        m_ArchiveJobSeqinFamArray[PS] = TempQCCEANoRapidPopulation[PS].m_JobSeqInFamArray;  //工件组序列
        m_ArchiveFacFamSeqArray[PS] = TempQCCEANoRapidPopulation[PS].m_FacFamSeqArray;   //工厂组序列
        m_ArchiveSpanArray[PS] = TempQCCEANoRapidPopulation[PS].MS;  // 最大完工时间
        m_ArchiveTECArray[PS] = TempQCCEANoRapidPopulation[PS].TEC;  //总能耗
        m_ArchiveSpeedVector[PS] = TempQCCEANoRapidPopulation[PS].m_SpeedVector;
        //更新组种群
        m_FacFamSeqArray[PS] = m_ArchiveFacFamSeqArray[PS];
        //更新工件种群
        m_JobSeqinFamArray[PS] = m_ArchiveJobSeqinFamArray[PS];
        m_SpeedMatrix = m_ArchiveSpeedVector[PS];
    }
}

void QCCEANoRapid::SaveTECandDeMS(vector<Individual> &CCEAPopulation) {
    vector<Individual> tempQCCEANoRapidPopulation;
    tempQCCEANoRapidPopulation.clear();
    vector<Individual> tempQCCEANoRapidPopulation2;
    tempQCCEANoRapidPopulation2.clear();
    vector<Individual> tempQCCEANoRapidPopulation3;
    tempQCCEANoRapidPopulation3.clear();
    vector<Individual> tempQCCEANoRapidPopulation4;
    tempQCCEANoRapidPopulation4.clear();

    int orgSize = CCEAPopulation.size();
    for (int PS = 0; PS < orgSize; PS++) {
        tempQCCEANoRapidPopulation.push_back(CCEAPopulation[PS]);
        bool flag2 = false;

        vector<int> FacSpan(m_Factories);//组，工厂完工时间
        vector<vector<int>> JFDTime(this->m_Jobs, vector<int>(this->m_Machines, 0));
        //节能策略1
        vector<vector<int>> DelayTime;
        DelayTime.clear();
        DelayTime.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
        GetDelayTime_Forward_NoAC(tempQCCEANoRapidPopulation[PS].m_FacFamSeqArray, tempQCCEANoRapidPopulation[PS].m_JobSeqInFamArray,
                            FacSpan, DelayTime);

        bool sign = false;

        for (int j = 0; j < m_Jobs; j++) {
            for (int i = 1; i < m_Machines; i++) {
                if ((DelayTime[j][i] > 0) && (i < m_Machines - 1)) {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] /m_Speed[tempQCCEANoRapidPopulation[PS].m_SpeedVector[j][i]]);

                    m_TureJobOpertime[j][i + 1] = static_cast<int>(m_JobOperPTime[j][i + 1] /m_Speed[tempQCCEANoRapidPopulation[PS].m_SpeedVector[j][i + 1]]);

                    if ((JFDTime[j][i] < (JFDTime[j][i + 1] - m_TureJobOpertime[j][i + 1]))) {
                        int Speedlevel = tempQCCEANoRapidPopulation[PS].m_SpeedVector[j][i] - 1;

                        for (int level = Speedlevel; level > 0; level--) {
                            if ((static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[level])) -
                                 m_TureJobOpertime[j][i]) < DelayTime[j][i]) {
                                if ((m_TempJFDTime[j][i] +
                                     (static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[level])) -
                                      m_TureJobOpertime[j][i])) <
                                    (JFDTime[j][i + 1] - m_TureJobOpertime[j][i + 1])) {
                                    tempQCCEANoRapidPopulation[PS].m_SpeedVector[j][i] = level;
                                    sign = true;
                                } else
                                    break;
                            } else
                                break;
                        }
                    }

                } else if ((DelayTime[j][i] > 0) && (i == m_Machines - 1)) {
                    int Speedlevel = tempQCCEANoRapidPopulation[PS].m_SpeedVector[j][i] - 1;
                    m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] /
                                                                      m_Speed[tempQCCEANoRapidPopulation[PS].m_SpeedVector[j][i]]));
                    for (int level = Speedlevel; level > 0; level--) {
                        if ((static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[level])) -
                             m_TureJobOpertime[j][i]) < DelayTime[j][i]) {
                            tempQCCEANoRapidPopulation[PS].m_SpeedVector[j][i] = level;
                            sign = true;
                        } else
                            break;
                    }
                }
            }//end machine
        }//end job

        if (sign == true) {
            int Makespan1 = 0;
            float TotalEC1 = 0.0;

            m_SpeedMatrix = tempQCCEANoRapidPopulation[PS].m_SpeedVector;

            //得到真正的处理时间 及单位加工时间能耗
            for (int j = 0; j < m_Jobs; j++) {
                for (int i = 0; i < m_Machines; i++) {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }


            Makespan1 = GetSpan(tempQCCEANoRapidPopulation[PS].m_FacFamSeqArray,tempQCCEANoRapidPopulation[PS].m_JobSeqInFamArray,FacSpan);

            TotalEC1 = GetTECForAllFac(tempQCCEANoRapidPopulation[PS].m_FacFamSeqArray,tempQCCEANoRapidPopulation[PS].m_JobSeqInFamArray );

            //检查
            //CheckSol(tempQCCEANoRapidPopulation[PS].m_FacFamSeqArray, tempQCCEANoRapidPopulation[PS].m_JobSeqInFamArray, Makespan1);
            //CheckSolTEC(tempQCCEANoRapidPopulation[PS].m_FacFamSeqArray, tempQCCEANoRapidPopulation[PS].m_JobSeqInFamArray, TotalEC1);
            tempQCCEANoRapidPopulation[PS].MS = Makespan1;
            tempQCCEANoRapidPopulation[PS].TEC = TotalEC1;

            //判断
            if ((CCEAPopulation[PS].MS >= tempQCCEANoRapidPopulation[PS].MS) &&
                (CCEAPopulation[PS].TEC > tempQCCEANoRapidPopulation[PS].TEC)) {
                CCEAPopulation.push_back(tempQCCEANoRapidPopulation[PS]);
                flag2 = true;
            }
        }


        //节能策略3：在1的基础上，在不超过Makespan的情况下，降低非关键工厂的v，来降低TEC
        vector<int> FacSpan3(m_Factories);
        if (flag2) {
            //float gapTEC = CCEAPopulation[PS].TEC - tempQCCEANoRapidPopulation[PS].TEC;
            tempQCCEANoRapidPopulation3.push_back(tempQCCEANoRapidPopulation[PS]);
            vector<int> notcirfac;
            notcirfac.clear();
            for (int fac = 0; fac < m_Factories; fac++) {
                if (tempQCCEANoRapidPopulation3.back().MS > FacSpan[fac]) {
                    notcirfac.push_back(fac);
                }
                FacSpan3[fac] = FacSpan[fac];
            }
            int diedai = 0;
            float tempTEC = 0.0;

            do {
                diedai++;
                if (diedai == 1)
                    tempTEC = tempQCCEANoRapidPopulation[PS].TEC;

                if (diedai == 2)
                    tempTEC = tempQCCEANoRapidPopulation3.back().TEC;
                if (diedai == 3)
                    break;

                bool judgeV = true;
                vector<vector<int>> tempSpeedMatrix;
                tempSpeedMatrix.clear();
                tempSpeedMatrix.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
                tempSpeedMatrix = tempQCCEANoRapidPopulation3.back().m_SpeedVector;

                //对非关键工厂进行变速
                for (int f = 0; f < notcirfac.size(); f++) {
                    judgeV = true;

                    for (int fam = 0; fam < tempQCCEANoRapidPopulation3.back().m_FacFamSeqArray[notcirfac[f]].size(); fam++) {
                        int CurFam = tempQCCEANoRapidPopulation3.back().m_FacFamSeqArray[notcirfac[f]][fam];
                        for (int job = 0; job < tempQCCEANoRapidPopulation3.back().m_JobSeqInFamArray[CurFam].size(); job++) {
                            int CurJob = tempQCCEANoRapidPopulation3.back().m_JobSeqInFamArray[CurFam][job];
                            for (int m = 0; m < m_Machines; m++) {
                                if (tempQCCEANoRapidPopulation3.back().m_SpeedVector[CurJob][m] > 0) {
                                    tempQCCEANoRapidPopulation3.back().m_SpeedVector[CurJob][m] =
                                            tempQCCEANoRapidPopulation3.back().m_SpeedVector[CurJob][m] - 1;
                                    judgeV = false;

                                    //真实的处理时间
                                    //m_TureJobOpertime[CurJob][m] = static_cast<int>(100 * (m_JobOperPTime[CurJob][m] / m_Speed[speedLevel]));
                                }

                            }
                        }
                    }
                    //若有速度有变化
                    if (judgeV == false) {
                        //faccount++;
                        m_SpeedMatrix = tempQCCEANoRapidPopulation3.back().m_SpeedVector;

                        //得到真正的处理时间 及 单位加工能耗PEC系数
                        for (int j = 0; j < m_Jobs; j++) {
                            for (int i = 0; i < m_Machines; i++) {
                                m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] /
                                                                                  m_Speed[m_SpeedMatrix[j][i]]));
                                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                            }
                        }

                        //得到工厂f变速后的Span
                        int tempnotcirSpan = GetSpan_Forward_InFactory(tempQCCEANoRapidPopulation3.back().m_FacFamSeqArray[notcirfac[f]],tempQCCEANoRapidPopulation3.back().m_JobSeqInFamArray);//得到每个工厂的完工时间
                        //如果改变速度后未超过Makspan，才真正变速
                        if (tempnotcirSpan > tempQCCEANoRapidPopulation3.back().MS) {
                            tempQCCEANoRapidPopulation3.back().m_SpeedVector = tempSpeedMatrix;
                            //countfac++;
                            //judgeV = true;
                        } else {
                            tempSpeedMatrix = tempQCCEANoRapidPopulation3.back().m_SpeedVector;
                        }

                    }

                }

                /*if (judgeV && (countfac == m_Factories))
                    break;*/

                int Makespan3 = 0;
                float TotalEC3 = 0.0;

                m_SpeedMatrix = tempQCCEANoRapidPopulation3.back().m_SpeedVector;

                //得到真正的处理时间 及 单位加工能耗PEC系数
                for (int j = 0; j < m_Jobs; j++) {
                    for (int i = 0; i < m_Machines; i++) {
                        m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] /m_Speed[m_SpeedMatrix[j][i]]);
                        UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                    }
                }

                Makespan3 = GetSpan(tempQCCEANoRapidPopulation3.back().m_FacFamSeqArray,tempQCCEANoRapidPopulation3.back().m_JobSeqInFamArray,FacSpan3);

                TotalEC3 = GetTECForAllFac(tempQCCEANoRapidPopulation3.back().m_FacFamSeqArray,tempQCCEANoRapidPopulation3.back().m_JobSeqInFamArray);

                //检查
                //CheckSol(tempQCCEANoRapidPopulation3.back().m_FacFamSeqArray, tempQCCEANoRapidPopulation3.back().m_JobSeqInFamArray, Makespan3);
                //CheckSolTEC(tempQCCEANoRapidPopulation3.back().m_FacFamSeqArray, tempQCCEANoRapidPopulation3.back().m_JobSeqInFamArray, TotalEC3);
                tempQCCEANoRapidPopulation3.back().MS = Makespan3;
                tempQCCEANoRapidPopulation3.back().TEC = TotalEC3;

                //判断
                if ((CCEAPopulation[PS].MS >= tempQCCEANoRapidPopulation3.back().MS) &&
                    (tempQCCEANoRapidPopulation[PS].TEC > tempQCCEANoRapidPopulation3.back().TEC)) {
                    CCEAPopulation.push_back(tempQCCEANoRapidPopulation3.back());
                }

            } while (tempTEC > tempQCCEANoRapidPopulation3.back().TEC);

        }

        //节能策略4：在3的基础上，在不超过原TEC的情况下，提高关键工厂的v，来降低MS
        vector<int> FacSpan4(m_Factories);
        if (flag2) {
            //float gapTEC = CCEAPopulation[PS].TEC - tempQCCEANoRapidPopulation[PS].TEC;
            tempQCCEANoRapidPopulation4.push_back(tempQCCEANoRapidPopulation3.back());
            vector<int> cirfac;
            cirfac.clear();
            for (int fac = 0; fac < m_Factories; fac++) {
                if (tempQCCEANoRapidPopulation4.back().MS == FacSpan[fac]) {
                    cirfac.push_back(fac);
                }
                FacSpan4[fac] = FacSpan[fac];
            }
            int diedai = 0;
            int tempMS = 0;

            do {
                diedai++;
                if (diedai == 1)
                    tempMS = CCEAPopulation[PS].MS;

                if (diedai == 2) {
                    tempMS = tempQCCEANoRapidPopulation4.back().MS;
                    //cirfac.clear();
                    vector<int> tempcirfac;
                    tempcirfac.clear();
                    //判断关键工厂是否改变
                    for (int fac = 0; fac < m_Factories; fac++) {
                        if (tempQCCEANoRapidPopulation4.back().MS == FacSpan4[fac]) {
                            tempcirfac.push_back(fac);
                        }
                    }
                    bool judgecirfacifchange = true;
                    for (int i = 0; i < tempcirfac.size(); i++) {
                        if (tempcirfac[i] != cirfac[i]) {
                            judgecirfacifchange = false;
                            break;
                        }
                    }
                    if (!judgecirfacifchange) {
                        diedai = 1;
                        cirfac.clear();
                        cirfac = tempcirfac;
                    }

                }

                if (diedai == 3)
                    break;

                bool judgeV = true;
                vector<vector<int>> tempSpeedMatrix;
                tempSpeedMatrix.clear();
                tempSpeedMatrix.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
                tempSpeedMatrix = tempQCCEANoRapidPopulation4.back().m_SpeedVector;

                //int countfac = 0;
                //int faccount = 0;
                //对非关键工厂进行变速
                for (int f = 0; f < cirfac.size(); f++) {
                    judgeV = true;
                    //int speedLevel;
                    for (int fam = 0; fam < tempQCCEANoRapidPopulation4.back().m_FacFamSeqArray[cirfac[f]].size(); fam++) {
                        int CurFam = tempQCCEANoRapidPopulation4.back().m_FacFamSeqArray[cirfac[f]][fam];
                        for (int job = 0; job < tempQCCEANoRapidPopulation4.back().m_JobSeqInFamArray[CurFam].size(); job++) {
                            int CurJob = tempQCCEANoRapidPopulation4.back().m_JobSeqInFamArray[CurFam][job];
                            for (int m = 0; m < m_Machines; m++) {
                                if (tempQCCEANoRapidPopulation4.back().m_SpeedVector[CurJob][m] < 2) {
                                    tempQCCEANoRapidPopulation4.back().m_SpeedVector[CurJob][m] =
                                            tempQCCEANoRapidPopulation4.back().m_SpeedVector[CurJob][m] + 1;
                                    judgeV = false;
                                }

                            }
                        }
                    }
                    //若有速度有变化
                    if (judgeV == false) {
                        //faccount++;
                        m_SpeedMatrix = tempQCCEANoRapidPopulation4.back().m_SpeedVector;

                        //得到真正的处理时间 及 单位加工能耗PEC系数
                        for (int j = 0; j < m_Jobs; j++) {
                            for (int i = 0; i < m_Machines; i++) {
                                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] /m_Speed[m_SpeedMatrix[j][i]]);
                                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                            }
                        }

                        //得到工厂f变速后的TEC
                      //  GetJFDTime_Forward(tempQCCEANoRapidPopulation4.back().m_FacFamSeqArray,tempQCCEANoRapidPopulation4.back().m_JobSeqInFamArray, m_TempJFDTime, FacSpan4);
                        float tempTotalEC = GetTECForAllFac(tempQCCEANoRapidPopulation4.back().m_FacFamSeqArray,
                                                                 tempQCCEANoRapidPopulation4.back().m_JobSeqInFamArray);


                        //如果改变速度后未超过原TEC，才真正变速
                        if (tempTotalEC > CCEAPopulation[PS].TEC) {
                            tempQCCEANoRapidPopulation4.back().m_SpeedVector = tempSpeedMatrix;
                            break;
                        }

                    }

                }

                int Makespan4 = 0;
                float TotalEC4 = 0.0;

                m_SpeedMatrix = tempQCCEANoRapidPopulation4.back().m_SpeedVector;

                //得到真正的处理时间 及 单位加工能耗PEC系数
                for (int j = 0; j < m_Jobs; j++) {
                    for (int i = 0; i < m_Machines; i++) {
                        m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] /
                                                                          m_Speed[m_SpeedMatrix[j][i]]));
                        UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                    }
                }

                Makespan4 = GetSpan(tempQCCEANoRapidPopulation4.back().m_FacFamSeqArray,tempQCCEANoRapidPopulation4.back().m_JobSeqInFamArray,  FacSpan4);

                TotalEC4 = GetTECForAllFac(tempQCCEANoRapidPopulation4.back().m_FacFamSeqArray,tempQCCEANoRapidPopulation4.back().m_JobSeqInFamArray);

                //检查
                //CheckSol(tempQCCEANoRapidPopulation4.back().m_FacFamSeqArray, tempQCCEANoRapidPopulation4.back().m_JobSeqInFamArray, Makespan4);
                //CheckSolTEC(tempQCCEANoRapidPopulation4.back().m_FacFamSeqArray, tempQCCEANoRapidPopulation4.back().m_JobSeqInFamArray, TotalEC4);
                tempQCCEANoRapidPopulation4.back().MS = Makespan4;
                tempQCCEANoRapidPopulation4.back().TEC = TotalEC4;

                //判断
                if ((CCEAPopulation[PS].MS > tempQCCEANoRapidPopulation4.back().MS) &&
                    (CCEAPopulation[PS].TEC >= tempQCCEANoRapidPopulation4.back().TEC)) {
                    CCEAPopulation.push_back(tempQCCEANoRapidPopulation4.back());
                }

            } while (tempMS > tempQCCEANoRapidPopulation4.back().MS);

        }
    }

}

//Q learning

int QCCEANoRapid::DetermineAction(int state, const vector<vector<double>> &Q) {

    double rande = wyt_rand_include_right(0.0f, 1.0f);
    if (rande < e) {
        vector<int> QI = findMaxIndices(Q[state]);
        return int(QI[wyt_rand(QI.size())]);
    } else {
        return wyt_rand(Q[state].size());
    }
}

//int QCCEANoRapid::DetermineState(int PreMS, int actMS, float PreTEC, float actTEC) {
//    if (actMS < PreMS && actTEC < PreTEC) {
//        return 0;//好
//    } else if (actTEC < PreTEC || actMS < PreMS) {
//        return 1;//改进一个解 较好
//    } else {
//        return 2;//没改进
//    }
//}


// 比较函数，用于按makespan排序
bool QCCEANoRapid::CompareByMakespan(const Individual &a, const Individual &b) {
    return a.MS < b.MS;
}

// 比较函数，用于按TEC排序
bool QCCEANoRapid::CompareByTEC(const Individual &a, const Individual &b) {
    return a.TEC < b.TEC;
}

// 确定个体的状态
int QCCEANoRapid::DetermineState(int ObjectMS, float ObjectTEC, const std::vector<Individual> &QCCEANoRapidPopulation) {

    int N = QCCEANoRapidPopulation.size();

    // 按makespan排序
    std::vector<Individual> SortedByMakespan = QCCEANoRapidPopulation;
    std::sort(SortedByMakespan.begin(), SortedByMakespan.end(), CompareByMakespan);
    int MedianMakeSpan = SortedByMakespan[N / 2].MS;

    // 按TEC排序
    std::vector<Individual> SortedByTEC = QCCEANoRapidPopulation;
    std::sort(SortedByTEC.begin(), SortedByTEC.end(), CompareByTEC);
    float MedianTEC = SortedByTEC[N / 2].TEC;

    // 确定状态
    if (ObjectMS <= MedianMakeSpan && ObjectTEC <= MedianTEC) {
        return 0; // 状态0：MS和TEC均较好
    } else if (ObjectMS <= MedianMakeSpan && ObjectTEC > MedianTEC) {
        return 1; // 状态1：MS较好，TEC较差
    } else if (ObjectMS > MedianMakeSpan && ObjectTEC <= MedianTEC) {
        return 1; // 状态2：MS较差，TEC较好
    } else {
        return 2; // 状态3：MS和TEC均较差
    }
}

void QCCEANoRapid::UpdateQ(int oldstate, int newstate, int act, double reward, vector<vector<double>> &Q) {
    //double e = 0.75;贪婪策略  discountFactor = 0.2;折扣因子 rewardFactor = 0.75;学习率
    if (newstate == 0 || newstate == 1) {
        Q[oldstate][act] = (1.0 - rewardFactor) * Q[oldstate][act] + rewardFactor * (reward + discountFactor *
                                                                                              *max_element(
                                                                                                      Q[newstate].begin(),
                                                                                                      Q[newstate].end()));
    }
}

double QCCEANoRapid::CaculateReward(int oldMS, int newMS, float oldTEC, float newTEC,
                             const std::vector<Individual> &population) {

    // 计算目标函数的最大值和最小值
    int MS_max = std::numeric_limits<int>::lowest();
    int MS_min = std::numeric_limits<int>::max();
    float TEC_max = std::numeric_limits<float>::lowest();
    float TEC_min = std::numeric_limits<float>::max();

    for (const auto &ind: population) {
        MS_max = std::max<int>(MS_max, ind.MS);
        MS_min = std::min<int>(MS_min, ind.MS);
        TEC_max = std::max<float>(TEC_max, ind.TEC);
        TEC_min = std::min<float>(TEC_min, ind.TEC);
    }

    // 防止除零错误
    if (MS_max == MS_min || TEC_max == TEC_min) {
        return 0.0;  // 或者根据情况返回一个合理的默认奖励值
    }

    // 计算奖励值
    double reward = 0.0;
    if (newMS <= oldMS && newTEC <= oldTEC) {   // 计算新解和旧解在每个目标上的支配关系
        reward = 2.0 + ((MS_max - newMS) / (MS_max - MS_min)) + ((TEC_max - newTEC) / (TEC_max - TEC_min));
    } else if
            (oldMS <= newMS && oldTEC <= newTEC) {
        reward = 0.0;
    } else {
        reward = 1.0 + ((MS_max - newMS) / (MS_max - MS_min)) + ((TEC_max - newTEC) / (TEC_max - TEC_min));
    }

    return reward;
}


void QCCEANoRapid::PerformAction_Fams(int act, vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqinFam,

                               vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC,
                               int &ObjectMS, float &ObjectTEC) {
    switch (act) {
        case 0:
            cout<<"0"<<endl;
            //对makespan最大的工厂中组在所有工厂中重新找最好的位置（基于收敛指标插入）
            BasedInd_RandFamInFacTobestPos(FacFamSeq, JobSeqinFam, m_NadirPointMS,
                                               m_NadirPointTEC, m_IdealPointMS, m_IdealPointTEC,
                                               CCEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
            cout<<"0 end "<<endl;
            break;
        case 1:
            cout<<"1"<<endl;
            BasedInd_SwapFam(FacFamSeq, JobSeqinFam,  m_NadirPointMS, m_NadirPointTEC,
                                 m_IdealPointMS, m_IdealPointTEC, CCEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
            cout<<"1 end "<<endl;
            break;
        case 2:
            cout<<"2"<<endl;
            BasedInd_DeAndConFams(FacFamSeq, JobSeqinFam,  m_NadirPointMS, m_NadirPointTEC,
                                      m_IdealPointMS, m_IdealPointTEC,
                                      CCEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
            cout<<"2 end "<<endl;
            break;
        case 3:
            cout<<"3  "<<endl;
            BasedIndLS_SetupSwap(FacFamSeq, JobSeqinFam, m_NadirPointMS, m_NadirPointTEC,
                                     m_IdealPointMS, m_IdealPointTEC, CCEAPopulation, FacSpan, FacEC, ObjectMS,
                                     ObjectTEC);
            cout<<"3 end "<<endl;
            break;
        case 4:
            cout<<"4  "<<endl;
            BasedInd_ShiftFam(FacFamSeq, JobSeqinFam,m_NadirPointMS, m_NadirPointTEC,
                                  m_IdealPointMS, m_IdealPointTEC,
                                  CCEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
            cout<<"4 end   "<<endl;
            break;

        default:
            std::cerr << "Error: Invalid action " << act << " in PerformAction_Fams." << std::endl;
    }
    // 重新计算总的 MS 和 TEC 确保最新解
    GetMSandTECForPerandToalFac(FacFamSeq, JobSeqinFam, FacSpan, FacEC, ObjectMS, ObjectTEC);
}

vector<int> QCCEANoRapid::findMaxIndices(const vector<double> &vec) {
    vector<int> maxIndices;
    if (vec.empty()) {
        return maxIndices; // 返回空的索引向量
    }
    double maxVal = vec[0];
    maxIndices.push_back(0);

    for (int i = 1; i < vec.size(); ++i) {
        if (vec[i] > maxVal) {
            maxVal = vec[i];
            maxIndices.clear();
            maxIndices.push_back(i);
        } else if (vec[i] == maxVal) {
            maxIndices.push_back(i);
        }
    }

    return maxIndices;
}

void QCCEANoRapid::InitialPop() {

    //初始化档案集 Initialize are set and Set flags
    for (int PS = 0; PS < m_RefSize; PS++) {

        //得到真正的处理时间 及单位时间加工能耗
        for (int j = 0; j < m_Jobs; j++) {
            for (int i = 0; i < m_Machines; i++) {
                //m_SpeedMatrix[j][i] = rand() % 3;
                m_SpeedMatrix[j][i] = 0;
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }

        vector<int> FamPrmu(m_Families), FacSpan(m_Factories);//组，工厂完工时间
        vector<float> FacEC(m_Factories);
        vector<vector<int>> JobSeqInFam, FacFamSeq; //工件在组的序列，组在工厂的序列，新组在工厂的序列
        if (PS == 0 || PS == 1) //第一个参考解用LPT，其它随机生成 Generate job sequence in each family and family sequence using LPT
        {
            this->SortJobsInFam(0, JobSeqInFam); //LPT
            this->CombineSortedSequences(2, 3, FamPrmu);

        } else if ( PS == 2|| PS == 3) {
            this->SortJobsInFam(1, JobSeqInFam); //SPT
            this->CombineSortedSequences(4, 5, FamPrmu);
        } else // 随机生成工件序列和组序列 Generate job sequence in each family and family sequence randomly
        {
            JobSeqInFam = this->m_JobsInEachFamily;
            for (int fam = 0; fam < JobSeqInFam.size(); fam++)
                shuffle(JobSeqInFam[fam].begin(), JobSeqInFam[fam].end(), std::mt19937(std::random_device()()));//打乱组内工件顺序
            for (int fam = 0; fam < FamPrmu.size(); fam++)
                FamPrmu[fam] = fam;
            shuffle(FamPrmu.begin(), FamPrmu.end(), std::mt19937(std::random_device()()));//打乱组顺序
            // 输出生成的组序列
        }

        FacFamSeq.clear();
        FacFamSeq.resize(this->m_Factories);

        int CurFam = -1;
        int BestFac = -1, BestPos = -1;

        int Fam = 0;
        for (int Fac = 0; Fac < FacFamSeq.size(); Fac++) {
            CurFam = FamPrmu[Fam];
            FacFamSeq[Fac].insert(FacFamSeq[Fac].begin() + 0, CurFam);
            Fam++;
        }
        if (PS == 0 || PS == 1) {
            for (; Fam < FamPrmu.size(); Fam++) {
                CurFam = FamPrmu[Fam];
                BestFac = -1;
                BestPos = -1;
                this-> FindBestPosToInsertFam_NoAC_Span(FacFamSeq, JobSeqInFam, CurFam, BestFac,BestPos);

                FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
            }
            // 输出每个工厂的组序列，并标明对应的种群
//            for (int fac = 0; fac < FacFamSeq.size(); ++fac) {
//                std::cout << "Ps " << PS << " Fac " << fac << " SequenceOfGroup: ";
//                for (int fam: FacFamSeq[fac]) {
//                    std::cout << fam << " ";
//                }
//                std::cout << std::endl;
//            }
        } else if (PS == 2 || PS == 3) {
            for (; Fam < FamPrmu.size(); Fam++) {
                CurFam = FamPrmu[Fam];
                BestFac = -1;
                BestPos = -1;
                this->FindBestPosToInsertFamForAllFac_TEC(FacFamSeq, JobSeqInFam, CurFam, BestFac, BestPos);
                FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
            }
            // 输出每个工厂的组序列，并标明对应的种群
//            for (int fac = 0; fac < FacFamSeq.size(); ++fac) {
//                std::cout << "Ps " << PS << " Fac " << fac << " SequenceOfGroup: ";
//                for (int fam: FacFamSeq[fac]) {
//                    std::cout << fam << " ";
//                }
//                std::cout << std::endl;
//            }
        } else {
            for (; Fam < FamPrmu.size(); Fam++) {
                CurFam = FamPrmu[Fam];
                BestFac = -1;
                BestPos = -1;
                this->FindBestPosToInsertFamForAllFacs_MC(FacFamSeq, JobSeqInFam, CurFam, BestFac, BestPos);
                FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
            }

            // 输出每个工厂的组序列，并标明对应的种群
//            for (int fac = 0; fac < FacFamSeq.size(); ++fac) {
//                std::cout << "Ps " << PS << " Fac " << fac << " SequenceOfGroup: ";
//                for (int fam: FacFamSeq[fac]) {
//                    std::cout << fam << " ";
//                }
//                std::cout << std::endl;
//            }
        }

        int MS = 0;
        float TEC = 0;
        // cout << endl << "种群中第" << PS << "个个体" << endl;
        GetMSandTECForPerandToalFac(FacFamSeq, JobSeqInFam, FacSpan, FacEC, MS, TEC);
        // cout << "Makespan：" << MS << "\t" << "TEC：" << TEC << endl;

        // 初始化档案集（赋值）
        m_ArchiveJobSeqinFamArray[PS] = JobSeqInFam;  //工件组序列
        m_ArchiveFacFamSeqArray[PS] = FacFamSeq;   //工厂组序列
        m_ArchiveFacSpanArray[PS] = FacSpan;  //所有工厂的完工时间
        m_ArchiveFacECArray[PS] = FacEC;  //所有工厂的能耗
        m_ArchiveSpanArray[PS] = MS;  // 最大完工时间
        m_ArchiveTECArray[PS] = TEC;  //总能耗
        m_ArchiveSpeedVector[PS] = m_SpeedMatrix;

        //检查makespan和TEC
//        this->CheckSol(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], m_ArchiveSpanArray[PS]);
//        this->CheckSolTEC(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], m_ArchiveTECArray[PS]);

        //更新参考集
        m_QCCEANoRapidPopulation[PS].m_FacFamSeqArray = FacFamSeq;
        m_QCCEANoRapidPopulation[PS].m_JobSeqInFamArray = JobSeqInFam;
        m_QCCEANoRapidPopulation[PS].MS = MS;
        m_QCCEANoRapidPopulation[PS].TEC = TEC;
        m_QCCEANoRapidPopulation[PS].m_SpeedVector = m_SpeedMatrix;  //速度矩阵
    }

    // 初始化组种群 Initilize Familiy-sequence population, i.e., PS1
    for (int PS = 0; PS < m_PS1; PS++) {
        //将档案集的解赋给组序列
        m_FacFamSeqArray[PS] = m_ArchiveFacFamSeqArray[PS];
        m_SpanArray1[PS] = m_ArchiveSpanArray[PS];
        m_TECArray1[PS] = m_ArchiveTECArray[PS];
        m_Map1[PS] = PS;  //标记档案集中对应的jobpop序号
    }

    // 初始化工件种群 Initialize Job-Sequence population, i.e., PS2
    for (int PS = 0; PS < m_PS2; PS++) {
        //前AS个解
        m_JobSeqinFamArray[PS] = m_ArchiveJobSeqinFamArray[PS];
        m_SpanArray2[PS] = m_ArchiveSpanArray[PS];
        m_TECArray2[PS] = m_ArchiveTECArray[PS];
        m_Map2[PS] = PS;  //标记对应的档案集中fampop序号
    }

    //归一化
    Normalize(m_QCCEANoRapidPopulation, m_NadirPointMS, m_NadirPointTEC, m_IdealPointMS, m_IdealPointTEC);
    //cout << "nadirpointMS：" << m_NadirPointMS << "\tnadirpointTEC：" << m_NadirPointTEC << endl;

//    //计算收敛指标
    calc_convergence_ind(m_QCCEANoRapidPopulation);
//
//    //计算分散指标
    calc_distribution_ind(m_QCCEANoRapidPopulation);
//      cout<<"InitialPop end"<<endl;

}

void QCCEANoRapid::EvolutionProcess(int mu) {

    vector<vector<int>> FacFamSeq(this->m_Factories);
    vector<vector<int>> JobSeqInFam;

    thr_zeta = 1.0;
    long InitTime = Base::GetElapsedProcessTime();

    vector<vector<vector<double>>> Qfam(m_PS1, vector<vector<double>>(3, vector<double>(5)));//q表
    // 定义破坏算子和重构算子的权重
    vector<vector<double>> scores_combined(m_PS1, vector<double>(3, 1.0));
    // 破坏和重构算子的初始权重
    //进化
    while (Base::GetElapsedProcessTime() - InitTime < m_TimeLimit) {
        //协同进化--组
        cout<<"group"<<endl;
        for (int PS = 0; PS < m_PS1; PS++) {
            m_SpeedMatrix = m_ArchiveSpeedVector[PS];
            //得到真正的处理时间 及单位加工能耗
            for (int j = 0; j < m_Jobs; j++) {
                for (int i = 0; i < m_Machines; i++) {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }
            vector<vector<int>> FacFamSeq = m_FacFamSeqArray[PS];
            int Map1 = rand() % m_RefSize;// 随机挑选一个合作者（工件序列）从档案集中 form a solution by randomly selecting a RefJobSeq;
            vector<int> FacSpan(m_Factories, 0);
            vector<float> FacEC(m_Factories, 0);
            int ObjectMS;
            float ObjectTEC;
            GetMSandTECForPerandToalFac(FacFamSeq, m_ArchiveJobSeqinFamArray[Map1], FacSpan, FacEC, ObjectMS, ObjectTEC);
            int OrgMS = ObjectMS;
            float OrgTEC = ObjectTEC;
            int NewState = 0;
            int OldState = 0;

            //使用q-learning进行组序列选择策略
            int Act = DetermineAction(OldState, Qfam[PS]);
            PerformAction_Fams(Act, FacFamSeq, m_ArchiveJobSeqinFamArray[Map1],
                               m_QCCEANoRapidPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
            double Reward = CaculateReward(OrgMS, ObjectMS, OrgTEC, ObjectTEC, m_QCCEANoRapidPopulation);

            NewState = DetermineState(ObjectMS, ObjectTEC, m_QCCEANoRapidPopulation);
            UpdateQ(OldState, NewState, Act, Reward, Qfam[PS]);
            OldState = NewState;
        }
        cout<<"group end "<<endl;
        cout<<"job"<<endl;
        //协同进化--工件
        for (int PS = 0; PS < m_PS2; PS++) {
            m_SpeedMatrix = m_ArchiveSpeedVector[PS];
            //得到真正的处理时间 及单位加工能耗
            for (int j = 0; j < m_Jobs; j++) {
                for (int i = 0; i < m_Machines; i++) {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }

            vector<vector<int>> JobSeqInFam = m_JobSeqinFamArray[PS];
            int Map2 = rand() % m_RefSize;// 随机挑选一个合作者（组序列） form a solution by randomly selecting a ReFacFamseq
            vector<int> FacSpan(m_Factories, 0);
            vector<float> FacTEC(m_Factories, 0);
            int ObjectMS = -1;
            float ObjectTEC = -1;

            GetMSandTECForPerandToalFac(m_ArchiveFacFamSeqArray[Map2], JobSeqInFam, FacSpan, FacTEC,
                                             ObjectMS, ObjectTEC);

            JobLocalSearch(m_ArchiveFacFamSeqArray[Map2], JobSeqInFam,
                           m_NadirPointMS, m_NadirPointTEC, m_IdealPointMS, m_IdealPointTEC,
                           m_QCCEANoRapidPopulation, FacSpan, FacTEC, ObjectMS, ObjectTEC);
            BasedInd_SwapJob(m_ArchiveFacFamSeqArray[Map2], JobSeqInFam,
                                 m_NadirPointMS, m_NadirPointTEC, m_IdealPointMS, m_IdealPointTEC, m_QCCEANoRapidPopulation,
                                 FacSpan, FacTEC,
                                 ObjectMS, ObjectTEC);
        }
        cout<<"job end "<<endl;
        //协同进化--档案集
        cout<<"dangan "<<endl;
        for (int PS = 0; PS < m_RefSize; PS++) {

            m_SpeedMatrix = m_ArchiveSpeedVector[PS];
            //得到真正的处理时间 及单位加工能耗
            for (int j = 0; j < m_Jobs; j++) {
                for (int i = 0; i < m_Machines; i++) {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }
            vector<int> FacSpan(m_Factories, 0);
            vector<float> FacTEC(m_Factories, 0);
            int ObjectMS = -1;
            float ObjectTEC = -1;
            GetMSandTECForPerandToalFac(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], FacSpan, FacTEC,   ObjectMS, ObjectTEC);
            //对组和工件执行破坏重构操作
            BasedInd_RefDeconAndCon_New(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS],
                                        m_NadirPointMS, m_NadirPointTEC, m_IdealPointMS,
                                        m_IdealPointTEC, m_QCCEANoRapidPopulation, FacSpan, FacTEC,
                                        m_ArchiveSpanArray[PS], m_ArchiveTECArray[PS], scores_combined[PS]);

        }
        cout<<"dangan end "<<endl;
        //从参考集选择解 更新档案集，组，工件；并且变速
        UpdateArchiveGroupJobSet(mu);

    }

    //节能和降低MS
    SaveTECandDeMS(m_QCCEANoRapidPopulation);
    cout<<"Evoltuon end "<<endl;
}


int QCCEANoRapid::RunEvolution(int CPUTime, vector<vector<Individual>> &QCCEANoRapidFinalAfterRepParetoSet, int AN, int mu) {
    ReadInstanceFileNameList("..\\Benchmark\\");
    int Instances = 405;
    double T = 0.8;
    int Reps = 5; //针对每个算例重复运行的次数

    string FileDirectory = "..\\Result\\";
    ostringstream str;
    str << "QCCEANoRapid_" << CPUTime << "_experiment" << ".txt"; //不同的算法
    ofstream ofile;
    ofile.open(FileDirectory + str.str());

    for (int ins = 174; ins < Instances; ins++) {

        this->ReadInstance(ins);
        this->GetJobTotalPTime();
        this->GetJobWeightTotalPTime();
        this->GetFamSumSetupTime();
        this->GetFamTotalSkewness();
        this->GetFamTotalPTime_QCCEA();
        this->GetFamAvgSetupTime_QCCEA();
        this->GetFamTotalPTimeOnLastMachine_QCCEA();
        this->GetFamTotalPTimeOnFirstMachine_QCCEA();
        this->GetFamWeightTotalPTime_QCCEA();

        //新
        vector<Individual> FinalQCCEANoRapidParetoSet;
        vector<Individual> AfterRepParetoSet;

        vector<Individual> TempAfterRepParetoSet;
        TempAfterRepParetoSet.clear();

        for (int r = 0; r < Reps; r++) {

            long TimeLimit = CPUTime * m_Machines * m_Families; //original: 20 * m_Jobs * m_Machines
            this->SetParameters(AN,  TimeLimit, T);

            //变速后的处理时间（真正的处理时间）
            m_TureJobOpertime.clear();
            m_TureJobOpertime.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));

            //速度矩阵
            m_SpeedMatrix.clear();
            m_SpeedMatrix.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));

            //处理时间P的单位能耗
            UnitPEC.clear();
            UnitPEC.resize(this->m_Jobs, vector<float>(this->m_Machines, 0));


            MachReadyTime.clear();
            MachReadyTime.resize(this->m_Machines,0);

            //初始化种群
            InitialPop();

            //进化过程
            EvolutionProcess(mu);

            //非支配解
            m_QCCEANoRapidParetoSet.clear();

            Pareto_relation(m_QCCEANoRapidPopulation);

            //弱化支配解 de-emphasize dominated solutions
            for (int j = 0; j < m_QCCEANoRapidPopulation.size(); j++) {
                for (int i = 0; i < m_QCCEANoRapidPopulation.size(); i++) {
                    if (m_QCCEANoRapidPopulation[i].flag == 0) {
                        if (m_QCCEANoRapidPopulation[i].pareto_rel[j] == 1) {
                            m_QCCEANoRapidPopulation[i].flag = 999;
                        }
                    }
                }
            }

            for (int i = 0; i < m_QCCEANoRapidPopulation.size(); i++) {
                if (m_QCCEANoRapidPopulation[i].flag == 0)
                    m_QCCEANoRapidParetoSet.push_back(m_QCCEANoRapidPopulation[i]);
            }

            //去除重复
            FinalQCCEANoRapidParetoSet.clear();
            for (int i = 0; i < m_QCCEANoRapidParetoSet.size(); i++) {
                bool fg = true;
                for (int j = 0; j < FinalQCCEANoRapidParetoSet.size(); j++) {
                    if ((FinalQCCEANoRapidParetoSet[j].MS == m_QCCEANoRapidParetoSet[i].MS) &&
                        (FinalQCCEANoRapidParetoSet[j].TEC == m_QCCEANoRapidParetoSet[i].TEC)) {
                        fg = false;
                        break;
                    }
                }
                if (fg) {
                    FinalQCCEANoRapidParetoSet.push_back(m_QCCEANoRapidParetoSet[i]);
                }

            }

            long EndTime_IG = Base::GetElapsedProcessTime();

            for (int PS = 0; PS < FinalQCCEANoRapidParetoSet.size(); PS++) {
                TempAfterRepParetoSet.push_back(FinalQCCEANoRapidParetoSet[PS]);
            }
        }//end rep

        //非支配解
        AfterRepParetoSet.clear();

        Pareto_relation(TempAfterRepParetoSet);

        //弱化支配解 de-emphasize dominated solutions
        for (int j = 0; j < TempAfterRepParetoSet.size(); j++) {
            for (int i = 0; i < TempAfterRepParetoSet.size(); i++) {
                if (TempAfterRepParetoSet[i].flag == 0) {
                    if (TempAfterRepParetoSet[i].pareto_rel[j] == 1) {
                        TempAfterRepParetoSet[i].flag = 999;
                    }
                }
            }
        }

        for (int i = 0; i < TempAfterRepParetoSet.size(); i++) {
            if (TempAfterRepParetoSet[i].flag == 0)
                AfterRepParetoSet.push_back(TempAfterRepParetoSet[i]);
        }

        //去除重复
        vector<Individual> FinalAfterRepParetoSet;
        FinalAfterRepParetoSet.clear();
        for (int i = 0; i < AfterRepParetoSet.size(); i++) {
            bool fg = true;
            for (int j = 0; j < FinalAfterRepParetoSet.size(); j++) {
                if ((FinalAfterRepParetoSet[j].MS == AfterRepParetoSet[i].MS) &&
                    (FinalAfterRepParetoSet[j].TEC == AfterRepParetoSet[i].TEC)) {
                    fg = false;
                    break;
                }
            }
            if (fg) {
                FinalAfterRepParetoSet.push_back(AfterRepParetoSet[i]);
            }
        }

        QCCEANoRapidFinalAfterRepParetoSet[ins].clear();
        cout << ins + 1 << "\t" << "Factories:" << this->m_Factories << "\t" << "Machines :" << this->m_Machines
             << "\t" << "Families:" << this->m_Families << "\t" << "Jobs :" << this->m_Jobs << "\t";

        for (int PS = 0; PS < FinalAfterRepParetoSet.size(); PS++) {
            // 输出 MS 和 TEC，所有解都放在同一行
            cout << "MS:" << FinalAfterRepParetoSet[PS].MS
                 << " " << "TEC:" << FinalAfterRepParetoSet[PS].TEC << "\t";
            ofile << FinalAfterRepParetoSet[PS].MS << "  " << FinalAfterRepParetoSet[PS].TEC << "," << "\t";
            QCCEANoRapidFinalAfterRepParetoSet[ins].push_back(FinalAfterRepParetoSet[PS]);
        }
        ofile << endl;  // 每个实例结果结束后换行
        cout << endl;   // 控制台输出换行，仅仅是在实例结束时换行

    }//end ins
    ofile.close();
    return 0;
}
