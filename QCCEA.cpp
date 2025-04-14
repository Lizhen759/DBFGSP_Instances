#include "QCCEA.h"

#include <random>

QCCEA::QCCEA() {

}

QCCEA::~QCCEA() {
}

void
QCCEA::SetParameters(int AN, long TimeLimit, double T) {

    this->m_TimeLimit = TimeLimit;

    this->m_PopSize = 1000;
    this->m_RefSize = AN;
    this->m_PS1 = AN;
    this->m_PS2 = AN;
    this->thr_zeta = 1.0;

    this->m_T = T;

    m_QCCEAPopulation.clear();
    m_QCCEAPopulation.resize(m_RefSize);

    m_NadirPointMS = -1;
    m_NadirPointTEC = -1;

    m_IdealPointMS = -1;
    m_IdealPointTEC = -1;

    //������ Reference set
    m_ArchiveSpanArray.clear();//�������깤ʱ��
    m_ArchiveSpanArray.resize(m_RefSize);

    m_ArchiveTECArray.clear();//���������ܺ�
    m_ArchiveTECArray.resize(m_RefSize);

    m_ArchiveSpeedVector.clear();
    m_ArchiveSpeedVector.resize(m_RefSize);

    m_ArchiveFacFamSeqArray.clear();//������������
    m_ArchiveFacFamSeqArray.resize(m_RefSize);

    m_ArchiveJobSeqinFamArray.clear();//�������鹤��
    m_ArchiveJobSeqinFamArray.resize(m_RefSize);

    m_ArchiveFacSpanArray.clear();//�����깤ʱ��
    m_ArchiveFacSpanArray.resize(m_RefSize);

    m_ArchiveFacECArray.clear();//�����ܺ�
    m_ArchiveFacECArray.resize(m_RefSize);

    //����Ⱥ
    m_SpanArray1.clear();//�깤ʱ��
    m_SpanArray1.resize(m_PS1);

    m_TECArray1.clear(); //�ܺ�
    m_TECArray1.resize(m_PS1);
    
    m_Map1.clear();//Э�����ڵ��������±�
    m_Map1.resize(m_PS1);
    
    m_FacFamSeqArray.clear();//������
    m_FacFamSeqArray.resize(m_PS1);

    //������Ⱥ
    m_SpanArray2.clear();//�깤ʱ��
    m_SpanArray2.resize(m_PS2);

    m_TECArray2.clear(); //�ܺ�
    m_TECArray2.resize(m_PS2);
    

    m_Map2.clear();//Э�����ڵ��������±�
    m_Map2.resize(m_PS2);


    m_JobSeqinFamArray.clear();//��������
    m_JobSeqinFamArray.resize(m_PS2);
}


void QCCEA::BasedInd_RandFamInFacTobestPos_New(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                               vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                               int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                               vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS,
                                               float &ObjectTEC) {


    int CriFac = 0;
    int max = INT_MIN;
    //�ؼ�����
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

        auto it = find(FacFamSeq[CriFac].begin(), FacFamSeq[CriFac].end(), CurFam); // ���������
        int j = static_cast<int>(it - FacFamSeq[CriFac].begin()); // ǿ��ת��Ϊ int ����
        FacFamSeq[CriFac].erase(it);
        RefreshJFDTimeHierarchy_InFactory_Erase(FacFamSeq[CriFac], JobSeqInFam, JFDTime, Hierarchy,
                                                j, 0);
        BestFac = -1;
        BestPos = -1;
        Ind = this->FindBestPosToInsertFamForAllFac_Ind_New(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy, CurFam, BestFac,
                                                            BestPos, nadirpointMS, nadirpointTEC,
                                                            idealpointMS, idealpointTEC, FacSpan, FacEC, ObjectMS, ObjectTEC,
                                                            CCEAPopulation);
        FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
        RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[BestFac], JobSeqInFam, JFDTime, Hierarchy, BestPos,
                                                 0,JobSeqInFam[CurFam].size() - 1);

        vector<int> AfterInsertSpanFac(FacFamSeq.size(), 0);
        vector<float> AfterInsertECFac(FacFamSeq.size(), 0);
        int AfterInsertMS = -1;
        float AfterInsertTEC = -1;
        GetMSandTECForPerandToalFacByJFD(FacFamSeq, JobSeqInFam, JFDTime, AfterInsertSpanFac, AfterInsertECFac, AfterInsertMS, AfterInsertTEC);
        ObjectMS = AfterInsertMS;
        ObjectTEC = AfterInsertTEC;
    }
}

void QCCEA::BasedInd_SwapFam_New(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                 vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy, int nadirpointMS,
                                 float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                 vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC,
                                 int &ObjectMS, float &ObjectTEC) {

    int OrgMS = ObjectMS;
    float OrgTEC = ObjectTEC;

    //�ж��Ƿ��滻��͵�������
    if (OrgMS > nadirpointMS)
        nadirpointMS = OrgMS;
    if (OrgMS < idealpointMS)
        idealpointMS = OrgMS;

    //�ж��Ƿ��滻��͵�������
    if (OrgTEC > nadirpointTEC)
        nadirpointTEC = OrgTEC;
    if (OrgTEC < idealpointTEC)
        idealpointTEC = OrgTEC;

    //��һ��
    float OrgnormalMS = -1;
    float OrgnormalTEC = -1;

    OrgnormalMS = (static_cast<float>(OrgMS - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    OrgnormalTEC = ((OrgTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));

    //���㽻��֮ǰָ��Ind
    float OrgConvergence_ind = 0;
    OrgConvergence_ind += (OrgnormalMS - 1.0f) * (OrgnormalMS - 1.0f);
    OrgConvergence_ind += (OrgnormalMS - 1.0f) * (OrgnormalTEC - 1.0f);
    OrgConvergence_ind = sqrt(OrgConvergence_ind);
    OrgConvergence_ind = 1 / (OrgConvergence_ind + 1);


    float  OrgDistribution_ind = 0.0;
    Individual::calc_distribution_ind(CCEAPopulation);
    for (int i = 0; i < CCEAPopulation.size(); i++) {
        for (int j = i + 1; j < CCEAPopulation.size(); j++) {
            OrgDistribution_ind += CCEAPopulation[i].distribution_ind[j];
        }
    }
    OrgDistribution_ind = OrgDistribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

    float OrgCombined_ind= 0.5f*OrgConvergence_ind+ 0.5f*OrgDistribution_ind;

    for (int i = 0; i < m_Families / m_Factories; i++) {
        
        CopyJFDHierarchy(JFDTime, Hierarchy, m_tempJFDTime1, m_tempHierarchy1);
        int CriFac = -1;
        int max = -1;
        int OptFac = -1;
        int min = INT_MAX;
        //�ؼ�������makespan��С�Ĺ���
        for (int facIndex = 0; facIndex < FacFamSeq.size(); facIndex++) {
            if (max < FacSpan[facIndex]) {
                max = FacSpan[facIndex];
                CriFac = facIndex;
            }
            if (min > FacSpan[facIndex]) {
                min = FacSpan[facIndex];
                OptFac = facIndex;
            }
        }
        // ȷ���ؼ����������Ź�������ͬ
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
        RefreshJFDTimeHierarchy_InFactory_Erase(FacFamSeq[CriFac], JobSeqInFam, m_tempJFDTime1, m_tempHierarchy1, Pos1,
                                                0);
        if (FacFamSeq[OptFac].empty()) {
            break;
        }

        // ��ʼ�����������ͷֲ�
        std::random_device rd;  // ����豸
        std::mt19937 gen(rd()); // Mersenne Twister ���棬�������������
        std::uniform_int_distribution<> dis(0, FacFamSeq[OptFac].size() - 1);  // ���ȷֲ������ɷ�Χ�ڵ������

        // ���������
        int Pos2 = dis(gen);  // ʹ���������ͷֲ������������

        int Fam2 = FacFamSeq[OptFac][Pos2];
        FacFamSeq[OptFac].erase(FacFamSeq[OptFac].begin() + Pos2);
        RefreshJFDTimeHierarchy_InFactory_Erase(FacFamSeq[OptFac], JobSeqInFam, m_tempJFDTime1, m_tempHierarchy1, Pos2,
                                                0);


        FacFamSeq[CriFac].insert(FacFamSeq[CriFac].begin() + Pos1, Fam2);
        RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[CriFac], JobSeqInFam, m_tempJFDTime1, m_tempHierarchy1, Pos1,
                                                 0, JobSeqInFam[Fam2].size() - 1);
        FacSpan[CriFac] = GetSpanPerFacByJFD(FacFamSeq[CriFac], JobSeqInFam, m_tempJFDTime1);
        FacEC[CriFac] = GetTECForPerFacByJFD(FacFamSeq[CriFac], JobSeqInFam, m_tempJFDTime1);

        FacFamSeq[OptFac].insert(FacFamSeq[OptFac].begin() + Pos2, Fam1);
        RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[OptFac], JobSeqInFam, m_tempJFDTime1, m_tempHierarchy1,
                                                 Pos2, 0, JobSeqInFam[Fam1].size() - 1);
        FacSpan[OptFac] = GetSpanPerFacByJFD(FacFamSeq[OptFac], JobSeqInFam, m_tempJFDTime1);
        FacEC[OptFac] = GetTECForPerFacByJFD(FacFamSeq[OptFac], JobSeqInFam, m_tempJFDTime1);

        // ���㽻����� Makespan/TEC
        int AfterMakespan = *max_element(FacSpan.begin(), FacSpan.end());
        float AfterTEC = accumulate(FacEC.begin(), FacEC.end(), 0.0f);

        //�ж��Ƿ��滻��͵�������
        if (AfterMakespan > nadirpointMS)
            nadirpointMS = AfterMakespan;
        if (AfterMakespan < idealpointMS)
            idealpointMS = AfterMakespan;

        //�ж��Ƿ��滻��͵�������
        if (AfterTEC > nadirpointTEC)
            nadirpointTEC = AfterTEC;
        if (AfterTEC < idealpointTEC)
            idealpointTEC = AfterTEC;

        //�ж��Ƿ�Ȳ���ǰ�Ľ������Ľ������ο���
        bool flag = true;
        if ((AfterMakespan < OrgMS && AfterTEC <= OrgTEC) || (AfterMakespan <= OrgMS && AfterTEC < OrgTEC)) {
            for (auto & Ind : CCEAPopulation) {
                if ((AfterMakespan == Ind.MS) && (AfterTEC == Ind.TEC)) {
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
                tempIndi.m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���
                CCEAPopulation.push_back(tempIndi);
            }
        }

        //��һ��
        float AfternormalMS = -1;
        float AfternormalTEC = -1;


        AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) /static_cast<float>(nadirpointMS - idealpointMS));
        AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));

        //����ָ��Ind
        float Afterconvergence_ind = 0;
        Afterconvergence_ind += (AfternormalMS - 1.0f) * (AfternormalMS - 1.0f);
        Afterconvergence_ind += (AfternormalTEC - 1.0f) * (AfternormalTEC - 1.0f);
        Afterconvergence_ind = sqrt(Afterconvergence_ind);
        Afterconvergence_ind = 1 / (Afterconvergence_ind + 1);

        float  AfterDistribution_ind = 0.0;
        Individual::calc_distribution_ind(CCEAPopulation);
        for (int ind = 0; ind < CCEAPopulation.size(); ind++) {
            for (int j = ind + 1; j < CCEAPopulation.size(); j++) {
                AfterDistribution_ind += CCEAPopulation[ind].distribution_ind[j];
            }
        }
        AfterDistribution_ind = AfterDistribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

        float AfterCombined_ind = 0.5f * AfterDistribution_ind + 0.5f * Afterconvergence_ind;

        if (AfterCombined_ind < OrgCombined_ind) {
            ObjectMS = AfterMakespan;
            ObjectTEC = AfterTEC;
            CopyJFDHierarchy(m_tempJFDTime1, m_tempHierarchy1, JFDTime, Hierarchy);
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

void QCCEA::BasedIndLS_SetupSwap_New(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                     vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                     int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                     vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC,
                                     int &ObjectMS, float &ObjectTEC) {


    int OrgMS = ObjectMS;
    float OrgTEC = ObjectTEC;

    //�ж��Ƿ��滻��͵�������
    if (OrgMS > nadirpointMS)
        nadirpointMS = OrgMS;
    if (OrgMS < idealpointMS)
        idealpointMS = OrgMS;

    //�ж��Ƿ��滻��͵�������
    if (OrgTEC > nadirpointTEC)
        nadirpointTEC = OrgTEC;
    if (OrgTEC < idealpointTEC)
        idealpointTEC = OrgTEC;

    //��һ��
    float OrgnormalMS = -1;
    float OrgnormalTEC = -1;

    OrgnormalMS = (static_cast<float>(OrgMS - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    OrgnormalTEC = ((OrgTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));

    //���㽻��֮ǰָ��Ind
    float OrgConvergence_ind = 0;
    OrgConvergence_ind += (OrgnormalMS - 1.0f) * (OrgnormalMS - 1.0f);
    OrgConvergence_ind += (OrgnormalTEC - 1.0f) * (OrgnormalTEC - 1.0f);
    OrgConvergence_ind = sqrt(OrgConvergence_ind);
    OrgConvergence_ind = 1 / (OrgConvergence_ind + 1);
    //cout << "���������ָ�꣺" << OrgConvergence_ind << endl;

    float  OrgDistribution_ind = 0.0;
    Individual::calc_distribution_ind(CCEAPopulation);
    for (int i = 0; i < CCEAPopulation.size(); i++) {
        for (int j = i + 1; j < CCEAPopulation.size(); j++) {
            OrgDistribution_ind += CCEAPopulation[i].distribution_ind[j];
        }
    }
    OrgDistribution_ind = OrgDistribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

    float OrgCombined_ind= 0.5f*OrgConvergence_ind+0.5f*OrgDistribution_ind;


    for (int i = 0; i < m_Families / m_Factories; i++) {
        CopyJFDHierarchy(JFDTime, Hierarchy, m_tempJFDTime1, m_tempHierarchy1);
        int CriMS = 0;
        int OptMS = INT_MAX;
        int CriFac = 0, OptFac = 0;
        //�ؼ�����
        for (int facIndex = 0; facIndex < FacFamSeq.size(); facIndex++) {
            if (CriMS < FacSpan[facIndex]) {
                CriFac = facIndex;  // �ؼ�����
                CriMS = FacSpan[facIndex];  // �������ֵ
            }
            if (OptMS >= FacSpan[facIndex])  // �ų��ؼ�����
            {
                OptFac = facIndex;  // ���Ź���
                OptMS = FacSpan[facIndex];  // ������Сֵ
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
        RefreshJFDTimeHierarchy_InFactory_Erase(FacFamSeq[CriFac], JobSeqInFam, m_tempJFDTime1, m_tempHierarchy1, Pos1,
                                                0);

        FacFamSeq[OptFac].erase(FacFamSeq[OptFac].begin() + Pos2);
        RefreshJFDTimeHierarchy_InFactory_Erase(FacFamSeq[OptFac], JobSeqInFam, m_tempJFDTime1, m_tempHierarchy1, Pos2,
                                                0);

        //��������crifac����������
        FacFamSeq[CriFac].insert(FacFamSeq[CriFac].begin() + Pos1, NowFam2);
        RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[CriFac], JobSeqInFam, m_tempJFDTime1, m_tempHierarchy1, Pos1,
                                                 0,
                                                 JobSeqInFam[NowFam2].size() - 1);
        FacSpan[CriFac] = GetSpanPerFacByJFD(FacFamSeq[CriFac], JobSeqInFam, m_tempJFDTime1);
        FacEC[CriFac] = GetTECForPerFacByJFD(FacFamSeq[CriFac], JobSeqInFam, m_tempJFDTime1);

        //��������optfac����������
        FacFamSeq[OptFac].insert(FacFamSeq[OptFac].begin() + Pos2, NowFam1);
        RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[OptFac], JobSeqInFam, m_tempJFDTime1, m_tempHierarchy1, Pos2,
                                                 0,
                                                 JobSeqInFam[NowFam1].size() - 1);
        FacSpan[OptFac] = GetSpanPerFacByJFD(FacFamSeq[OptFac], JobSeqInFam, m_tempJFDTime1);
        FacEC[OptFac] = GetTECForPerFacByJFD(FacFamSeq[OptFac], JobSeqInFam, m_tempJFDTime1);

        // ���㽻����� Makespan/TEC
        int AfterMakespan = *max_element(FacSpan.begin(), FacSpan.end());
        float AfterTEC = accumulate(FacEC.begin(), FacEC.end(), 0.0f);
        //�ж��Ƿ��滻��͵�������
        if (AfterMakespan > nadirpointMS)
            nadirpointMS = AfterMakespan;
        if (AfterMakespan < idealpointMS)
            idealpointMS = AfterMakespan;
        if (AfterTEC > nadirpointTEC)
            nadirpointTEC = AfterTEC;
        if (AfterTEC < idealpointTEC)
            idealpointTEC = AfterTEC;

        //�ж��Ƿ�Ȳ���ǰ�Ľ������Ľ������ο���
        bool flag = true;
        if (((AfterMakespan < OrgMS) && (AfterTEC < OrgTEC)) || ((AfterMakespan < OrgMS) && (AfterTEC == OrgTEC)) ||
            ((AfterMakespan == OrgMS) && (AfterTEC < OrgTEC)) || (AfterMakespan < idealpointMS) ||
            (AfterTEC < idealpointTEC)) {
            for (auto & ind : CCEAPopulation) {
                if ((AfterMakespan == ind.MS) && (AfterTEC == ind.TEC)) {
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
                tempIndi.m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���
                CCEAPopulation.push_back(tempIndi);
            }
        }

        //��һ��
        float AfternormalMS = -1;
        float AfternormalTEC = -1;
        //cout << endl << "normalize" << endl;

        AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) /
                         static_cast<float>(nadirpointMS - idealpointMS));
        AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
        //cout << "�����һ�����MS��" << AfternormalMS << "\tTEC��" << AfternormalTEC << endl;

        //����ָ��Ind
        float Afterconvergence_ind = 0;
        Afterconvergence_ind += (AfternormalMS - 1.0f) * (AfternormalMS - 1.0f);
        Afterconvergence_ind += (AfternormalTEC - 1.0f) * (AfternormalTEC - 1.0f);
        Afterconvergence_ind = sqrt(Afterconvergence_ind);
        Afterconvergence_ind = 1 / (Afterconvergence_ind + 1);

        float  AfterDistribution_ind = 0;
        Individual::calc_distribution_ind(CCEAPopulation);
        for (int ind = 0; ind < CCEAPopulation.size(); ind++) {
            for (int j = ind + 1; j < CCEAPopulation.size(); j++) {
                AfterDistribution_ind += CCEAPopulation[ind].distribution_ind[j];
            }
        }
        AfterDistribution_ind = AfterDistribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

        float AfterCombined_ind = 0.5f * AfterDistribution_ind + 0.5f * Afterconvergence_ind;


        if (AfterCombined_ind < OrgCombined_ind) {
            ObjectMS = AfterMakespan;
            ObjectTEC = AfterTEC;
            CopyJFDHierarchy(m_tempJFDTime1, m_tempHierarchy1, JFDTime, Hierarchy);
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

void QCCEA::BasedInd_ShiftFam_New(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                  vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                  int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                  vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC,
                                  int &ObjectMS, float &ObjectTEC) {
    int OrgMS = ObjectMS;
    float OrgTEC = ObjectTEC;

    //�ж��Ƿ��滻��͵�������
    if (OrgMS > nadirpointMS)
        nadirpointMS = OrgMS;
    if (OrgMS < idealpointMS)
        idealpointMS = OrgMS;

    //�ж��Ƿ��滻��͵�������
    if (OrgTEC > nadirpointTEC)
        nadirpointTEC = OrgTEC;
    if (OrgTEC < idealpointTEC)
        idealpointTEC = OrgTEC;

    //��һ��
    float OrgnormalMS = -1;
    float OrgnormalTEC = -1;

    OrgnormalMS = (static_cast<float>(OrgMS - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    OrgnormalTEC = ((OrgTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));

    //���㽻��֮ǰָ��Ind
    float OrgConvergence_ind = 0;
    OrgConvergence_ind += (OrgnormalMS - 1.0f) * (OrgnormalMS - 1.0f);
    OrgConvergence_ind += (OrgnormalTEC - 1.0f) * (OrgnormalTEC - 1.0f);
    OrgConvergence_ind = sqrt(OrgConvergence_ind);
    OrgConvergence_ind = 1 / (OrgConvergence_ind + 1);
    //cout << "���������ָ�꣺" << OrgConvergence_ind << endl;

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

    CopyJFDHierarchy(JFDTime, Hierarchy, m_tempJFDTime1, m_tempHierarchy1);

    int Old_Span1 =FacSpan[facIndex1];
    float Old_EC1 = FacEC[facIndex1] ;
    int Old_Span2=FacSpan[facIndex2] ;
    float Old_EC2 =FacEC[facIndex2]  ;

    int CurFam = FacFamSeq[facIndex1][pos];
    FacFamSeq[facIndex1].erase(FacFamSeq[facIndex1].begin() + pos);
    RefreshJFDTimeHierarchy_InFactory_Erase(FacFamSeq[facIndex1], JobSeqInFam, m_tempJFDTime1, m_tempHierarchy1, pos,
                                            0);
    FacSpan[facIndex1] = GetSpanPerFacByJFD(FacFamSeq[facIndex1], JobSeqInFam, m_tempJFDTime1);
    FacEC[facIndex1] = GetTECForPerFacByJFD(FacFamSeq[facIndex1], JobSeqInFam, m_tempJFDTime1);


    int bestPos = -1;
    float Ind = FindBestPosToInsertFamForPerFac_Ind_New(facIndex2, FacFamSeq, FacFamSeq[facIndex2], JobSeqInFam,
                                                        m_tempJFDTime1,
                                                        m_tempHierarchy1, CurFam, bestPos,
                                                        nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, FacSpan, FacEC, OrgMS,
                                                        OrgTEC, CCEAPopulation);

    FacFamSeq[facIndex2].insert(FacFamSeq[facIndex2].begin() + bestPos, CurFam);
    RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[facIndex2], JobSeqInFam, m_tempJFDTime1, m_tempHierarchy1,
                                             bestPos, 0,JobSeqInFam[CurFam].size() - 1);

    FacSpan[facIndex2] = GetSpanPerFacByJFD(FacFamSeq[facIndex2], JobSeqInFam, m_tempJFDTime1);
    FacEC[facIndex2] = GetTECForPerFacByJFD(FacFamSeq[facIndex2], JobSeqInFam, m_tempJFDTime1);

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

    //��һ��
    float AfternormalMS = -1;
    float AfternormalTEC = -1;

    AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) /static_cast<float>(nadirpointMS - idealpointMS));
    AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));

    //����ָ��Ind
    float Afterconvergence_ind = 0;
    Afterconvergence_ind += (AfternormalMS - 1.0f) * (AfternormalMS - 1.0f);
    Afterconvergence_ind += (AfternormalTEC - 1.0f) * (AfternormalTEC - 1.0f);
    Afterconvergence_ind = sqrt(Afterconvergence_ind);
    Afterconvergence_ind = 1 / (Afterconvergence_ind + 1);

    float  AfterDistribution_ind = 0.0;
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
        CopyJFDHierarchy(m_tempJFDTime1, m_tempHierarchy1, JFDTime, Hierarchy);
    } else {
        // ���û�и��ƣ���Ҫ�ع�
        FacFamSeq[facIndex2].erase(FacFamSeq[facIndex2].begin() + bestPos);
        FacFamSeq[facIndex1].insert(FacFamSeq[facIndex1].begin() + pos, CurFam);
        FacSpan[facIndex1]= Old_Span1 ;
        FacEC[facIndex1] = Old_EC1 ;
        FacSpan[facIndex2]= Old_Span2 ;
        FacEC[facIndex2] = Old_EC2 ;
    }

}


void QCCEA::BasedInd_DeAndConFams_New(vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                      vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
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

    size_t Len = FacFamSeq[CriFac].size();
    vector<int> ExtractFamSeq;
    bool groupSelected = false;
    do {
        int Fac;
        if (ExtractFamSeq.size() < Len / 2) //�ӹؼ�������ѡd/2����
            Fac = CriFac;
        else  //���������ѡʣ�µ�d/2����
            Fac = rand() % m_Factories;
        if (FacFamSeq[Fac].size() > 1) {
            int Pos = rand() % FacFamSeq[Fac].size();
            ExtractFamSeq.push_back(FacFamSeq[Fac][Pos]);  //����ѡ��������ӵ�ExtractFamSeq
            FacFamSeq[Fac].erase(FacFamSeq[Fac].begin() + Pos);  //��ԭ������ɾ��posλ�õ���
            groupSelected = true;
        }

    } while (ExtractFamSeq.size() < Len && groupSelected);

    GetJFDHierarchy_Forward_New( FacFamSeq,JobSeqInFam,JFDTime,Hierarchy);

    // ����ѡ������d������뵽��õ�λ�� Insert the extracted Families into the best Positions
    for (int Fam = 0; Fam < ExtractFamSeq.size(); Fam++) {

        int CurFam = ExtractFamSeq[Fam];  //��ǰ��
        int bestFac = -1, bestPos = -1;  //��õĹ�����λ��
        if (Fam == ExtractFamSeq.size() - 1) {
            float Ind = this->FindBestPosToInsertFamForAllFac_Ind_DR_New(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy,
                                                                         CurFam, bestFac, bestPos, nadirpointMS,
                                                                         nadirpointTEC, idealpointMS, idealpointTEC, CCEAPopulation);
        } else {
            float randValue = wyt_rand_include_right(0.0f,1.0f);
            if (randValue < 0.5) {

                FindBestPosToInsertFamForAllFacs_MakeSpan_New(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy, CurFam,
                                                              bestFac, bestPos);

            } else {
                FindBestPosInsertFamAllFacs_TEC_New(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy, CurFam, bestFac,
                                                    bestPos);

            }
        }
        // ����õĹ��������λ���в��뵱ǰ�� Insert CurFam to bestPos at bestFac
        FacFamSeq[bestFac].insert(FacFamSeq[bestFac].begin() + bestPos, CurFam);
        RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[bestFac], JobSeqInFam, JFDTime, Hierarchy, bestPos, 0,
                                                 JobSeqInFam[CurFam].size() - 1);
    }

}

void QCCEA::BasedInd_SwapJob_New(const vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                 vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                 int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                 vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC,
                                 int &ObjectMS, float &ObjectTEC) {

    int OrgMS = ObjectMS;
    float OrgTEC = ObjectTEC;

    //�ж��Ƿ��滻��͵�������
    if (OrgMS > nadirpointMS)
        nadirpointMS = OrgMS;
    if (OrgMS < idealpointMS)
        idealpointMS = OrgMS;

    //�ж��Ƿ��滻��͵�������
    if (OrgTEC > nadirpointTEC)
        nadirpointTEC = OrgTEC;
    if (OrgTEC < idealpointTEC)
        idealpointTEC = OrgTEC;

    //��һ��
    float OrgnormalMS = -1;
    float OrgnormalTEC = -1;
    //cout << endl << "normalize" << endl;

    OrgnormalMS = (static_cast<float>(OrgMS - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    OrgnormalTEC = ((OrgTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
    //cout << "ԭʼ�����һ�����MS��" << OrgnormalMS << "\tTEC��" << OrgnormalTEC << endl;

    //���㽻��֮ǰָ��Ind
    float OrgConvergence_ind = 0;
    OrgConvergence_ind += (OrgnormalMS - 1.0f) * (OrgnormalMS - 1.0f);
    OrgConvergence_ind += (OrgnormalTEC - 1.0f) * (OrgnormalTEC - 1.0f);
    OrgConvergence_ind = sqrt(OrgConvergence_ind);
    OrgConvergence_ind = 1 / (OrgConvergence_ind + 1);
    //cout << "���������ָ�꣺" << OrgConvergence_ind << endl;

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

    //�ؼ�����
    for (int i = 0; i < FacFamSeq.size(); i++) {
        if (min < FacSpan[i]) {
            min = FacSpan[i];
            fac = i;
        }

    }

    for (int Fam = 0; Fam < FacFamSeq[fac].size(); Fam++) {
        int CurFam = FacFamSeq[fac][Fam]; //��ǰ��
        if (JobSeqInFam[CurFam].size() > 1) //�鹤��������1�����оֲ�����
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

            CopyJFDHierarchy(JFDTime, Hierarchy, m_tempJFDTime1, m_tempHierarchy1);
            vector<vector<int>> TempJobSeqInFam = JobSeqInFam;

            auto it = find(begin(FacFamSeq[fac]), end(FacFamSeq[fac]), CurFam);
            if (it != end(FacFamSeq[fac])) {
                std::vector<int>::size_type famPos = it - begin(FacFamSeq[fac]);  // ʹ����ȷ������
                int Job = JobSeqInFam[CurFam][pt1];
                JobSeqInFam[CurFam][pt1] = JobSeqInFam[CurFam][pt2];
                JobSeqInFam[CurFam][pt2] = Job;
            }

            GetJFDHierarchy_Forward_New(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy);

            FacSpan[fac] = GetSpanPerFacByJFD(FacFamSeq[fac], JobSeqInFam, JFDTime);

            int AfterMakespan = *max_element(FacSpan.begin(), FacSpan.end());

            //�ж��Ƿ��滻��͵�������
            if (AfterMakespan > nadirpointMS)
                nadirpointMS = AfterMakespan;
            if (AfterMakespan < idealpointMS)
                idealpointMS = AfterMakespan;


            float AfterTEC = 0;
            FacEC[fac] = GetTECForPerFacByJFD(FacFamSeq[fac], JobSeqInFam, JFDTime);
            //cout << "�ؼ�����" << fac << "���ܺģ�" << FacEC[fac] << endl;

            AfterTEC = accumulate(FacEC.begin(), FacEC.end(), 0.0f);
            //cout << "���ܺģ�" << AfterTEC << endl;

            //�ж��Ƿ��滻��͵�������
            if (AfterTEC > nadirpointTEC)
                nadirpointTEC = AfterTEC;
            if (AfterTEC < idealpointTEC)
                idealpointTEC = AfterTEC;

            //�ж��Ƿ�Ȳ���ǰ�Ľ������Ľ������ο���
            bool flag = true;
            if (((AfterMakespan < OrgMS) && (AfterTEC < OrgTEC)) || ((AfterMakespan < OrgMS) && (AfterTEC == OrgTEC)) ||
                ((AfterMakespan == OrgMS) && (AfterTEC < OrgTEC)) || (AfterMakespan < idealpointMS) ||
                (AfterTEC < idealpointTEC)) {
                for (auto & i : CCEAPopulation) {
                    if ((AfterMakespan == i.MS) && (AfterTEC == i.TEC)) {
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
                    tempIndi.m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���
                    CCEAPopulation.push_back(tempIndi);
                }
            }

            //��һ��
            float AfternormalMS = -1;
            float AfternormalTEC = -1;
            //cout << endl << "normalize" << endl;

            AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) /
                             static_cast<float>(nadirpointMS - idealpointMS));
            AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
            //cout << "�����һ�����MS��" << AfternormalMS << "\tTEC��" << AfternormalTEC << endl;

            //����ָ��Ind
            float Afterconvergence_ind = 0;
            Afterconvergence_ind += (AfternormalMS - 1.0f) * (AfternormalMS - 1.0f);
            Afterconvergence_ind += (AfternormalTEC - 1.0f) * (AfternormalTEC - 1.0f);
            Afterconvergence_ind = sqrt(Afterconvergence_ind);
            Afterconvergence_ind = 1 / (Afterconvergence_ind + 1);
            //cout << "������ĸ��������ָ�꣺" << Afterconvergence_ind << endl;

            float  AfterDistribution_ind = 0.0;
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
                CopyJFDHierarchy(m_tempJFDTime1, m_tempHierarchy1, JFDTime, Hierarchy);
            }
        }
        //��麯��
        //CheckSol(FacFamSeq, JobSeqInFam, ObjectMS);
        //CheckSolTEC(FacFamSeq, JobSeqInFam, ObjectTEC);
    }
}

void QCCEA::JobLocalSearch(const vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                           vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                           int NadirPointMS, float NadirPointTEC, int IdealPointMS, float IdealPointTEC,
                           vector<Individual> &CCEAPopulation, vector<int>& FacSpan, vector<float>& FacEC, int &ObjectMS, float &ObjectTEC)
{

    for (int FacIndex = 0; FacIndex < m_Factories; ++FacIndex)
    {
        for (int FamIndex = 0; FamIndex < FacFamSeq[FacIndex].size(); ++FamIndex)
        {
            int CurFam = FacFamSeq[FacIndex][FamIndex];
            //�������ֻ��һ�������������
            if (JobSeqInFam[CurFam].size() == 1)
            {
                continue;
            }

            vector<int> SeqForExtracting = JobSeqInFam[CurFam];//��������
            shuffle(SeqForExtracting.begin(), SeqForExtracting.end(), rand_generator());//�����������
            int Counter = 0;
            int Pos = 0;
            while (Counter < ceil(SeqForExtracting.size() / 2))
            {

                vector<int> tempFacSpan(this->m_Factories, 0);
                vector<int> tempFacEC(this->m_Factories, 0);
                int CurJob = SeqForExtracting[Pos];
                auto it = find(JobSeqInFam[CurFam].begin(), JobSeqInFam[CurFam].end(), CurJob);
                int OrgPos = static_cast<int>(std::distance(JobSeqInFam[CurFam].begin(), it));

                JobSeqInFam[CurFam].erase(it);
                RefreshJFDTimeHierarchy_InFactory_Erase(FacFamSeq[FacIndex], JobSeqInFam, JFDTime,
                                                                     Hierarchy, FamIndex, OrgPos);
                int bestPos = -1;
                auto re = FindBestPosToInsertJobForPerFac_Ind_New(FacIndex,FacFamSeq,FacFamSeq[FacIndex], JobSeqInFam, JFDTime, Hierarchy,
                                                             CurFam,FamIndex, CurJob, bestPos, NadirPointMS,
                                                            NadirPointTEC, IdealPointMS,  IdealPointTEC,ObjectMS, ObjectTEC, CCEAPopulation );

                JobSeqInFam[CurFam].insert(JobSeqInFam[CurFam].begin() + bestPos, CurJob);
                RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[FacIndex], JobSeqInFam, JFDTime,
                                                                          Hierarchy, FamIndex, bestPos, bestPos);
                FacSpan[FacIndex] = GetSpanPerFacByJFD(FacFamSeq[FacIndex], JobSeqInFam,JFDTime);
                FacEC[FacIndex] = GetTECForPerFacByJFD(FacFamSeq[FacIndex], JobSeqInFam,JFDTime);

                Pos = (Pos + 1) % (int) ceil(JobSeqInFam[CurFam].size() / 2);
                Counter += 1;
            }
        }
    }
}

void QCCEA::BasedInd_RefDeconAndCon_New(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                        vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                        int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                        vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC,
                                        int &ObjectMS, float &ObjectTEC, vector<double> &scores_combined) {

    int OrgMS = ObjectMS;
    float OrgTEC = ObjectTEC;

    //�ж��Ƿ��滻��͵�������
    if (OrgMS > nadirpointMS)
        nadirpointMS = OrgMS;
    if (OrgMS < idealpointMS)
        idealpointMS = OrgMS;

    //�ж��Ƿ��滻��͵�������
    if (OrgTEC > nadirpointTEC)
        nadirpointTEC = OrgTEC;
    if (OrgTEC < idealpointTEC)
        idealpointTEC = OrgTEC;

    // ��һ��
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
    CopyJFDHierarchy(JFDTime, Hierarchy, m_tempJFDTime1, m_tempHierarchy1);
    vector<int> FamsExtracted;
    unordered_map<int, vector<int>> JobsExtracted;

    // ���̶�ѡ���ƻ�+�ع����
    auto select_operator = [](const vector<double> &weights) -> int {
        double total_weight = accumulate(weights.begin(), weights.end(), 0.0);  // ʹ�� double ���ͱ��⾫�ȶ�ʧ
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, total_weight);
        double r = dis(gen);  // ʹ�� double ����������һ��
        double cumulative = 0.0;  // ʹ�� double ����
        for (size_t i = 0; i < weights.size(); ++i) {  // ʹ�� size_t �������������Ͳ�ƥ��
            cumulative += weights[i];
            if (r <= cumulative) {
                return static_cast<int>(i);  // ȷ������ֵ�� int ����
            }
        }
        return 0;  // Ĭ�Ϸ���ֵ
    };


    // ѡ���ƻ��ع��������
    int DeconConIndex = select_operator(scores_combined);

    // ִ�ж�Ӧ���ƻ��ع�����
    ApplyDeconAndCon_Ref(DeconConIndex, FacFamSeq, JobSeqInFam, JFDTime, Hierarchy, FamsExtracted,
                         JobsExtracted,nadirpointMS, nadirpointTEC,
                         idealpointMS, idealpointTEC, FacSpan, FacEC,CCEAPopulation);

    // ������Ŀ��ֵ
    int AfterMakespan = GetSpanForAllFacByJFD(FacFamSeq, JobSeqInFam, JFDTime);
    float AfterTEC = GetTECForAllFacByJFD(FacFamSeq, JobSeqInFam, JFDTime);

    //�ж��Ƿ��滻��͵�������
    if (AfterTEC > nadirpointTEC)
        nadirpointTEC = AfterTEC;
    if (AfterTEC < idealpointTEC)
        idealpointTEC = AfterTEC;

    //�ж��½��Ƿ����ڵ�ǰ�⣬���Ľ������ο���
    bool flag = true;
    if (((AfterMakespan < OrgMS) && (AfterTEC < OrgTEC)) ||
        ((AfterMakespan < OrgMS) && (AfterTEC == OrgTEC)) ||
        ((AfterMakespan == OrgMS) && (AfterTEC < OrgTEC)) ||
        (AfterMakespan < idealpointMS) || (AfterTEC < idealpointTEC)) {

        // ����½��Ƿ��Ѵ����ڲο���
        for (const auto &individual: CCEAPopulation) {
            if (AfterMakespan == individual.MS && AfterTEC == individual.TEC) {
                flag = false;
                break;
            }
        }

        if (flag) {
            // �½�δ�ڲο����У���ӵ��ο���
            Individual tempIndi;
            tempIndi.m_FacFamSeqArray = FacFamSeq;
            tempIndi.m_JobSeqInFamArray = JobSeqInFam;
            tempIndi.MS = AfterMakespan;
            tempIndi.TEC = AfterTEC;
            tempIndi.m_SpeedVector = m_SpeedMatrix;  // �ٶȾ���
            CCEAPopulation.push_back(tempIndi);
            // ����ʽ�ӷ�
            scores_combined[DeconConIndex] += 1.5; // �ƻ����Ӽӷ�
            scores_combined[DeconConIndex] += 1.5;  // �ع����Ӽӷ�
        } else {
            scores_combined[DeconConIndex] += 1.0; // �ƻ����Ӽӷ�
            scores_combined[DeconConIndex] += 1.0;  // �ع����Ӽӷ�
        }
    } else {
        // �½�δ�Ľ���ִ��������������ģ���˻�Ľ���׼��
        double Temperature = this->m_T * accumulate(this->m_JobTotalPTime.begin(), this->m_JobTotalPTime.end(), 0) /
                             (this->m_Factories * this->m_Jobs * this->m_Machines);
        bool accept = SA(OrgMS, OrgTEC, AfterMakespan, AfterTEC, Temperature);
        if (accept) {
            // ���ܽϲ�⣬����ϵ͵ļӷ�
            scores_combined[DeconConIndex] += 0.5; // �ƻ����Ӽӷ�
            scores_combined[DeconConIndex] += 0.5;  // �ع����Ӽӷ�
        } else {
            // �����ܽϲ�⣬���ӷ�
            scores_combined[DeconConIndex] += 0.0; // �ƻ����Ӽӷ�
            scores_combined[DeconConIndex] += 0.0;  // �ع����Ӽӷ�
        }
    }

    //��һ��
    float AfternormalMS = -1;
    float AfternormalTEC = -1;

    AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) /static_cast<float>(nadirpointMS - idealpointMS));
    AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));


    //����ָ��Ind
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
        CopyJFDHierarchy(m_tempJFDTime1, m_tempHierarchy1, JFDTime, Hierarchy);
    }

}

// temperature����ǰ�¶�
// ����true��ʾ�����½⣬����false��ʾ������
bool QCCEA::SA(int OrgMS, float OrgTEC, int AfterMS, float AfterTEC, double Temperature) {

    // Metropolis ׼��
    int deltaMS = OrgMS - AfterMS;
    float deltaTEC = OrgTEC - AfterTEC;
    float delta = static_cast<float>(deltaMS) + deltaTEC;
    double r1 = wyt_rand_include_right(0.1f, 1.0f);
    if (r1 < exp(delta) / Temperature) {
        return true;
    }
    return false;
}

void QCCEA::ApplyDeconAndCon_Ref(int DeAndConIndex, vector<vector<int>> &FacFamSeq,
                                 vector<vector<int>> &JobSeqInFam,
                                 vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                 vector<int> &FamsExtracted, unordered_map<int, vector<int>> &JobsExtracted,
                                 int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                 vector<int> &FacSpan, vector<float> &FacEC,vector<Individual> &CCEAPopulation) {
    // ����ѡ�������ִ���ƻ����ع�����
    if (DeAndConIndex == 0) {
//        cout<<"00"<<endl;
        BasedInd_Destruction_FamsAndJobs_New_Opeator0(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy, FamsExtracted,
                                                      JobsExtracted);
        BasedInd_Construction_FamsAndJobs_New_Opeator0(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy, FamsExtracted,
                                                       JobsExtracted,
                                                       nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC,CCEAPopulation);
//        cout<<"00 end "<<endl;
    } else if (DeAndConIndex== 1) {
//        cout<<"11"<<endl;
        BasedInd_Destruction_FamsAndJobs_New_Opeator1(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy, FamsExtracted,
                                                      JobsExtracted, FacSpan, FacEC);
        BasedInd_Construction_FamsAndJobs_New_Opeator1(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy, FamsExtracted,
                                                       JobsExtracted,
                                                       nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC,CCEAPopulation);
//        cout<<"11 end "<<endl;
    } else {
//        cout<<"22"<<endl;
        BasedInd_Destruction_FamsAndJobs_New_Opeator2(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy, FamsExtracted,
                                                      JobsExtracted, FacSpan, FacEC);
        BasedInd_Construction_FamsAndJobs_New_Opeator2(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy, FamsExtracted,
                                                       JobsExtracted,
                                                       nadirpointMS, nadirpointTEC, idealpointMS,
                                                       idealpointTEC,CCEAPopulation);
//        cout<<"22 end"<<endl;
    }
}

void
QCCEA::BasedInd_Destruction_FamsAndJobs_New_Opeator0(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                     vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
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

    // �ȼ���Ƿ������� 2 ��Ĺ���
    std::vector<int> availableFactories;
    for (int i = 0; i < this->m_Factories; ++i) {
        if (FacFamSeq[i].size() > 1) {
            availableFactories.push_back(i);  // ֻ�洢���� 2 ��Ĺ���
        }
    }

    // ���û�з��������Ĺ�����ֱ�ӷ���
    if (availableFactories.empty()) return;

    // ���ѡ��һ�����ʵĹ���
    int Fac = availableFactories[rand() % availableFactories.size()];

    // ���� D_extract �������ù���������
    D_extract = std::min(D_extract, (int) FacFamSeq[Fac].size());

    while (FamsExtracted.size() < D_extract) {
        int Pos = rand() % FacFamSeq[Fac].size();
        int FamExt = FacFamSeq[Fac][Pos];

        FacFamSeq[Fac].erase(FacFamSeq[Fac].begin() + Pos);
        FamsExtracted.push_back(FamExt);
        UpdateFamErasedFromFac(Pos, Fac);
    }


    //�ӹ�����������ȡ����
    for (int Fam : FamsExtracted) //����������Ĺ�����ȡһ��
    {
        if (JobSeqInFam[Fam].size() >= 3) {
            for (int j = 0; j <= JobSeqInFam[Fam].size() / 2; j++) {
                int JobPos = rand() % JobSeqInFam[Fam].size();
                int Job = JobSeqInFam[Fam][JobPos];
                JobSeqInFam[Fam].erase(JobSeqInFam[Fam].begin() + JobPos);
                JobsExtracted[Fam].push_back(Job);
            }
        }
    }
    GetJFDHierarchy_Forward_New(FacFamSeq,JobSeqInFam,JFDTime,Hierarchy);
}


void
QCCEA::BasedInd_Destruction_FamsAndJobs_New_Opeator1(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                     vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
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
    float max_val = -1.0f;
    for (int i = 0; i < FacFamSeq.size(); i++) {
        float p = wyt_rand_include_right(0.0f, 1.0f);
        if (p < 0.5) {
            if (max_val < static_cast<float>(FacSpan[i])) {  // ǿ��ת��
                max_val = static_cast<float>(FacSpan[i]);
                CriFac = i;  // ����� i �� int ���ͣ�������� narrowing conversion
            }
        } else {

            if (max_val < FacEC[i]) {
                max_val = FacEC[i];
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
    GetJFDHierarchy_Forward_New(FacFamSeq,JobSeqInFam,JFDTime,Hierarchy);
}


void
QCCEA::BasedInd_Destruction_FamsAndJobs_New_Opeator2(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                     vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
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

    int CriFac = 0; // Ĭ��ȡ��һ������������δ��ʼ������
    float max_val = -1.0f;  // ���� FacSpan �� FacEC ������С�� 0
    for (int i = 0; i < FacFamSeq.size(); i++) {
        float p = wyt_rand_include_right(0.0f, 1.0f);

        // �� FacSpan[i] ת��Ϊ float ���бȽ�
        if (p < 0.5) {
            if (max_val < static_cast<float>(FacSpan[i])) {  // ǿ��ת��
                max_val = static_cast<float>(FacSpan[i]);
                CriFac = i;  // ����� i �� int ���ͣ�������� narrowing conversion
            }
        } else {
            // FacEC[i] �Ѿ��� float ���ͣ�����ת��
            if (max_val < FacEC[i]) {
                max_val = FacEC[i];
                CriFac = i;
            }
        }
    }

    int FamD1 = wyt_rand(2, 7);
    int CriFacFamD = FamD1 / 2;   // �ؼ�����Ĩȥһ��
    CriFacFamD = std::min(CriFacFamD, (int) FacFamSeq[CriFac].size()); // ���ⳬ����Χ

    // 1. **�ӹؼ����� CriFac ɾ����**
    for (int i = 0; i < CriFacFamD; i++) {
        int Pos = rand() % FacFamSeq[CriFac].size();
        int FamExt = FacFamSeq[CriFac][Pos];
        FacFamSeq[CriFac].erase(FacFamSeq[CriFac].begin() + Pos);
        FamsExtracted.push_back(FamExt);
        UpdateFamErasedFromFac(Pos, CriFac);
    }

    // 2. **����������ɾ����**
    vector<int> OtherFacs;
    for (int i = 0; i < FacFamSeq.size(); i++) {
        if (i != CriFac && !FacFamSeq[i].empty()) {
            OtherFacs.push_back(i);
        }
    }

    // ȷ�� OtherFacs ���ǿյ�
    if (!OtherFacs.empty()) {
        int OtherFacFamD = FamD1 - CriFacFamD; // ��������Ĩȥʣ�µ���
        int deletedCount = 0;

        while (deletedCount < OtherFacFamD && !OtherFacs.empty()) {
            // �������ѡ��һ������
            int RandomFacIndex = rand() % OtherFacs.size();
            int SelectedFac = OtherFacs[RandomFacIndex];

            // ���㵱ǰ���������ɾ��������
            int extractNum = std::min(OtherFacFamD - deletedCount, (int) FacFamSeq[SelectedFac].size());

            for (int i = 0; i < extractNum; ++i) {
                int Pos = rand() % FacFamSeq[SelectedFac].size();
                int FamExt = FacFamSeq[SelectedFac][Pos];
                FacFamSeq[SelectedFac].erase(FacFamSeq[SelectedFac].begin() + Pos);
                FamsExtracted.push_back(FamExt);
                UpdateFamErasedFromFac(Pos, SelectedFac);
            }

            deletedCount += extractNum;
            // ���ĳ�������Ѿ�û�����ˣ��ʹ� OtherFacs �Ƴ�
            if (FacFamSeq[SelectedFac].empty()) {
                OtherFacs.erase(OtherFacs.begin() + RandomFacIndex);
            }
        }
    }
    // �ӹ�����������ȡ����
    for (int Fam : FamsExtracted) // ����������Ĺ�����ȡһ��
    {
        if (JobSeqInFam[Fam].size() >= 3) {
            //            cout << "���� " << Fam << " ɾ���Ĺ���: ";
            for (int j = 0; j <= JobSeqInFam[Fam].size() / 2; j++) {
                int JobPos = rand() % JobSeqInFam[Fam].size();
                int Job = JobSeqInFam[Fam][JobPos];
                //                cout << Job << " ";  // �����ɾ���Ĺ���
                JobSeqInFam[Fam].erase(JobSeqInFam[Fam].begin() + JobPos);
                JobsExtracted[Fam].push_back(Job);
            }
        }
    }

    GetJFDHierarchy_Forward_New(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy);

}


void
QCCEA::BasedInd_Construction_FamsAndJobs_New_Opeator0(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                      vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                                      vector<int> &FamsExtracted,
                                                      unordered_map<int, vector<int>> &JobsExtracted, int nadirpointMS,
                                                      float nadirpointTEC, int idealpointMS, float idealpointTEC,vector<Individual> &CCEAPopulation) {

    while (!FamsExtracted.empty()) {
        int BestFac = -1;
        int BestPos = -1;
        int Pos = rand() % FamsExtracted.size();
        int CurFam = FamsExtracted[Pos];
        FamsExtracted.erase(FamsExtracted.begin() + Pos);
        float Ind = this->FindBestPosToInsertFamForAllFac_Ind_DR_New(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy, CurFam,
                                                                     BestFac, BestPos, nadirpointMS, nadirpointTEC,
                                                                     idealpointMS, idealpointTEC, CCEAPopulation);
        FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
        RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[BestFac], JobSeqInFam, JFDTime, Hierarchy, BestPos, 0,
                                                 JobSeqInFam[CurFam].size() - 1);
        //cout << "�ڹ���" << BestFac << "��" << BestPos << "λ�ò���, ָ��Ϊ��" << Ind << endl;

        int FamPos = find(begin(FacFamSeq[BestFac]), end(FacFamSeq[BestFac]), CurFam) - begin(FacFamSeq[BestFac]);

        for (int i = 0; i < JobsExtracted[CurFam].size(); i++) {
            int BestJobPos = -1;
            double minFacInd = INT_MAX;
            int CurJob = JobsExtracted[CurFam][i];
            for (int JobPos = 0; JobPos <= JobSeqInFam[CurFam].size(); JobPos++) {
                float FacInd = GetIndForPerFacAfterInsertJob_DR_New(BestFac, FacFamSeq, FacFamSeq[BestFac], JobSeqInFam,
                                                                    JFDTime, Hierarchy, CurFam, FamPos, CurJob,
                                                                    JobPos, nadirpointMS, nadirpointTEC, idealpointMS,
                                                                    idealpointTEC, CCEAPopulation);

                if (FacInd < minFacInd) {
                    minFacInd = FacInd;
                    BestPos = JobPos;
                }

            }
            JobSeqInFam[CurFam].insert(JobSeqInFam[CurFam].begin() + BestPos, JobsExtracted[CurFam][i]);
            RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[BestFac], JobSeqInFam, JFDTime, Hierarchy, FamPos,
                                                     BestPos,BestPos);

        }
    }
}


void
QCCEA::BasedInd_Construction_FamsAndJobs_New_Opeator1(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                      vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                                      vector<int> &FamsExtracted,
                                                      unordered_map<int, vector<int>> &JobsExtracted, int nadirpointMS,
                                                      float nadirpointTEC, int idealpointMS, float idealpointTEC,vector<Individual> &CCEAPopulation) {
    // �������б��Ƴ��Ĺ�����
    for (auto &CurFam: FamsExtracted) {
        int BestFac = -1;
        int BestPos = -1;
        double minFacInd = INT_MAX;
        // �����й������ҵ����ŵĲ���λ��
        for (int Fac = 0; Fac < this->m_Factories; ++Fac) {
            int maxPos = FacFamSeq[Fac].size(); // ���Ŀ��ܲ���λ��
            // �������п��ܵ�λ��
            for (int Pos = 0; Pos <= maxPos; ++Pos) {
                // ���㵱ǰ����λ�õ�ָ��
                float tempFacInd = this->GetIndForPerFacAfterInsertFam_DR_New(Fac, FacFamSeq, FacFamSeq[Fac],
                                                                              JobSeqInFam, JFDTime, Hierarchy, CurFam,
                                                                              Pos, nadirpointMS, nadirpointTEC,
                                                                              idealpointMS, idealpointTEC, CCEAPopulation);
                // ��������λ��
                if (tempFacInd < minFacInd) {
                    minFacInd = tempFacInd;
                    BestPos = Pos;
                    BestFac = Fac;
                    // ��������λ�ã��ֲ�������
                    for (int offset: {-1, 1}) {
                        int neighborPos = Pos + offset;
                        if (neighborPos >= 0 && neighborPos <= maxPos) {
                            float neighborFacInd = this->GetIndForPerFacAfterInsertFam_DR_New(Fac, FacFamSeq,
                                                                                              FacFamSeq[Fac],
                                                                                              JobSeqInFam, JFDTime,
                                                                                              Hierarchy, CurFam,
                                                                                              neighborPos, nadirpointMS,
                                                                                              nadirpointTEC,
                                                                                              idealpointMS,
                                                                                              idealpointTEC, CCEAPopulation);
                            // ��������λ��
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

        // ���빤���嵽����λ��
        FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
        // ���¼ӹ�ʱ��Ͳ㼶��ϵ
        RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[BestFac], JobSeqInFam, JFDTime, Hierarchy, BestPos, 0,
                                                 JobSeqInFam[CurFam].size() - 1);
    }

}


void
QCCEA::BasedInd_Construction_FamsAndJobs_New_Opeator2(vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                      vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                                                      vector<int> &FamsExtracted,
                                                      unordered_map<int, vector<int>> &JobsExtracted, int nadirpointMS,
                                                      float nadirpointTEC, int idealpointMS, float idealpointTEC,vector<Individual> &CCEAPopulation) {

    // ����ѡ������d������뵽��õ�λ�� Insert the extracted Families into the best Positions
    for (int Fam = 0; Fam < FamsExtracted.size(); Fam++) {

        int CurFam = FamsExtracted[Fam];  //��ǰ��
        int BestFac = -1;
        int BestPos = -1;  //��õĹ�����λ��

        if (Fam == FamsExtracted.size() - 1) {
            float Ind = this->FindBestPosToInsertFamForAllFac_Ind_DR_New(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy,
                                                                         CurFam, BestFac, BestPos, nadirpointMS,
                                                                         nadirpointTEC, idealpointMS, idealpointTEC, CCEAPopulation);
        } else {
            float randValue = wyt_rand_include_right(0.0f, 1.0f);

            if (randValue < 0.5) {
                FindBestPosToInsertFamForAllFacs_MakeSpan_New(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy, CurFam,
                                                              BestFac, BestPos);
            } else {
                FindBestPosInsertFamAllFacs_TEC_New(FacFamSeq, JobSeqInFam, JFDTime, Hierarchy, CurFam, BestFac,
                                                    BestPos);
            }
        }
        // ����õĹ��������λ���в��뵱ǰ�� Insert CurFam to bestPos at bestFac
        FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
        RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[BestFac], JobSeqInFam, JFDTime, Hierarchy, BestPos, 0,
                                                 JobSeqInFam[CurFam].size() - 1);

        int FamPos = find(begin(FacFamSeq[BestFac]), end(FacFamSeq[BestFac]), CurFam) - begin(FacFamSeq[BestFac]);
        for (int i = 0; i < JobsExtracted[CurFam].size(); i++) {
            int BestJobPos = -1;
            double minFacInd = INT_MAX;
            int CurJob = JobsExtracted[CurFam][i];
            for (int Pos = 0; Pos <= JobSeqInFam[CurFam].size(); Pos++) {

                float FacInd = GetIndForPerFacAfterInsertJob_DR_New(BestFac, FacFamSeq, FacFamSeq[BestFac], JobSeqInFam,
                                                                    JFDTime,
                                                                    Hierarchy, CurFam, FamPos, CurJob,
                                                                    Pos, nadirpointMS, nadirpointTEC, idealpointMS,
                                                                    idealpointTEC, CCEAPopulation);

                if (FacInd < minFacInd) {
                    minFacInd = FacInd;
                    BestPos = Pos;
                }

            }
            JobSeqInFam[CurFam].insert(JobSeqInFam[CurFam].begin() + BestPos, JobsExtracted[CurFam][i]);
            RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[BestFac], JobSeqInFam, JFDTime, Hierarchy, FamPos,
                                                     BestPos,
                                                     BestPos);
            //cout << "�ڹ���" << BestFac << "��" << CurFam << "���" << BestPos << "λ�ò���, ָ��Ϊ��" << minFacInd << endl;
        }
    }
}


void QCCEA::UpdateArchiveGroupJobSet(int Muu)
{

    vector<Individual> Temp;
    Temp.clear();

    //ռ�Ź�ϵ
    Pareto_relation(m_QCCEAPopulation);

    //����֧��� de-emphasize dominated solutions
    for (int j = 0; j < m_QCCEAPopulation.size(); j++)
    {
        for (auto & i : m_QCCEAPopulation)
        {
            if (i.flag == 0)
            {
                if (i.pareto_rel[j] == 1)
                {
                    i.flag = 999;
                }
            }
        }
    }


    for (auto & i : m_QCCEAPopulation)
    {
        if (i.flag == 0)
            Temp.push_back(i);
    }

    //������mu(Archive��Ⱥ����)
    int mu = Muu;
    if (Temp.size() < mu * m_RefSize)
    {
        while (true)
        {
            for (auto & i : m_QCCEAPopulation)
            {
                if (i.flag == 999)
                {
                    Temp.push_back(i);
                    if (Temp.size() == 10 * m_RefSize)
                        break;
                }
            }
            break;
        }
    }

    m_QCCEAPopulation.clear();
    for (auto & i : Temp)
        m_QCCEAPopulation.push_back(i);


    //��һ��
    Normalize(m_QCCEAPopulation, m_NadirPointMS, m_NadirPointTEC, m_IdealPointMS, m_IdealPointTEC);

    //��������ָ��
    calc_convergence_ind(m_QCCEAPopulation);

    //�����ɢָ��
    calc_distribution_ind(m_QCCEAPopulation);

    int n, n1, n2, n3, nrank;
    int the_one;
    n = 0;   //Qs.size  Ԥѡ��ĸ��弯��
    n1 = 0;  //Q.size   �ο����ĸ�������
    n2 = 0;  //Qth.size ��distribution threshold�����ļ���
    n3 = 0;  //Qd.size  ��֧�����ļ���

    vector<Individual> TempQCCEAPopulation;
    TempQCCEAPopulation.clear();

    //ռ�Ź�ϵ
    Pareto_relation(m_QCCEAPopulation);

    n1 = m_QCCEAPopulation.size();

    while (n1 > 0)
    {
        for (int i = 0; i < m_QCCEAPopulation.size(); i++)
        {
            if (m_QCCEAPopulation[i].flag == 0)
            {
                the_one = i;
                break;
            }
        }

        m_QCCEAPopulation[the_one].flag = 1;
        TempQCCEAPopulation.push_back(m_QCCEAPopulation[the_one]);
        n++;


        //���������� de-emphasize neighbors
        for (auto & i : m_QCCEAPopulation)
        {
            if (i.flag == 0)
            {
                if (i.distribution_ind[the_one] < thr_zeta)
                {
                    i.flag = 999;
                    n2++;
                }
            }
        }

        //����֧��� de-emphasize dominated solutions
        for (auto & i : m_QCCEAPopulation)
        {
            if (i.flag == 0)
            {
                if (i.pareto_rel[the_one] == 1)
                {
                    i.flag = 999;
                    n3++;
                }
            }
        }

        //ʣ������ number of the rest
        n1 = 0;
        for (auto & i : m_QCCEAPopulation)
        {
            if (i.flag == 0)
            {
                n1++;
            }
        }
    }

    //���� threshold
    float oldthr_zeta = thr_zeta;
    float ratio = static_cast<float>(n) / static_cast<float>(m_RefSize);
    if (n3 < (mu - 1) * m_RefSize)
        thr_zeta = thr_zeta * exp((ratio - 1.0f) / (2 * 1.0f));
    else
        thr_zeta = oldthr_zeta;


    while (n < m_RefSize)
    {
        //ʣ������ number of the rest
        n1 = 0;
        for (auto & i : m_QCCEAPopulation)
        {
            if (i.flag == 0)
            {
                n1++;
            }
        }

        if (n1 == 0)
        {
            for (auto & i : m_QCCEAPopulation)
            {
                if (i.flag == 999)
                {
                    i.flag = 0;
                }
            }
        }

        //ѡ��һ��
        for (int i = 0; i < m_QCCEAPopulation.size(); i++)
        {
            if (m_QCCEAPopulation[i].flag == 0)
            {
                the_one = i;
                break;
            }
        }
        m_QCCEAPopulation[the_one].flag = 1;
        TempQCCEAPopulation.push_back(m_QCCEAPopulation[the_one]);
        n++;

        //���������� de-emphasize neighbors
        for (auto & i : m_QCCEAPopulation)
        {
            if (i.flag == 0)
            {
                if (i.distribution_ind[the_one] < oldthr_zeta)
                {
                    i.flag = 999;
                    n2++;
                }
            }
        }

        //����֧��� de-emphasize dominated solutions
        for (auto & i : m_QCCEAPopulation)
        {
            if (i.flag == 0)
            {
                if (i.pareto_rel[the_one] == 1)
                {
                    i.flag = 999;
                    n3++;
                }
            }
        }
    }

    vector<Individual> Temp1;
    Temp1.clear();
    // m_RefSize
    for (auto & i : TempQCCEAPopulation)
    {
        Temp1.push_back(i);
    }

    //����
    Speed_mutation(TempQCCEAPopulation, Temp1, m_QCCEAPopulation);

    for (int PS = 0; PS < m_RefSize; PS++)
    {
        // ���µ���������ֵ��
        m_ArchiveJobSeqinFamArray[PS] = TempQCCEAPopulation[PS].m_JobSeqInFamArray;  //����������
        m_ArchiveFacFamSeqArray[PS] = TempQCCEAPopulation[PS].m_FacFamSeqArray;   //����������
        m_ArchiveSpanArray[PS] = TempQCCEAPopulation[PS].MS;  // ����깤ʱ��
        m_ArchiveTECArray[PS] = TempQCCEAPopulation[PS].TEC;  //���ܺ�
        m_ArchiveSpeedVector[PS] = TempQCCEAPopulation[PS].m_SpeedVector;

        //��������Ⱥ
        m_FacFamSeqArray[PS] = m_ArchiveFacFamSeqArray[PS];

        //���¹�����Ⱥ
        m_JobSeqinFamArray[PS] = m_ArchiveJobSeqinFamArray[PS];

        m_SpeedMatrix = m_ArchiveSpeedVector[PS];
        //�õ������Ĵ���ʱ�� ����λʱ��ӹ��ܺ�
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }
        GetJFDHierarchy_Forward_New(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], m_TempJFDTimePop[PS], m_TempHierarchyPop[PS]);
    }
}
void QCCEA::SaveTECandDeMS(vector<Individual> &CCEAPopulation) {
    vector<Individual> tempQCCEAPopulation;
    tempQCCEAPopulation.clear();
    vector<Individual> tempQCCEAPopulation2;
    tempQCCEAPopulation2.clear();
    vector<Individual> tempQCCEAPopulation3;
    tempQCCEAPopulation3.clear();
    vector<Individual> tempQCCEAPopulation4;
    tempQCCEAPopulation4.clear();

    int orgSize = CCEAPopulation.size();
    for (int PS = 0; PS < orgSize; PS++) {
        tempQCCEAPopulation.push_back(CCEAPopulation[PS]);
        bool flag2 = false;

        vector<int> FacSpan(m_Factories);//�飬�����깤ʱ��

        //���ܲ���1
        vector<vector<int>> DelayTime;
        DelayTime.clear();
        DelayTime.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
        GetDelayTime_Forward(tempQCCEAPopulation[PS].m_FacFamSeqArray, tempQCCEAPopulation[PS].m_JobSeqInFamArray,
                             m_TempJFDTime, FacSpan, DelayTime);
        //int makespannn = GetJFDTime_Forward(tempQCCEAPopulation[PS].m_FacFamSeqArray, tempQCCEAPopulation[PS].m_JobSeqInFamArray, m_TempJFDTime, FacSpan);
        bool sign = false;

        for (int j = 0; j < m_Jobs; j++) {
            for (int i = 1; i < m_Machines; i++) {
                if ((DelayTime[j][i] > 0) && (i < m_Machines - 1)) {

                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);

                    m_TureJobOpertime[j][i + 1] = static_cast<int>(m_JobOperPTime[j][i + 1] /m_Speed[tempQCCEAPopulation[PS].m_SpeedVector[j][
                                                                                  i + 1]]);

                    if ((m_TempJFDTime[j][i] < (m_TempJFDTime[j][i + 1] - m_TureJobOpertime[j][i + 1]))) {
                        int Speedlevel = tempQCCEAPopulation[PS].m_SpeedVector[j][i] - 1;

                        for (int level = Speedlevel; level > 0; level--) {
                            if ((static_cast<int>(m_JobOperPTime[j][i] / m_Speed[level]) -
                                 m_TureJobOpertime[j][i]) < DelayTime[j][i]) {
                                if ((m_TempJFDTime[j][i] +
                                     (static_cast<int>(m_JobOperPTime[j][i] / m_Speed[level]) -
                                      m_TureJobOpertime[j][i])) <
                                    (m_TempJFDTime[j][i + 1] - m_TureJobOpertime[j][i + 1])) {
                                    tempQCCEAPopulation[PS].m_SpeedVector[j][i] = level;
                                    sign = true;
                                } else
                                    break;
                            } else
                                break;
                        }
                    }

                } else if ((DelayTime[j][i] > 0) && (i == m_Machines - 1)) {
                    int Speedlevel = tempQCCEAPopulation[PS].m_SpeedVector[j][i] - 1;
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] /
                                                                      m_Speed[tempQCCEAPopulation[PS].m_SpeedVector[j][i]]);
                    for (int level = Speedlevel; level > 0; level--) {
                        if ((static_cast<int>(m_JobOperPTime[j][i] / m_Speed[level]) -
                             m_TureJobOpertime[j][i]) < DelayTime[j][i]) {
                            tempQCCEAPopulation[PS].m_SpeedVector[j][i] = level;
                            sign = true;
                        } else
                            break;
                    }
                }
            }//end machine
        }//end job

        if (sign) {
            int Makespan1 = 0;
            float TotalEC1 = 0.0;

            m_SpeedMatrix = tempQCCEAPopulation[PS].m_SpeedVector;

            //�õ������Ĵ���ʱ�� ����λ�ӹ�ʱ���ܺ�
            for (int j = 0; j < m_Jobs; j++) {
                for (int i = 0; i < m_Machines; i++) {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }


            Makespan1 = GetJFDTime_Forward(tempQCCEAPopulation[PS].m_FacFamSeqArray,
                                           tempQCCEAPopulation[PS].m_JobSeqInFamArray, m_TempJFDTime, FacSpan);

            TotalEC1 = GetTECForAllFacByJFD(tempQCCEAPopulation[PS].m_FacFamSeqArray,
                                            tempQCCEAPopulation[PS].m_JobSeqInFamArray, m_TempJFDTime);

            //���
            //CheckSol(tempQCCEAPopulation[PS].m_FacFamSeqArray, tempQCCEAPopulation[PS].m_JobSeqInFamArray, Makespan1);
            //CheckSolTEC(tempQCCEAPopulation[PS].m_FacFamSeqArray, tempQCCEAPopulation[PS].m_JobSeqInFamArray, TotalEC1);
            tempQCCEAPopulation[PS].MS = Makespan1;
            tempQCCEAPopulation[PS].TEC = TotalEC1;

            //�ж�
            if ((CCEAPopulation[PS].MS >= tempQCCEAPopulation[PS].MS) &&
                (CCEAPopulation[PS].TEC > tempQCCEAPopulation[PS].TEC)) {
                CCEAPopulation.push_back(tempQCCEAPopulation[PS]);
                flag2 = true;
            }
        }


        //���ܲ���3����1�Ļ����ϣ��ڲ�����Makespan������£����ͷǹؼ�������v��������TEC
        vector<int> FacSpan3(m_Factories);
        if (flag2) {
            //float gapTEC = CCEAPopulation[PS].TEC - tempQCCEAPopulation[PS].TEC;
            tempQCCEAPopulation3.push_back(tempQCCEAPopulation[PS]);
            vector<int> notcirfac;
            notcirfac.clear();
            for (int fac = 0; fac < m_Factories; fac++) {
                if (tempQCCEAPopulation3.back().MS > FacSpan[fac]) {
                    notcirfac.push_back(fac);
                }
                FacSpan3[fac] = FacSpan[fac];
            }
            int diedai = 0;
            float tempTEC = 0.0;

            do {
                diedai++;
                if (diedai == 1)
                    tempTEC = tempQCCEAPopulation[PS].TEC;

                if (diedai == 2)
                    tempTEC = tempQCCEAPopulation3.back().TEC;
                if (diedai == 3)
                    break;

                bool judgeV = true;
                vector<vector<int>> tempSpeedMatrix;
                tempSpeedMatrix.clear();
                tempSpeedMatrix.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
                tempSpeedMatrix = tempQCCEAPopulation3.back().m_SpeedVector;

                //�Էǹؼ��������б���
                for (int f : notcirfac) {
                    judgeV = true;

                    for (int fam = 0; fam < tempQCCEAPopulation3.back().m_FacFamSeqArray[f].size(); fam++) {
                        int CurFam = tempQCCEAPopulation3.back().m_FacFamSeqArray[f][fam];
                        for (int job = 0; job < tempQCCEAPopulation3.back().m_JobSeqInFamArray[CurFam].size(); job++) {
                            int CurJob = tempQCCEAPopulation3.back().m_JobSeqInFamArray[CurFam][job];
                            for (int m = 0; m < m_Machines; m++) {
                                if (tempQCCEAPopulation3.back().m_SpeedVector[CurJob][m] > 0) {
                                    tempQCCEAPopulation3.back().m_SpeedVector[CurJob][m] =
                                            tempQCCEAPopulation3.back().m_SpeedVector[CurJob][m] - 1;
                                    judgeV = false;
                                }

                            }
                        }
                    }
                    //�����ٶ��б仯
                    if (!judgeV) {
                        //faccount++;
                        m_SpeedMatrix = tempQCCEAPopulation3.back().m_SpeedVector;

                        //�õ������Ĵ���ʱ�� �� ��λ�ӹ��ܺ�PECϵ��
                        for (int j = 0; j < m_Jobs; j++) {
                            for (int i = 0; i < m_Machines; i++) {

                                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);

                                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                            }
                        }

                        //�õ�����f���ٺ��Span
                        int tempnotcirSpan = GetJFDTime_Forward_InFactory(
                                tempQCCEAPopulation3.back().m_FacFamSeqArray[f],
                                tempQCCEAPopulation3.back().m_JobSeqInFamArray, m_TempJFDTime);//�õ�ÿ���������깤ʱ��
                        //����ı��ٶȺ�δ����Makspan������������
                        if (tempnotcirSpan > tempQCCEAPopulation3.back().MS) {
                            tempQCCEAPopulation3.back().m_SpeedVector = tempSpeedMatrix;
                            //countfac++;
                            //judgeV = true;
                        } else {
                            tempSpeedMatrix = tempQCCEAPopulation3.back().m_SpeedVector;
                        }

                    }

                }

                /*if (judgeV && (countfac == m_Factories))
                    break;*/

                int Makespan3 = 0;
                float TotalEC3 = 0.0;

                m_SpeedMatrix = tempQCCEAPopulation3.back().m_SpeedVector;

                //�õ������Ĵ���ʱ�� �� ��λ�ӹ��ܺ�PECϵ��
                for (int j = 0; j < m_Jobs; j++) {
                    for (int i = 0; i < m_Machines; i++) {
                        m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] /
                                                                          m_Speed[m_SpeedMatrix[j][i]]);
                        UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                    }
                }

                Makespan3 = GetJFDTime_Forward(tempQCCEAPopulation3.back().m_FacFamSeqArray,
                                               tempQCCEAPopulation3.back().m_JobSeqInFamArray, m_TempJFDTime, FacSpan3);

                TotalEC3 = GetTECForAllFacByJFD(tempQCCEAPopulation3.back().m_FacFamSeqArray,
                                                tempQCCEAPopulation3.back().m_JobSeqInFamArray, m_TempJFDTime);

                //���
                //CheckSol(tempQCCEAPopulation3.back().m_FacFamSeqArray, tempQCCEAPopulation3.back().m_JobSeqInFamArray, Makespan3);
                //CheckSolTEC(tempQCCEAPopulation3.back().m_FacFamSeqArray, tempQCCEAPopulation3.back().m_JobSeqInFamArray, TotalEC3);
                tempQCCEAPopulation3.back().MS = Makespan3;
                tempQCCEAPopulation3.back().TEC = TotalEC3;

                //�ж�
                if ((CCEAPopulation[PS].MS >= tempQCCEAPopulation3.back().MS) &&
                    (tempQCCEAPopulation[PS].TEC > tempQCCEAPopulation3.back().TEC)) {
                    CCEAPopulation.push_back(tempQCCEAPopulation3.back());
                }

            } while (tempTEC > tempQCCEAPopulation3.back().TEC);

        }

        //���ܲ���4����3�Ļ����ϣ��ڲ�����ԭTEC������£���߹ؼ�������v��������MS
        vector<int> FacSpan4(m_Factories);
        if (flag2) {
            //float gapTEC = CCEAPopulation[PS].TEC - tempQCCEAPopulation[PS].TEC;
            tempQCCEAPopulation4.push_back(tempQCCEAPopulation3.back());
            vector<int> cirfac;
            cirfac.clear();
            for (int fac = 0; fac < m_Factories; fac++) {
                if (tempQCCEAPopulation4.back().MS == FacSpan[fac]) {
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
                    tempMS = tempQCCEAPopulation4.back().MS;
                    //cirfac.clear();
                    vector<int> tempcirfac;
                    tempcirfac.clear();
                    //�жϹؼ������Ƿ�ı�
                    for (int fac = 0; fac < m_Factories; fac++) {
                        if (tempQCCEAPopulation4.back().MS == FacSpan4[fac]) {
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
                tempSpeedMatrix = tempQCCEAPopulation4.back().m_SpeedVector;

                //int countfac = 0;
                //int faccount = 0;
                //�Էǹؼ��������б���
                for (int f : cirfac) {
                    judgeV = true;
                    //int speedLevel;
                    for (int fam = 0; fam < tempQCCEAPopulation4.back().m_FacFamSeqArray[f].size(); fam++) {
                        int CurFam = tempQCCEAPopulation4.back().m_FacFamSeqArray[f][fam];
                        for (int job = 0; job < tempQCCEAPopulation4.back().m_JobSeqInFamArray[CurFam].size(); job++) {
                            int CurJob = tempQCCEAPopulation4.back().m_JobSeqInFamArray[CurFam][job];
                            for (int m = 0; m < m_Machines; m++) {
                                if (tempQCCEAPopulation4.back().m_SpeedVector[CurJob][m] < 2) {
                                    tempQCCEAPopulation4.back().m_SpeedVector[CurJob][m] =
                                            tempQCCEAPopulation4.back().m_SpeedVector[CurJob][m] + 1;
                                    judgeV = false;
                                }

                            }
                        }
                    }
                    //�����ٶ��б仯
                    if (!judgeV) {
                        //faccount++;
                        m_SpeedMatrix = tempQCCEAPopulation4.back().m_SpeedVector;

                        //�õ������Ĵ���ʱ�� �� ��λ�ӹ��ܺ�PECϵ��
                        for (int j = 0; j < m_Jobs; j++) {
                            for (int i = 0; i < m_Machines; i++) {
                                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] /
                                                                                  m_Speed[m_SpeedMatrix[j][i]]);
                                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                            }
                        }

                        //�õ�����f���ٺ��TEC
                        GetJFDTime_Forward(tempQCCEAPopulation4.back().m_FacFamSeqArray,
                                           tempQCCEAPopulation4.back().m_JobSeqInFamArray, m_TempJFDTime, FacSpan4);
                        float tempTotalEC = GetTECForAllFacByJFD(tempQCCEAPopulation4.back().m_FacFamSeqArray,
                                                                 tempQCCEAPopulation4.back().m_JobSeqInFamArray,
                                                                 m_TempJFDTime);


                        //����ı��ٶȺ�δ����ԭTEC������������
                        if (tempTotalEC > CCEAPopulation[PS].TEC) {
                            tempQCCEAPopulation4.back().m_SpeedVector = tempSpeedMatrix;
                            break;
                        }

                    }

                }

                int Makespan4 = 0;
                float TotalEC4 = 0.0;

                m_SpeedMatrix = tempQCCEAPopulation4.back().m_SpeedVector;

                //�õ������Ĵ���ʱ�� �� ��λ�ӹ��ܺ�PECϵ��
                for (int j = 0; j < m_Jobs; j++) {
                    for (int i = 0; i < m_Machines; i++) {
                        m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] /
                                                                          m_Speed[m_SpeedMatrix[j][i]]);
                        UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                    }
                }

                Makespan4 = GetJFDTime_Forward(tempQCCEAPopulation4.back().m_FacFamSeqArray,
                                               tempQCCEAPopulation4.back().m_JobSeqInFamArray, m_TempJFDTime, FacSpan4);

                TotalEC4 = GetTECForAllFacByJFD(tempQCCEAPopulation4.back().m_FacFamSeqArray,
                                                tempQCCEAPopulation4.back().m_JobSeqInFamArray, m_TempJFDTime);

                //���
                //CheckSol(tempQCCEAPopulation4.back().m_FacFamSeqArray, tempQCCEAPopulation4.back().m_JobSeqInFamArray, Makespan4);
                //CheckSolTEC(tempQCCEAPopulation4.back().m_FacFamSeqArray, tempQCCEAPopulation4.back().m_JobSeqInFamArray, TotalEC4);
                tempQCCEAPopulation4.back().MS = Makespan4;
                tempQCCEAPopulation4.back().TEC = TotalEC4;

                //�ж�
                if ((CCEAPopulation[PS].MS > tempQCCEAPopulation4.back().MS) &&
                    (CCEAPopulation[PS].TEC >= tempQCCEAPopulation4.back().TEC)) {
                    CCEAPopulation.push_back(tempQCCEAPopulation4.back());
                }

            } while (tempMS > tempQCCEAPopulation4.back().MS);

        }
    }

}

//Q learning

int QCCEA::DetermineAction(int state, const vector<vector<double>> &Q) const {

    double rande = wyt_rand_include_right(0.0f, 1.0f);
    if (rande < e) {
        vector<int> QI = findMaxIndices(Q[state]);
        return int(QI[wyt_rand(QI.size())]);
    } else {
        return wyt_rand(Q[state].size());
    }
}

//int QCCEA::DetermineState(int PreMS, int actMS, float PreTEC, float actTEC) {
//    if (actMS < PreMS && actTEC < PreTEC) {
//        return 0;//��
//    } else if (actTEC < PreTEC || actMS < PreMS) {
//        return 1;//�Ľ�һ���� �Ϻ�
//    } else {
//        return 2;//û�Ľ�
//    }
//}


// �ȽϺ��������ڰ�makespan����
bool QCCEA::CompareByMakespan(const Individual &a, const Individual &b) {
    return a.MS < b.MS;
}

// �ȽϺ��������ڰ�TEC����
bool QCCEA::CompareByTEC(const Individual &a, const Individual &b) {
    return a.TEC < b.TEC;
}

// ȷ�������״̬
int QCCEA::DetermineState(int ObjectMS, float ObjectTEC, const std::vector<Individual> &QCCEAPopulation) {

    int N = QCCEAPopulation.size();

    // ��makespan����
    std::vector<Individual> SortedByMakespan = QCCEAPopulation;
    std::sort(SortedByMakespan.begin(), SortedByMakespan.end(), CompareByMakespan);
    int MedianMakeSpan = SortedByMakespan[N / 2].MS;

    // ��TEC����
    std::vector<Individual> SortedByTEC = QCCEAPopulation;
    std::sort(SortedByTEC.begin(), SortedByTEC.end(), CompareByTEC);
    float MedianTEC = SortedByTEC[N / 2].TEC;

    // ȷ��״̬
//    if (ObjectMS <= MedianMakeSpan && ObjectTEC <= MedianTEC) {
//        return 0; // ״̬0��MS��TEC���Ϻ�
//    } else if (ObjectMS <= MedianMakeSpan && ObjectTEC > MedianTEC) {
//        return 1; // ״̬1��MS�Ϻã�TEC�ϲ�
//    } else if (ObjectMS > MedianMakeSpan && ObjectTEC <= MedianTEC) {
//        return 1; // ״̬2��MS�ϲTEC�Ϻ�
//    } else {
//        return 2; // ״̬3��MS��TEC���ϲ�
//    }

    if (ObjectMS <= MedianMakeSpan && ObjectTEC <= MedianTEC) {
        return 0; // ״̬0��MS��TEC���Ϻ�
    } else if (ObjectMS <= MedianMakeSpan || ObjectTEC <= MedianTEC) {
        return 1; // ״̬1��MS�Ϻã�TEC�ϲ� �� MS�ϲTEC�Ϻ�
    } else {
        return 2; // ״̬2��MS��TEC���ϲ�
    }

}

void QCCEA::UpdateQ(int oldstate, int newstate, int act, double reward, vector<vector<double>> &Q) const {
    //double e = 0.75;̰������  discountFactor = 0.2;�ۿ����� rewardFactor = 0.75;ѧϰ��
    if (newstate == 0 || newstate == 1) {
        Q[oldstate][act] = (1.0 - rewardFactor) * Q[oldstate][act] + rewardFactor * (reward + discountFactor **max_element(Q[newstate].begin(),Q[newstate].end()));
    }
}

double QCCEA::CaculateReward(int oldMS, int newMS, float oldTEC, float newTEC,
                             const std::vector<Individual> &population) {

    // ����Ŀ�꺯�������ֵ����Сֵ
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

    // ��ֹ�������
    if (MS_max == MS_min || TEC_max == TEC_min) {
        return 0.0;  // ���߸����������һ�������Ĭ�Ͻ���ֵ
    }

    // ���㽱��ֵ
    double reward = 0.0;
    if (newMS <= oldMS && newTEC <= oldTEC) {   // �����½�;ɽ���ÿ��Ŀ���ϵ�֧���ϵ
        reward = 2.0 + ((MS_max - newMS) / static_cast<double>(MS_max - MS_min)) +
                 ((TEC_max - newTEC) / static_cast<double>(TEC_max - TEC_min));
    } else if (oldMS <= newMS && oldTEC <= newTEC) {
        reward = 0.0;
    } else {
        reward = 1.0 + ((MS_max - newMS) / static_cast<double>(MS_max - MS_min)) +
                 ((TEC_max - newTEC) / static_cast<double>(TEC_max - TEC_min));
    }

    return reward;
}


void QCCEA::PerformAction_Fams(int act, vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqinFam,
                               vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy,
                               vector<Individual> &CCEAPopulation, vector<int> &FacSpan, vector<float> &FacEC,
                               int &ObjectMS, float &ObjectTEC) {
    switch (act) {
        case 0:
            cout<<"0"<<endl;
            //��makespan���Ĺ������������й�������������õ�λ�ã���������ָ����룩
            BasedInd_RandFamInFacTobestPos_New(FacFamSeq, JobSeqinFam, JFDTime, Hierarchy, m_NadirPointMS,
                                               m_NadirPointTEC, m_IdealPointMS, m_IdealPointTEC,
                                               CCEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
            cout<<"0 end "<<endl;
            break;
        case 1:
            cout<<"1"<<endl;
           BasedInd_SwapFam_New(FacFamSeq, JobSeqinFam, JFDTime, Hierarchy, m_NadirPointMS, m_NadirPointTEC,
                                 m_IdealPointMS, m_IdealPointTEC, CCEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
           cout<<"1 end "<<endl;
            break;
        case 2:
//            cout<<"2"<<endl;
//            BasedInd_DeAndConFams_New(FacFamSeq, JobSeqinFam, JFDTime, Hierarchy, m_NadirPointMS, m_NadirPointTEC,
//                                      m_IdealPointMS, m_IdealPointTEC,
//                                      CCEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
//           cout<<"2 end "<<endl;
            break;
        case 3:
//           cout<<"3  "<<endl;
//            BasedIndLS_SetupSwap_New(FacFamSeq, JobSeqinFam, JFDTime, Hierarchy, m_NadirPointMS, m_NadirPointTEC,
//                                     m_IdealPointMS, m_IdealPointTEC, CCEAPopulation, FacSpan, FacEC, ObjectMS,
//                                     ObjectTEC);
//            cout<<"3 end "<<endl;
            break;
        case 4:
            cout<<"4  "<<endl;
            BasedInd_ShiftFam_New(FacFamSeq, JobSeqinFam, JFDTime, Hierarchy, m_NadirPointMS, m_NadirPointTEC,
                                  m_IdealPointMS, m_IdealPointTEC,
                                  CCEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
            cout<<"4 end   "<<endl;
            break;

        default:
            std::cerr << "Error: Invalid action " << act << " in PerformAction_Fams." << std::endl;
    }
    // ���¼����ܵ� MS �� TEC ȷ�����½�
    GetMSandTECForPerandToalFacByJFD(FacFamSeq, JobSeqinFam, JFDTime, FacSpan, FacEC, ObjectMS, ObjectTEC);
}

vector<int> QCCEA::findMaxIndices(const vector<double> &vec) {
    vector<int> maxIndices;
    if (vec.empty()) {
        return maxIndices; // ���ؿյ���������
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

void QCCEA::InitialPop(vector<vector<vector<int>>> &JFDTimePop, vector<vector<vector<int>>> &HierarchyPop) {

    //��ʼ�������� Initialize are set and Set flags
    for (int PS = 0; PS < m_RefSize; PS++) {

        //�õ������Ĵ���ʱ�� ����λʱ��ӹ��ܺ�
        for (int j = 0; j < m_Jobs; j++) {
            for (int i = 0; i < m_Machines; i++) {
                //m_SpeedMatrix[j][i] = rand() % 3;
                m_SpeedMatrix[j][i] = 0;
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }

        vector<int> FamPrmu(m_Families), FacSpan(m_Factories);//�飬�����깤ʱ��
        vector<float> FacEC(m_Factories);
        vector<vector<int>> JobSeqInFam, FacFamSeq; //������������У����ڹ��������У������ڹ���������
        if (PS == 0 || PS == 1) //��һ���ο�����LPT������������� Generate job sequence in each family and family sequence using LPT
        {
            this->SortJobsInFam(0, JobSeqInFam); //LPT
            this->CombineSortedSequences(2, 3, FamPrmu);

        } else if ( PS == 2|| PS == 3) {
            this->SortJobsInFam(1, JobSeqInFam); //SPT
            this->CombineSortedSequences(4, 5, FamPrmu);
        } else // ������ɹ������к������� Generate job sequence in each family and family sequence randomly
        {
            JobSeqInFam = this->m_JobsInEachFamily;
            for (auto & fam : JobSeqInFam)
                shuffle(fam.begin(), fam.end(), std::mt19937(std::random_device()()));//�������ڹ���˳��
            for (int fam = 0; fam < FamPrmu.size(); fam++)
                FamPrmu[fam] = fam;
            shuffle(FamPrmu.begin(), FamPrmu.end(), std::mt19937(std::random_device()()));//������˳��
            // ������ɵ�������
        }

        FacFamSeq.clear();
        FacFamSeq.resize(this->m_Factories);

        int CurFam = -1;
        int BestFac = -1, BestPos = -1;

        int Fam = 0;
        for (auto & Fac : FacFamSeq) {
            CurFam = FamPrmu[Fam];
            Fac.insert(Fac.begin() + 0, CurFam);
            // ˢ�²����빤��-����Ϣ
            RefreshJFDTimeHierarchy_InFactory_Insert(Fac, JobSeqInFam, JFDTimePop[PS], HierarchyPop[PS], 0,
                                                     0, JobSeqInFam[CurFam].size() - 1);
            Fam++;
        }
        if (PS == 0 || PS == 1) {
            for (; Fam < FamPrmu.size(); Fam++) {
                CurFam = FamPrmu[Fam];
                BestFac = -1;
                BestPos = -1;
                this->FindBestPosToInsertFamForAllFacs_MakeSpan_New(FacFamSeq, JobSeqInFam, JFDTimePop[PS], HierarchyPop[PS],
                                                              CurFam, BestFac, BestPos);
                FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
                // ˢ�²����빤��-����Ϣ
                RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[BestFac], JobSeqInFam, JFDTimePop[PS],
                                                         HierarchyPop[PS],
                                                         BestPos, 0, JobSeqInFam[CurFam].size() - 1);
            }
            // ���ÿ�������������У���������Ӧ����Ⱥ
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
                this->FindBestPosInsertFamAllFacs_TEC_New(FacFamSeq, JobSeqInFam, JFDTimePop[PS], HierarchyPop[PS],
                                                          CurFam, BestFac, BestPos);
                FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
                // ˢ�²����빤��-����Ϣ
                RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[BestFac], JobSeqInFam, JFDTimePop[PS],
                                                         HierarchyPop[PS],
                                                         BestPos, 0, JobSeqInFam[CurFam].size() - 1);
            }
        } else {
                for (; Fam < FamPrmu.size(); Fam++) {
                    CurFam = FamPrmu[Fam];
                    BestFac = -1;
                    BestPos = -1;
                    this->FindBestPosToInsertFamForAllFacs_MC_New(FacFamSeq, JobSeqInFam, JFDTimePop[PS],
                                                                  HierarchyPop[PS],
                                                                  CurFam, BestFac, BestPos);
                    FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
                    // ˢ�²����빤��-����Ϣ
                    RefreshJFDTimeHierarchy_InFactory_Insert(FacFamSeq[BestFac], JobSeqInFam, JFDTimePop[PS],
                                                             HierarchyPop[PS],
                                                             BestPos, 0, JobSeqInFam[CurFam].size() - 1);
                }
        }

        int MS = 0;
        float TEC = 0;
        // cout << endl << "��Ⱥ�е�" << PS << "������" << endl;
        GetMSandTECForPerandToalFacByJFD(FacFamSeq, JobSeqInFam, JFDTimePop[PS], FacSpan, FacEC, MS, TEC);
        // cout << "Makespan��" << MS << "\t" << "TEC��" << TEC << endl;

        // ��ʼ������������ֵ��
        m_ArchiveJobSeqinFamArray[PS] = JobSeqInFam;  //����������
        m_ArchiveFacFamSeqArray[PS] = FacFamSeq;   //����������
        m_ArchiveFacSpanArray[PS] = FacSpan;  //���й������깤ʱ��
        m_ArchiveFacECArray[PS] = FacEC;  //���й������ܺ�
        m_ArchiveSpanArray[PS] = MS;  // ����깤ʱ��
        m_ArchiveTECArray[PS] = TEC;  //���ܺ�
        m_ArchiveSpeedVector[PS] = m_SpeedMatrix;

        //���makespan��TEC
//        this->CheckSol(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], m_ArchiveSpanArray[PS]);
//        this->CheckSolTEC(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], m_ArchiveTECArray[PS]);

        //���²ο���
        m_QCCEAPopulation[PS].m_FacFamSeqArray = FacFamSeq;
        m_QCCEAPopulation[PS].m_JobSeqInFamArray = JobSeqInFam;
        m_QCCEAPopulation[PS].MS = MS;
        m_QCCEAPopulation[PS].TEC = TEC;
        m_QCCEAPopulation[PS].m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���
    }

    // ��ʼ������Ⱥ Initilize Familiy-sequence population, i.e., PS1
    for (int PS = 0; PS < m_PS1; PS++) {
        //���������Ľ⸳��������
        m_FacFamSeqArray[PS] = m_ArchiveFacFamSeqArray[PS];
        m_SpanArray1[PS] = m_ArchiveSpanArray[PS];
        m_TECArray1[PS] = m_ArchiveTECArray[PS];
        m_Map1[PS] = PS;  //��ǵ������ж�Ӧ��jobpop���
    }

    // ��ʼ��������Ⱥ Initialize Job-Sequence population, i.e., PS2
    for (int PS = 0; PS < m_PS2; PS++) {
        //ǰAS����
        m_JobSeqinFamArray[PS] = m_ArchiveJobSeqinFamArray[PS];
        m_SpanArray2[PS] = m_ArchiveSpanArray[PS];
        m_TECArray2[PS] = m_ArchiveTECArray[PS];
        m_Map2[PS] = PS;  //��Ƕ�Ӧ�ĵ�������fampop���
    }

    //��һ��
    Normalize(m_QCCEAPopulation, m_NadirPointMS, m_NadirPointTEC, m_IdealPointMS, m_IdealPointTEC);
    //cout << "nadirpointMS��" << m_NadirPointMS << "\tnadirpointTEC��" << m_NadirPointTEC << endl;

//    //��������ָ��
    calc_convergence_ind(m_QCCEAPopulation);
//
//    //�����ɢָ��
    calc_distribution_ind(m_QCCEAPopulation);
//      cout<<"InitialPop end"<<endl;

}

void QCCEA::EvolutionProcess(int mu) {

    vector<vector<int>> FacFamSeq(this->m_Factories);
    vector<vector<int>> JobSeqInFam;
    m_IterNum = 0;
    thr_zeta = 1.0;
    long InitTime = ::GetElapsedProcessTime();

    vector<vector<vector<double>>> QFam(m_PS1, vector<vector<double>>(3, vector<double>(5)));//q��
    // �����ƻ����Ӻ��ع����ӵ�Ȩ��
    vector<vector<double>> scores_combined(m_PS1, vector<double>(3, 1.0));
    // �ƻ����ع����ӵĳ�ʼȨ��
    //����
    while (::GetElapsedProcessTime() - InitTime < m_TimeLimit) {
        //Эͬ����--��
        cout<<"group"<<endl;
        for (int PS = 0; PS < m_PS1; PS++) {
            m_SpeedMatrix = m_ArchiveSpeedVector[PS];
            //�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
            for (int j = 0; j < m_Jobs; j++) {
                for (int i = 0; i < m_Machines; i++) {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }
            FacFamSeq = m_FacFamSeqArray[PS];
            int Map1 = rand() % m_RefSize;// �����ѡһ�������ߣ��������У��ӵ������� form a solution by randomly selecting a RefJobSeq;
            vector<int> FacSpan(m_Factories, 0);
            vector<float> FacEC(m_Factories, 0);
            int ObjectMS;
            float ObjectTEC;

            GetJFDHierarchy_Forward_New(FacFamSeq, m_ArchiveJobSeqinFamArray[Map1], m_TempJFDTime,m_TempHierarchy);
            GetMSandTECForPerandToalFacByJFD(FacFamSeq, m_ArchiveJobSeqinFamArray[Map1], m_TempJFDTime, FacSpan, FacEC,
                                             ObjectMS, ObjectTEC);
            int OrgMS = ObjectMS;
            float OrgTEC = ObjectTEC;
            int NewState = 0;
            int OldState = 0;

            //ʹ��q-learning����������ѡ�����
            int Act = DetermineAction(OldState, QFam[PS]);
            PerformAction_Fams(Act, FacFamSeq, m_ArchiveJobSeqinFamArray[Map1], m_TempJFDTime, m_TempHierarchy,
                               m_QCCEAPopulation, FacSpan, FacEC, ObjectMS, ObjectTEC);
            double Reward = CaculateReward(OrgMS, ObjectMS, OrgTEC, ObjectTEC, m_QCCEAPopulation);

            NewState = DetermineState(ObjectMS, ObjectTEC, m_QCCEAPopulation);
            UpdateQ(OldState, NewState, Act, Reward, QFam[PS]);
            OldState = NewState;
        }
        cout<<"group end "<<endl;
       cout<<"job"<<endl;
        //Эͬ����--����
        for (int PS = 0; PS < m_PS2; PS++) {
            m_SpeedMatrix = m_ArchiveSpeedVector[PS];
            //�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
            for (int j = 0; j < m_Jobs; j++) {
                for (int i = 0; i < m_Machines; i++) {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }

            JobSeqInFam = m_JobSeqinFamArray[PS];
            int Map2 = rand() % m_RefSize;// �����ѡһ�������ߣ������У� form a solution by randomly selecting a ReFacFamseq
            vector<int> FacSpan(m_Factories, 0);
            vector<float> FacTEC(m_Factories, 0);
            int ObjectMS = -1;
            float ObjectTEC = -1;
            GetJFDHierarchy_Forward_New(m_ArchiveFacFamSeqArray[Map2], JobSeqInFam, m_TempJFDTime, m_TempHierarchy);
            GetMSandTECForPerandToalFacByJFD(m_ArchiveFacFamSeqArray[Map2], JobSeqInFam, m_TempJFDTime, FacSpan, FacTEC,
                                             ObjectMS, ObjectTEC);

            JobLocalSearch(m_ArchiveFacFamSeqArray[Map2], JobSeqInFam, m_TempJFDTime, m_TempHierarchy,
                                              m_NadirPointMS, m_NadirPointTEC, m_IdealPointMS, m_IdealPointTEC,
                                              m_QCCEAPopulation, FacSpan, FacTEC, ObjectMS, ObjectTEC);

            BasedInd_SwapJob_New(m_ArchiveFacFamSeqArray[Map2], JobSeqInFam, m_TempJFDTime, m_TempHierarchy,
                                 m_NadirPointMS, m_NadirPointTEC, m_IdealPointMS, m_IdealPointTEC, m_QCCEAPopulation,
                                 FacSpan, FacTEC,
                                 ObjectMS, ObjectTEC);
        }
        cout<<"job end "<<endl;
        //Эͬ����--������
        cout<<"dangan "<<endl;
        for (int PS = 0; PS < m_RefSize; PS++) {

            m_SpeedMatrix = m_ArchiveSpeedVector[PS];
            //�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
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
            GetJFDTime_Forward(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], m_TempJFDTimePop[PS], FacSpan);
            GetJFDHierarchy_Forward_New(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS],
                                        m_TempJFDTimePop[PS],m_TempHierarchyPop[PS]);
            GetMSandTECForPerandToalFacByJFD(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], m_TempJFDTimePop[PS], FacSpan, FacTEC,
                                             ObjectMS, ObjectTEC);
            //����͹���ִ���ƻ��ع�����
            BasedInd_RefDeconAndCon_New(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], m_TempJFDTimePop[PS],
                                        m_TempHierarchyPop[PS],
                                        m_NadirPointMS, m_NadirPointTEC, m_IdealPointMS,
                                        m_IdealPointTEC, m_QCCEAPopulation, FacSpan, FacTEC,
                                        m_ArchiveSpanArray[PS], m_ArchiveTECArray[PS], scores_combined[PS]);

        }
        cout<<"dangan end "<<endl;
        //�Ӳο���ѡ��� ���µ��������飬���������ұ���
        UpdateArchiveGroupJobSet(mu);

        m_IterNum++;

    }

    //���ܺͽ���MS
    SaveTECandDeMS(m_QCCEAPopulation);
    cout<<"Evoltuon end "<<endl;
}


int QCCEA::RunEvolution(int CPUTime, vector<vector<Individual>> &QCCEAFinalAfterRepParetoSet, int AN, int mu) {
    ReadInstanceFileNameList("..\\MILPInstances\\");
    int Instances = 405;
    double T = 0.8;
    int Reps = 10; //���ÿ�������ظ����еĴ���

    string FileDirectory = "..\\Result\\";
    ostringstream str;
    str << "QCCEA_" << CPUTime << "_experiment" << ".txt"; //��ͬ���㷨
    ofstream ofile;
    ofile.open(FileDirectory + str.str());

//    ofstream ofileIter;
//    ostringstream strIter;
//    strIter << "QCCEA_" << CPUTime << "_iterNum" << ".txt";
//    ofileIter.open(FileDirectory + strIter.str());

    for (int ins = 4; ins < Instances; ins++) {

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

        //��
        vector<Individual> FinalQCCEAParetoSet;
        vector<Individual> AfterRepParetoSet;

        vector<Individual> TempAfterRepParetoSet;
        TempAfterRepParetoSet.clear();
        int totalIterNum = 0;
        for (int r = 0; r < Reps; r++) {

            long TimeLimit = CPUTime * m_Machines * m_Families; //original: 20 * m_Machines * m_Families
            this->SetParameters(AN,  TimeLimit, T);

            m_TempJFDTime.clear();
            m_TempJFDTime.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));

            m_tempJFDTime1.clear();
            m_tempJFDTime1.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));

            m_TempJFDTimePop.clear();
            m_TempJFDTimePop.resize(this->m_PopSize,
                                    vector<vector<int>>(this->m_Jobs, vector<int>(this->m_Machines, 0)));

            m_TempHierarchy.clear();
            m_TempHierarchy.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));

            m_tempHierarchy1.clear();
            m_tempHierarchy1.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));

            m_TempHierarchyPop.clear();
            m_TempHierarchyPop.resize(this->m_PopSize,
                                      vector<vector<int>>(this->m_Jobs, vector<int>(this->m_Machines, 0)));

            //���ٺ�Ĵ���ʱ�䣨�����Ĵ���ʱ�䣩
            m_TureJobOpertime.clear();
            m_TureJobOpertime.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));

            //�ٶȾ���
            m_SpeedMatrix.clear();
            m_SpeedMatrix.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));

            //����ʱ��P�ĵ�λ�ܺ�
            UnitPEC.clear();
            UnitPEC.resize(this->m_Jobs, vector<float>(this->m_Machines, 0));


            MachReadyTime.clear();
            MachReadyTime.resize(this->m_Machines);

            //��ʼ����Ⱥ
            InitialPop(m_TempJFDTimePop, m_TempHierarchyPop);

            //��������
            EvolutionProcess(mu);

            totalIterNum += m_IterNum;

            //��֧���
            m_QCCEAParetoSet.clear();

            Pareto_relation(m_QCCEAPopulation);

            //����֧��� de-emphasize dominated solutions
            for (int j = 0; j < m_QCCEAPopulation.size(); j++) {
                for (auto & i : m_QCCEAPopulation) {
                    if (i.flag == 0) {
                        if (i.pareto_rel[j] == 1) {
                            i.flag = 999;
                        }
                    }
                }
            }

            for (auto & i : m_QCCEAPopulation) {
                if (i.flag == 0)
                    m_QCCEAParetoSet.push_back(i);
            }

            //ȥ���ظ�
            FinalQCCEAParetoSet.clear();
            for (auto & i : m_QCCEAParetoSet) {
                bool fg = true;
                for (auto & j : FinalQCCEAParetoSet) {
                    if ((j.MS == i.MS) &&
                        (j.TEC == i.TEC)) {
                        fg = false;
                        break;
                    }
                }
                if (fg) {
                    FinalQCCEAParetoSet.push_back(i);
                }

            }

            long EndTime_IG = GetElapsedProcessTime();

            for (auto & PS : FinalQCCEAParetoSet) {
                TempAfterRepParetoSet.push_back(PS);
            }
        }//end rep

     //   ofileIter  << totalIterNum << endl;

        //��֧���
        AfterRepParetoSet.clear();

        Pareto_relation(TempAfterRepParetoSet);

        //����֧��� de-emphasize dominated solutions
        for (int j = 0; j < TempAfterRepParetoSet.size(); j++) {
            for (auto & i : TempAfterRepParetoSet) {
                if (i.flag == 0) {
                    if (i.pareto_rel[j] == 1) {
                        i.flag = 999;
                    }
                }
            }
        }

        for (auto & i : TempAfterRepParetoSet) {
            if (i.flag == 0)
                AfterRepParetoSet.push_back(i);
        }

        //ȥ���ظ�
        vector<Individual> FinalAfterRepParetoSet;
        FinalAfterRepParetoSet.clear();
        for (auto & i : AfterRepParetoSet) {
            bool fg = true;
            for (auto & j : FinalAfterRepParetoSet) {
                if ((j.MS == i.MS) &&
                    (j.TEC == i.TEC)) {
                    fg = false;
                    break;
                }
            }
            if (fg) {
                FinalAfterRepParetoSet.push_back(i);
            }
        }

        QCCEAFinalAfterRepParetoSet[ins].clear();
        cout << ins + 1 << "\t" << "Factories:" << this->m_Factories << "\t" << "Machines :" << this->m_Machines
             << "\t" << "Families:" << this->m_Families << "\t" << "Jobs :" << this->m_Jobs << "\t";

        for (auto & PS : FinalAfterRepParetoSet) {
            // ��� MS �� TEC�����нⶼ����ͬһ��
            cout << "MS:" << PS.MS
                 << " " << "TEC:" << PS.TEC << "\t";
            ofile << PS.MS << "  " << PS.TEC << "," << "\t";
            QCCEAFinalAfterRepParetoSet[ins].push_back(PS);
        }
        ofile << endl;  // ÿ��ʵ�������������
        cout << endl;   // ����̨������У���������ʵ������ʱ����

    }//end ins
    ofile.close();
    return 0;
}
