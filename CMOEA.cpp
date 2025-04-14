#include "CMOEA.h"
using namespace Base;

CMOEA::CMOEA()
{
}

CMOEA::~CMOEA()
{
}

void CMOEA::SetParameters(int AN, int NoImproveNumber, int NumberOfExtractedFams, long TimeLimit, int AllJobTotalPTime)
{
    m_BlockLengthChoice = AN;
    m_NoImproveNumber = NoImproveNumber;
    this->m_FamD = NumberOfExtractedFams;
    this->m_TimeLimit = TimeLimit;
    this->m_Temperature = 0.6 * AllJobTotalPTime / (this->m_Factories * this->m_Jobs * this->m_Machines);

    this->m_Popsize = 1000;
    this->m_RefSize = AN;
    this->m_PS1 = AN;
    this->m_PS2 = AN;
    this->thr_zeta = 1.0;

    CMOEAPopulation.clear();
    CMOEAPopulation.resize(m_RefSize);

    m_nadirpointMS = -1;
    m_nadirpointTEC = -1;

    m_idealpointMS = -1;
    m_idealpointTEC = -1;

    //������ Reference set
    m_RefSpanArray.clear();//�������깤ʱ��
    m_RefSpanArray.resize(m_RefSize);

    m_RefTECArray.clear();//���������ܺ�
    m_RefTECArray.resize(m_RefSize);

    m_RefSpeedVector.clear();
    m_RefSpeedVector.resize(m_RefSize);

    m_nRefCriFacArray.clear();//�������ؼ�����
    m_nRefCriFacArray.resize(m_RefSize, 0);

    m_RefFacFamSeqArray.clear();//������������
    m_RefFacFamSeqArray.resize(m_RefSize);

    m_RefJobSeqinFamArray.clear();//�������鹤��
    m_RefJobSeqinFamArray.resize(m_RefSize);

    m_RefFacSpanArray.clear();//�����깤ʱ��
    m_RefFacSpanArray.resize(m_RefSize);

    m_RefFacECArray.clear();//�����ܺ�
    m_RefFacECArray.resize(m_RefSize);

    m_bFlag1.clear();//
    m_bFlag1.resize(m_RefSize);

    m_bFlag2.clear();//
    m_bFlag2.resize(m_RefSize);

    //����Ⱥ
    m_SpanArray1.clear();//�깤ʱ��
    m_SpanArray1.resize(m_PS1);

    m_TECArray1.clear(); //�ܺ�
    m_TECArray1.resize(m_PS1);

    m_nCriFacArray1.clear();//�ؼ�����
    m_nCriFacArray1.resize(m_PS1, 0);

    m_Map1.clear();//Э�����ڵ��������±�
    m_Map1.resize(m_PS1);

    m_Age1.clear();//����
    m_Age1.resize(m_PS1);

    m_FacFamSeqArray.clear();//������
    m_FacFamSeqArray.resize(m_PS1);

    //������Ⱥ
    m_SpanArray2.clear();//�깤ʱ��
    m_SpanArray2.resize(m_PS2);

    m_TECArray2.clear(); //�ܺ�
    m_TECArray2.resize(m_PS2);

    m_nCriFacArray2.clear();//�ؼ�����
    m_nCriFacArray2.resize(m_PS2, 0);

    m_Map2.clear();//Э�����ڵ��������±�
    m_Map2.resize(m_PS2);

    m_Age2.clear();//����
    m_Age2.resize(m_PS2);

    m_JobSeqinFamArray.clear();//��������
    m_JobSeqinFamArray.resize(m_PS2);
}


void CMOEA::BasedindRandFamInFacTobestPos(vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, vector<vector<int>>& JFDTime, vector<vector<int>>& JBDTime,
                                          int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual>& CMOEAPopulation, int& ObjectMS, float& ObjectTEC)
{
    vector<int> FacSpan(m_Factories, 0);

    GetJFDTime_Forward(FacFamSeq, JobSeqInFam, JFDTime, FacSpan);
    GetJBDTime_Backward(FacFamSeq, JobSeqInFam, JBDTime, FacSpan);

    vector<int> SpanFac(FacFamSeq.size(), 0);
    vector<float> ECFac(FacFamSeq.size(), 0);
    int OrgMS = -1;
    float OrgTEC = -1;
    GetMSandTECForPerandToalFacByJFD(FacFamSeq, JobSeqInFam, JFDTime, SpanFac, ECFac, OrgMS, OrgTEC);

    int fac = -1;
    int min = -1;

    //�ؼ�����
    for (int i = 0; i < FacFamSeq.size(); i++)
    {
        if (min < FacSpan[i])
        {
            min = FacSpan[i];
            fac = i;
        }

    }

    int BestFac;
    int BestPos;

    float Ind = -1;

    vector<vector<int>> TempFacFamSeq;
    TempFacFamSeq.clear();
    TempFacFamSeq = FacFamSeq;

    for (int f = 0; f < TempFacFamSeq[fac].size(); f++)
    {
        if (FacFamSeq[fac].size() < 2)
            continue;

        int CurFam = TempFacFamSeq[fac][f];

        vector<int>::iterator it = find(FacFamSeq[fac].begin(), FacFamSeq[fac].end(), CurFam); //���������
        int j = it - FacFamSeq[fac].begin();
        FacFamSeq[fac].erase(it);

        int BackwardJobPos;
        if (j - 1 == -1)
        {
            BackwardJobPos = -1;
        }
        else
        {
            BackwardJobPos = JobSeqInFam[FacFamSeq[fac][j - 1]].size() - 1;
        }

        RefreshJFDJBDTime_InFactory(FacFamSeq[fac], JobSeqInFam, JFDTime, JBDTime, j, j - 1, 0, BackwardJobPos);

        BestFac = -1;
        BestPos = -1;
        Ind = this->FindBestPosToInsertFamForAllFac_Ind(FacFamSeq, JobSeqInFam, JFDTime, JBDTime, CurFam,
                                                        BestFac, BestPos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, OrgMS, OrgTEC, CMOEAPopulation);
        FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
        RefreshJFDJBDTime_InFactory(FacFamSeq[BestFac], JobSeqInFam, JFDTime, JBDTime, BestPos, BestPos, 0, JobSeqInFam[CurFam].size() - 1);

        vector<int> AfterInsertSpanFac(FacFamSeq.size(), 0);
        vector<float> AfterInsertECFac(FacFamSeq.size(), 0);
        int AfterInsertMS = -1;
        float AfterInsertTEC = -1;
        GetMSandTECForPerandToalFacByJFD(FacFamSeq, JobSeqInFam, JFDTime, AfterInsertSpanFac, AfterInsertECFac, AfterInsertMS, AfterInsertTEC);

        ObjectMS = AfterInsertMS;
        ObjectTEC = AfterInsertTEC;

        //��麯��
        //CheckSol(FacFamSeq, JobSeqInFam, AfterInsertMS);
        //CheckSolTEC(FacFamSeq, JobSeqInFam, AfterInsertTEC);
    }
}

void CMOEA::BasedindSwapFam(vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, vector<vector<int>>& JFDTime, vector<vector<int>>& JBDTime,
                            int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual>& CMOEAPopulation, int& ObjectMS, float& ObjectTEC)
{

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
    float Orgconvergence_ind = 0;
    Orgconvergence_ind += (OrgnormalMS - 1.0) * (OrgnormalMS - 1.0);
    Orgconvergence_ind += (OrgnormalTEC - 1.0) * (OrgnormalTEC - 1.0);
    Orgconvergence_ind = sqrt(Orgconvergence_ind);
    Orgconvergence_ind = 1 / (Orgconvergence_ind + 1);
    //cout << "���������ָ�꣺" << Orgconvergence_ind << endl;
    //cout << endl;

    //cout << "����ǰ��" << endl;
    vector<int> FacSpan(FacFamSeq.size(), 0);
    for (int Fac = 0; Fac < m_Factories; Fac++)
    {
        int LastFam = *FacFamSeq[Fac].rbegin();
        int LastJobInLstFam = *JobSeqInFam[LastFam].rbegin(); //������һ������
        FacSpan[Fac] = JFDTime[LastJobInLstFam][this->m_Machines - 1];
        //cout << "����" << Fac << "��Span��" << FacSpan[Fac] << endl;
    }

    vector<float> FacEC(FacFamSeq.size(), 0);

    //�����ܺ�
    for (int fac = 0; fac < FacFamSeq.size(); fac++)
    {
        //�����ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
            //int CurFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
            int Fam = FacFamSeq[fac][0];
            for (int g = 0; g < FacFamSeq[fac].size(); g++)
            {
                Fam = FacFamSeq[fac][g];

                if (preFam == -1)
                {
                    TSetupTimeEC += m_SetupTime[m][Fam][Fam] * UnitSetupEC;
                    Tsetuptime += m_SetupTime[m][Fam][Fam];
                }
                else
                {
                    TSetupTimeEC += m_SetupTime[m][preFam][Fam] * UnitSetupEC;
                    Tsetuptime += m_SetupTime[m][preFam][Fam];
                }
                preFam = Fam;
                for (int j = 0; j < JobSeqInFam[Fam].size(); j++)
                {
                    int Job = JobSeqInFam[Fam][j];
                    TpTimeEC += m_TureJobOpertime[Job][m] * UnitPEC[Job][m];
                    Tptime += m_TureJobOpertime[Job][m];
                }

            }
            Tidletime += (FacSpan[fac] - Tptime - Tsetuptime);
            TIdleTimeEC += Tidletime * UnitIdleEC;
        }
        FacEC[fac] = TSetupTimeEC + TIdleTimeEC + TpTimeEC;
        //cout << "����" << fac << "���ܺģ�" << FacEC[fac] << endl;
    }


    for (int i = 0; i < (m_Families / m_Factories); i++)
    {

        CopyJFDJBDTime(JFDTime, JBDTime, m_tempJFDTime1, m_tempJBDTime1);


        if (i != 0)
        {
            for (int j = 0; j < FacFamSeq.size(); j++)
                FacSpan[j] = 0;

            for (int Fac = 0; Fac < m_Factories; Fac++)
            {
                int LastFam = *FacFamSeq[Fac].rbegin();
                int LastJobInLstFam = *JobSeqInFam[LastFam].rbegin(); //������һ������
                FacSpan[Fac] = JFDTime[LastJobInLstFam][this->m_Machines - 1];
                //cout << "����" << Fac << "��Span��" << FacSpan[Fac] << endl;
            }
        }

        int cirfac = -1;
        int min = -1;
        int optfac = -1;
        int max = INT_MAX;

        //�ؼ�������makespan��С�Ĺ���
        for (int i = 0; i < FacFamSeq.size(); i++)
        {
            if (min < FacSpan[i])
            {
                min = FacSpan[i];
                cirfac = i;
            }
        }

        for (int i = FacFamSeq.size() - 1; i >= 0; i--)
        {
            if (max > FacSpan[i])
            {
                max = FacSpan[i];
                optfac = i;
            }
        }

        //���ѡ��һ����ӹ�����
        int Pos1 = rand() % FacFamSeq[cirfac].size();
        int Fam1 = FacFamSeq[cirfac][Pos1];
        FacFamSeq[cirfac].erase(FacFamSeq[cirfac].begin() + Pos1);

        int BackwardJobPos;
        if (Pos1 - 1 == -1)
        {
            BackwardJobPos = -1;
        }
        else
        {
            BackwardJobPos = JobSeqInFam[FacFamSeq[cirfac][Pos1 - 1]].size() - 1;
        }
        RefreshJFDJBDTime_InFactory(FacFamSeq[cirfac], JobSeqInFam, JFDTime, JBDTime, Pos1, Pos1 - 1, 0, BackwardJobPos);

        int Pos2 = rand() % FacFamSeq[optfac].size();
        int Fam2 = FacFamSeq[optfac][Pos2];
        FacFamSeq[optfac].erase(FacFamSeq[optfac].begin() + Pos2);
        if (Pos2 - 1 == -1)
        {
            BackwardJobPos = -1;
        }
        else
        {
            BackwardJobPos = JobSeqInFam[FacFamSeq[optfac][Pos2 - 1]].size() - 1;
        }
        RefreshJFDJBDTime_InFactory(FacFamSeq[optfac], JobSeqInFam, JFDTime, JBDTime, Pos2, Pos2 - 1, 0, BackwardJobPos);

        //��������optfac����������
        vector<int> ForwardMachReadyTime(this->m_Machines, 0);

        if (Pos2 == 0)
        {
            for (int m = 0; m < this->m_Machines; m++)
            {
                ForwardMachReadyTime[m] = this->m_SetupTime[m][Fam1][Fam1]; //MachReadyTime ��Ϊ׼��ʱ��
            }
        }
        else
        {
            int PreFam = FacFamSeq[optfac][Pos2 - 1]; //ǰһ����
            int LastJobinPreFam = JobSeqInFam[PreFam][JobSeqInFam[PreFam].size() - 1]; //ǰһ��������һ������
            for (int m = 0; m < this->m_Machines; m++)
            {
                //ǰһ�������һ���������깤ʱ��+׼��ʱ��
                ForwardMachReadyTime[m] = JFDTime[LastJobinPreFam][m] + this->m_SetupTime[m][PreFam][Fam1];
            }
        }

        // ��������﹤������ Scheduling jobs in InsFam
        for (int j = 0; j < JobSeqInFam[Fam1].size(); j++)
        {
            int CurJob = JobSeqInFam[Fam1][j]; //��ǰ����
            ForwardMachReadyTime[0] = std::max(ForwardMachReadyTime[0] + this->m_TureJobOpertime[CurJob][0], ForwardMachReadyTime[1]);

            for (int m = 1; m < this->m_Machines; m++)
            {
                if (m == this->m_Machines - 1)
                {
                    ForwardMachReadyTime[m] =ForwardMachReadyTime[m - 1] + this->m_TureJobOpertime[CurJob][m];
                }
                else
                {
                    ForwardMachReadyTime[m] = std::max( ForwardMachReadyTime[m - 1] + this-> m_TureJobOpertime[CurJob][m], ForwardMachReadyTime[m + 1]);
                }
            }
        }

        int Span2 = 0;
        if (Pos2 < FacFamSeq[optfac].size())
        {
            int NextFam = FacFamSeq[optfac][Pos2]; //��һ����

            int FirstJobinNextFam = JobSeqInFam[NextFam][0]; //��һ����ĵ�һ������

            for (int m = 0; m < this->m_Machines; m++)
            {
                if (ForwardMachReadyTime[m] + this->m_SetupTime[m][Fam1][NextFam] + JFDTime[FirstJobinNextFam][m] >= Span2)
                {
                    Span2 = ForwardMachReadyTime[m] + this->m_SetupTime[m][Fam1][NextFam] + JBDTime[FirstJobinNextFam][m];
                }
            }
        }
        else
        {
            Span2 = ForwardMachReadyTime[this->m_Machines - 1];
        }


        //��cirfac������
        vector<int> ForwardMachReadyTime2(this->m_Machines, 0);

        if (Pos1 == 0)
        {
            for (int m = 0; m < this->m_Machines; m++)
            {
                ForwardMachReadyTime2[m] = this->m_SetupTime[m][Fam2][Fam2]; //MachReadyTime ��Ϊ׼��ʱ��
            }
        }
        else
        {
            int PreFam = FacFamSeq[cirfac][Pos1 - 1]; //ǰһ����
            int LastJobinPreFam = JobSeqInFam[PreFam][JobSeqInFam[PreFam].size() - 1]; //ǰһ��������һ������
            for (int m = 0; m < this->m_Machines; m++)
            {
                //ǰһ�������һ���������깤ʱ��+׼��ʱ��
                ForwardMachReadyTime2[m] = JFDTime[LastJobinPreFam][m] + this->m_SetupTime[m][PreFam][Fam2];
            }
        }

        // ��������﹤������ Scheduling jobs in InsFam
        for (int j = 0; j < JobSeqInFam[Fam2].size(); j++)
        {
            int CurJob = JobSeqInFam[Fam2][j]; //��ǰ����
            ForwardMachReadyTime2[0] = std::max(ForwardMachReadyTime2[0] + this->m_TureJobOpertime[CurJob][0], ForwardMachReadyTime2[1]);

            for (int m = 1; m < this->m_Machines; m++)
            {
                if (m == this->m_Machines - 1)
                {
                    ForwardMachReadyTime2[m] =ForwardMachReadyTime2[m - 1] + this->m_TureJobOpertime[CurJob][m];
                }
                else
                {
                    ForwardMachReadyTime2[m] = std::max( ForwardMachReadyTime2[m - 1] + this-> m_TureJobOpertime[CurJob][m], ForwardMachReadyTime2[m + 1]);
                }
            }
        }

        int Span1 = 0;
        if (Pos1 < FacFamSeq[cirfac].size())
        {
            int NextFam = FacFamSeq[cirfac][Pos1]; //��һ����

            int FirstJobinNextFam = JobSeqInFam[NextFam][0]; //��һ����ĵ�һ������

            for (int m = 0; m < this->m_Machines; m++)
            {
                if (ForwardMachReadyTime2[m] + this->m_SetupTime[m][Fam2][NextFam] + JBDTime[FirstJobinNextFam][m] >= Span1)
                {
                    Span1 = ForwardMachReadyTime2[m] + this->m_SetupTime[m][Fam2][NextFam] + JBDTime[FirstJobinNextFam][m];
                }
            }
        }
        else
        {
            Span1 = ForwardMachReadyTime2[this->m_Machines - 1];
        }


        //�������ֵ
        FacSpan[cirfac] = Span1;
        FacSpan[optfac] = Span2;
        int AfterMakespan = *max_element(FacSpan.begin(), FacSpan.end());

        //�ж��Ƿ��滻��͵�������
        if (AfterMakespan > nadirpointMS)
            nadirpointMS = AfterMakespan;
        if (AfterMakespan < idealpointMS)
            idealpointMS = AfterMakespan;


        float AfterTEC = 0;
        //vector<float> FacEC(FacFamSeq.size(), 0);

        //�����ܺ�
        vector<vector<int>> TempFacFamSeq;
        TempFacFamSeq.clear();
        TempFacFamSeq.resize(FacFamSeq.size());
        TempFacFamSeq = FacFamSeq;
        TempFacFamSeq[optfac].insert(TempFacFamSeq[optfac].begin() + Pos2, Fam1);
        TempFacFamSeq[cirfac].insert(TempFacFamSeq[cirfac].begin() + Pos1, Fam2);

        //���㹤��2�ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
            //int CurFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
            int Fam = TempFacFamSeq[optfac][0];
            for (int g = 0; g < TempFacFamSeq[optfac].size(); g++)
            {
                Fam = TempFacFamSeq[optfac][g];

                if (preFam == -1)
                {
                    TSetupTimeEC += m_SetupTime[m][Fam][Fam] * UnitSetupEC;
                    Tsetuptime += m_SetupTime[m][Fam][Fam];
                }
                else
                {
                    TSetupTimeEC += m_SetupTime[m][preFam][Fam] * UnitSetupEC;
                    Tsetuptime += m_SetupTime[m][preFam][Fam];
                }
                preFam = Fam;
                for (int j = 0; j < JobSeqInFam[Fam].size(); j++)
                {
                    int Job = JobSeqInFam[Fam][j];
                    TpTimeEC += m_TureJobOpertime[Job][m] * UnitPEC[Job][m];
                    Tptime += m_TureJobOpertime[Job][m];
                }

            }
            Tidletime += (FacSpan[optfac] - Tptime - Tsetuptime);
            TIdleTimeEC += Tidletime * UnitIdleEC;
        }
        FacEC[optfac] = TSetupTimeEC + TIdleTimeEC + TpTimeEC;
        //cout << "�����󹤳�" << optfac << "���ܺģ�" << FacEC[optfac] << endl;

        //���㹤��1�ܺ�
        TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
            //int CurFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
            int Fam = TempFacFamSeq[cirfac][0];
            for (int g = 0; g < TempFacFamSeq[cirfac].size(); g++)
            {
                Fam = TempFacFamSeq[cirfac][g];

                if (preFam == -1)
                {
                    TSetupTimeEC += m_SetupTime[m][Fam][Fam] * UnitSetupEC;
                    Tsetuptime += m_SetupTime[m][Fam][Fam];
                }
                else
                {
                    TSetupTimeEC += m_SetupTime[m][preFam][Fam] * UnitSetupEC;
                    Tsetuptime += m_SetupTime[m][preFam][Fam];
                }
                preFam = Fam;
                for (int j = 0; j < JobSeqInFam[Fam].size(); j++)
                {
                    int Job = JobSeqInFam[Fam][j];
                    TpTimeEC += m_TureJobOpertime[Job][m] * UnitPEC[Job][m];
                    Tptime += m_TureJobOpertime[Job][m];
                }

            }
            Tidletime += (FacSpan[cirfac] - Tptime - Tsetuptime);
            TIdleTimeEC += Tidletime * UnitIdleEC;
        }
        FacEC[cirfac] = TSetupTimeEC + TIdleTimeEC + TpTimeEC;
        //cout << "�����󹤳�" << cirfac << "���ܺģ�" << FacEC[cirfac] << endl;

        AfterTEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);
        //cout << "���ܺģ�" << AfterTEC << endl;

        //�ж��Ƿ��滻��͵�������
        if (AfterTEC > nadirpointTEC)
            nadirpointTEC = AfterTEC;
        if (AfterTEC < idealpointTEC)
            idealpointTEC = AfterTEC;

        //�ж��Ƿ�Ȳ���ǰ�Ľ������Ľ������ο���
        bool flag = true;
        if (((AfterMakespan < OrgMS) && (AfterTEC < OrgTEC)) || ((AfterMakespan < OrgMS) && (AfterTEC == OrgTEC)) || ((AfterMakespan == OrgMS) && (AfterTEC < OrgTEC)) || (AfterMakespan < idealpointMS) || (AfterTEC < idealpointTEC))
        {
            for (int i = 0; i < CMOEAPopulation.size(); i++)
            {
                if ((AfterMakespan == CMOEAPopulation[i].MS) && (AfterTEC == CMOEAPopulation[i].TEC))
                {
                    flag = false;
                    break;
                }
            }
            if (flag)
            {
                //cout << "*********************" << endl;
                Individual tempIndi;
                tempIndi.m_FacFamSeqArray = TempFacFamSeq;
                tempIndi.m_JobSeqInFamArray = JobSeqInFam;
                tempIndi.MS = AfterMakespan;
                tempIndi.TEC = AfterTEC;
                tempIndi.m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���
                CMOEAPopulation.push_back(tempIndi);

            }
        }

        //��һ��
        float AfternormalMS = -1;
        float AfternormalTEC = -1;
        //cout << endl << "normalize" << endl;

        AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
        AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
        //cout << "�����һ�����MS��" << AfternormalMS << "\tTEC��" << AfternormalTEC << endl;

        //����ָ��Ind
        float Afterconvergence_ind = 0;
        Afterconvergence_ind += (AfternormalMS - 1.0) * (AfternormalMS - 1.0);
        Afterconvergence_ind += (AfternormalTEC - 1.0) * (AfternormalTEC - 1.0);
        Afterconvergence_ind = sqrt(Afterconvergence_ind);
        Afterconvergence_ind = 1 / (Afterconvergence_ind + 1);
        //cout << "������ĸ��������ָ�꣺" << Afterconvergence_ind << endl;
        //cout << endl;

        if (Afterconvergence_ind < Orgconvergence_ind)
        {
            FacFamSeq[optfac].insert(FacFamSeq[optfac].begin() + Pos2, Fam1);
            RefreshJFDJBDTime_InFactory(FacFamSeq[optfac], JobSeqInFam, JFDTime, JBDTime, Pos2, Pos2, 0, JobSeqInFam[Fam1].size() - 1);

            FacFamSeq[cirfac].insert(FacFamSeq[cirfac].begin() + Pos1, Fam2);
            RefreshJFDJBDTime_InFactory(FacFamSeq[cirfac], JobSeqInFam, JFDTime, JBDTime, Pos1, Pos1, 0, JobSeqInFam[Fam2].size() - 1);

            OrgMS = AfterMakespan;
            OrgTEC = AfterTEC;
            Orgconvergence_ind = Afterconvergence_ind;
            ObjectMS = AfterMakespan;
            ObjectTEC = AfterTEC;
        }
        else
        {
            FacFamSeq[optfac].insert(FacFamSeq[optfac].begin() + Pos2, Fam2);
            FacFamSeq[cirfac].insert(FacFamSeq[cirfac].begin() + Pos1, Fam1);
            CopyJFDJBDTime(m_tempJFDTime1, m_tempJBDTime1, JFDTime, JBDTime);
        }

        //���
        /*cout << "swapFam��" << endl;
        CheckSol(FacFamSeq, JobSeqInFam, ObjectMS);
        CheckSolTEC(FacFamSeq, JobSeqInFam, ObjectTEC);*/
    }

}

void CMOEA::BasedindCirJobInFamTobestPos(const vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, vector<vector<int>>& JFDTime, vector<vector<int>>& JBDTime,
                                         int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual>& CMOEAPopulation, int& ObjectMS, float& ObjectTEC)
{
    vector<int> FacSpan(m_Factories, 0);

    GetJFDTime_Forward(FacFamSeq, JobSeqInFam, JFDTime, FacSpan);
    GetJBDTime_Backward(FacFamSeq, JobSeqInFam, JBDTime, FacSpan);

    vector<int> SpanFac(FacFamSeq.size(), 0);
    vector<float> ECFac(FacFamSeq.size(), 0);
    int OrgMS = -1;
    float OrgTEC = -1;
    GetMSandTECForPerandToalFacByJFD(FacFamSeq, JobSeqInFam, JFDTime, SpanFac, ECFac, OrgMS, OrgTEC);
    //cout << "����֮ǰ��MS��" << OrgMS << "\tTEC��" << OrgTEC << endl;

    int fac = -1;
    int min = -1;

    //�ؼ�����
    for (int i = 0; i < FacFamSeq.size(); i++)
    {
        if (min < FacSpan[i])
        {
            min = FacSpan[i];
            fac = i;
        }

    }

    for (int Fam = 0; Fam < FacFamSeq[fac].size(); Fam++)
    {
        int CurFam = FacFamSeq[fac][Fam]; //��ǰ��
        if (JobSeqInFam[CurFam].size() > 1) //�鹤��������1�����оֲ�����
        {
            Basedind_JobsInFam(fac, FacFamSeq, FacFamSeq[fac], JobSeqInFam, JFDTime, JBDTime, CurFam, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, CMOEAPopulation, OrgMS, OrgTEC, ObjectMS, ObjectTEC);

        }
    }
}

void CMOEA::BasedindSwapJob(const vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, vector<vector<int>>& JFDTime, vector<vector<int>>& JBDTime,
                            int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual>& CMOEAPopulation, int& ObjectMS, float& ObjectTEC)
{
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
    float Orgconvergence_ind = 0;
    Orgconvergence_ind += (OrgnormalMS - 1.0) * (OrgnormalMS - 1.0);
    Orgconvergence_ind += (OrgnormalTEC - 1.0) * (OrgnormalTEC - 1.0);
    Orgconvergence_ind = sqrt(Orgconvergence_ind);
    Orgconvergence_ind = 1 / (Orgconvergence_ind + 1);
    //cout << "���������ָ�꣺" << Orgconvergence_ind << endl;
    //cout << endl;


    vector<int> FacSpan(FacFamSeq.size(), 0);
    for (int Fac = 0; Fac < m_Factories; Fac++)
    {
        int LastFam = *FacFamSeq[Fac].rbegin();
        int LastJobInLstFam = *JobSeqInFam[LastFam].rbegin(); //������һ������
        FacSpan[Fac] = JFDTime[LastJobInLstFam][this->m_Machines - 1];
        //cout << "����" << Fac << "��Span��" << FacSpan[Fac] << endl;
    }

    vector<float> FacEC(FacFamSeq.size(), 0);

    //�����ܺ�
    for (int fac = 0; fac < FacFamSeq.size(); fac++)
    {
        //�����ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
            //int CurFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
            int Fam = FacFamSeq[fac][0];
            for (int g = 0; g < FacFamSeq[fac].size(); g++)
            {
                Fam = FacFamSeq[fac][g];

                if (preFam == -1)
                {
                    TSetupTimeEC += m_SetupTime[m][Fam][Fam] * UnitSetupEC;
                    Tsetuptime += m_SetupTime[m][Fam][Fam];
                }
                else
                {
                    TSetupTimeEC += m_SetupTime[m][preFam][Fam] * UnitSetupEC;
                    Tsetuptime += m_SetupTime[m][preFam][Fam];
                }
                preFam = Fam;
                for (int j = 0; j < JobSeqInFam[Fam].size(); j++)
                {
                    int Job = JobSeqInFam[Fam][j];
                    TpTimeEC += m_TureJobOpertime[Job][m] * UnitPEC[Job][m];
                    Tptime += m_TureJobOpertime[Job][m];
                }

            }
            Tidletime += (FacSpan[fac] - Tptime - Tsetuptime);
            TIdleTimeEC += Tidletime * UnitIdleEC;
        }
        FacEC[fac] = TSetupTimeEC + TIdleTimeEC + TpTimeEC;
        //cout << "����" << fac << "���ܺģ�" << FacEC[fac] << endl;
    }

    int fac = -1;
    int min = -1;

    //�ؼ�����
    for (int i = 0; i < FacFamSeq.size(); i++)
    {
        if (min < FacSpan[i])
        {
            min = FacSpan[i];
            fac = i;
        }

    }

    for (int Fam = 0; Fam < FacFamSeq[fac].size(); Fam++)
    {
        int CurFam = FacFamSeq[fac][Fam]; //��ǰ��
        if (JobSeqInFam[CurFam].size() > 1) //�鹤��������1�����оֲ�����
        {
            //Basedind_JobsInFam(fac, FacFamSeq, FacFamSeq[fac], JobSeqInFam, JFDTime, JBDTime, CurFam, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, CMOEAPopulation, OrgMS, OrgTEC, ObjectMS, ObjectTEC);
            //Basedind_JobsSwap(fac, FacFamSeq, FacFamSeq[fac], JobSeqInFam, JFDTime, JBDTime, CurFam, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, CMOEAPopulation, OrgMS, OrgTEC, ObjectMS, ObjectTEC);

            int pt1 = 0;
            int pt2 = 0;

            if (JobSeqInFam[CurFam].size() > 2)
            {

                pt1 = rand() % JobSeqInFam[CurFam].size();
                do
                {
                    pt2 = rand() % JobSeqInFam[CurFam].size();
                } while (pt1 == pt2);

                if (pt1 > pt2)
                {
                    int t = pt1;
                    pt1 = pt2;
                    pt2 = t;
                }
            }
            else if (JobSeqInFam[CurFam].size() == 2)
            {
                pt1 = 0;
                pt2 = 1;
            }

            CopyJFDJBDTime( JFDTime, JBDTime, m_tempJFDTime1, m_tempJBDTime1);
            vector<vector<int>> TempJobSeqInFam = JobSeqInFam;

            auto it = find(begin(FacFamSeq[fac]), end(FacFamSeq[fac]), CurFam);
            if (it != end(FacFamSeq[fac]))
            {
                int famPos = it - begin(FacFamSeq[fac]);
                int Job = JobSeqInFam[CurFam][pt1];
                JobSeqInFam[CurFam][pt1] = JobSeqInFam[CurFam][pt2];
                JobSeqInFam[CurFam][pt2] = Job;
                RefreshJFDJBDTime_InFactory(FacFamSeq[fac], JobSeqInFam, JFDTime, JBDTime, famPos, famPos, pt1, pt2);
            }

            FacSpan[fac] =GetSpanPerFacByJFD(FacFamSeq[fac], JobSeqInFam,  JFDTime);

            int AfterMakespan = *max_element(FacSpan.begin(), FacSpan.end());

            //�ж��Ƿ��滻��͵�������
            if (AfterMakespan > nadirpointMS)
                nadirpointMS = AfterMakespan;
            if (AfterMakespan < idealpointMS)
                idealpointMS = AfterMakespan;


            float AfterTEC = 0;
            //vector<float> FacEC(FacFamSeq.size(), 0);

            //�����ܺ�
            float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

            for (int m = 0; m < m_Machines; m++)
            {
                int preFam = -1;
                //int CurFam = -1;
                int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
                int Fam = FacFamSeq[fac][0];
                for (int g = 0; g < FacFamSeq[fac].size(); g++)
                {
                    Fam = FacFamSeq[fac][g];

                    if (preFam == -1)
                    {
                        TSetupTimeEC += m_SetupTime[m][Fam][Fam] * UnitSetupEC;
                        Tsetuptime += m_SetupTime[m][Fam][Fam];
                    }
                    else
                    {
                        TSetupTimeEC += m_SetupTime[m][preFam][Fam] * UnitSetupEC;
                        Tsetuptime += m_SetupTime[m][preFam][Fam];
                    }
                    preFam = Fam;
                    for (int j = 0; j < JobSeqInFam[Fam].size(); j++)
                    {
                        int Job = JobSeqInFam[Fam][j];
                        TpTimeEC += m_TureJobOpertime[Job][m] * UnitPEC[Job][m];
                        Tptime += m_TureJobOpertime[Job][m];
                    }

                }
                Tidletime += (FacSpan[fac] - Tptime - Tsetuptime);
                TIdleTimeEC += Tidletime * UnitIdleEC;
            }
            FacEC[fac] = TSetupTimeEC + TIdleTimeEC + TpTimeEC;
            //cout << "�ؼ�����" << fac << "���ܺģ�" << FacEC[fac] << endl;

            AfterTEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);
            //cout << "���ܺģ�" << AfterTEC << endl;

            //�ж��Ƿ��滻��͵�������
            if (AfterTEC > nadirpointTEC)
                nadirpointTEC = AfterTEC;
            if (AfterTEC < idealpointTEC)
                idealpointTEC = AfterTEC;

            //�ж��Ƿ�Ȳ���ǰ�Ľ������Ľ������ο���
            bool flag = true;
            if (((AfterMakespan < OrgMS) && (AfterTEC < OrgTEC)) || ((AfterMakespan < OrgMS) && (AfterTEC == OrgTEC)) || ((AfterMakespan == OrgMS) && (AfterTEC < OrgTEC)) || (AfterMakespan < idealpointMS) || (AfterTEC < idealpointTEC))
            {
                for (int i = 0; i < CMOEAPopulation.size(); i++)
                {
                    if ((AfterMakespan == CMOEAPopulation[i].MS) && (AfterTEC == CMOEAPopulation[i].TEC))
                    {
                        flag = false;
                        break;
                    }
                }
                if (flag)
                {
                    //cout << "*********************" << endl;
                    Individual tempIndi;
                    tempIndi.m_FacFamSeqArray = FacFamSeq;
                    tempIndi.m_JobSeqInFamArray = JobSeqInFam;
                    tempIndi.MS = AfterMakespan;
                    tempIndi.TEC = AfterTEC;
                    tempIndi.m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���
                    CMOEAPopulation.push_back(tempIndi);
                }
            }

            //��һ��
            float AfternormalMS = -1;
            float AfternormalTEC = -1;
            //cout << endl << "normalize" << endl;

            AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
            AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
            //cout << "�����һ�����MS��" << AfternormalMS << "\tTEC��" << AfternormalTEC << endl;

            //����ָ��Ind
            float Afterconvergence_ind = 0;
            Afterconvergence_ind += (AfternormalMS - 1.0) * (AfternormalMS - 1.0);
            Afterconvergence_ind += (AfternormalTEC - 1.0) * (AfternormalTEC - 1.0);
            Afterconvergence_ind = sqrt(Afterconvergence_ind);
            Afterconvergence_ind = 1 / (Afterconvergence_ind + 1);
            //cout << "������ĸ��������ָ�꣺" << Afterconvergence_ind << endl;
            //cout << endl;

            if (Afterconvergence_ind < Orgconvergence_ind)
            {
                OrgMS = AfterMakespan;
                OrgTEC = AfterTEC;
                Orgconvergence_ind = Afterconvergence_ind;
                ObjectMS = AfterMakespan;
                ObjectTEC = AfterTEC;
            }
            else
            {
                JobSeqInFam = TempJobSeqInFam;
                CopyJFDJBDTime(m_tempJFDTime1, m_tempJBDTime1, JFDTime, JBDTime);
            }
        }
        //��麯��
        //CheckSol(FacFamSeq, JobSeqInFam, ObjectMS);
        //CheckSolTEC(FacFamSeq, JobSeqInFam, ObjectTEC);
    }
}

void CMOEA::BasedindRefInsert(vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqInFam, vector<vector<int>>& JFDTime, vector<vector<int>>& JBDTime,
                              int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual>& CMOEAPopulation, int& ObjectMS, float& ObjectTEC)
{
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
    //cout << endl << "�ƻ��ع�ǰ��MS��" << OrgMS << "\tTEC��" << OrgTEC << endl;
    //cout << endl << "normalize" << endl;

    OrgnormalMS = (static_cast<float>(OrgMS - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    OrgnormalTEC = ((OrgTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
    //cout << "ԭʼ�����һ�����MS��" << OrgnormalMS << "\tTEC��" << OrgnormalTEC << endl;

    //���㽻��֮ǰָ��Ind
    float Orgconvergence_ind = 0;
    Orgconvergence_ind += (OrgnormalMS - 1.0) * (OrgnormalMS - 1.0);
    Orgconvergence_ind += (OrgnormalTEC - 1.0) * (OrgnormalTEC - 1.0);
    Orgconvergence_ind = sqrt(Orgconvergence_ind);
    Orgconvergence_ind = 1 / (Orgconvergence_ind + 1);
    //cout << "�ƻ��ع�ǰ���������ָ�꣺" << Orgconvergence_ind << endl;
    //cout << endl;

    vector<int> FamsExtracted;
    unordered_map<int, vector<int>> JobsExtracted;

    vector<vector<int>> TempFacFamSeq = FacFamSeq;
    vector<vector<int>> TempJobSeqInFam = JobSeqInFam;
    CopyJFDJBDTime(JFDTime, JBDTime, m_tempJFDTime1, m_tempJBDTime1);

    BasedindDestruction_FamsAndJobs(FacFamSeq, JobSeqInFam, JFDTime, JBDTime, FamsExtracted, JobsExtracted, 0);

    Basedind_Construction_FamsAndJobs(FacFamSeq, JobSeqInFam, JFDTime, JBDTime, FamsExtracted, JobsExtracted, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);

    vector<int> FacSpan(FacFamSeq.size(), 0);
    for (int Fac = 0; Fac < m_Factories; Fac++)
    {
        int LastFam = *FacFamSeq[Fac].rbegin();
        int LastJobInLstFam = *JobSeqInFam[LastFam].rbegin(); //������һ������
        FacSpan[Fac] = JFDTime[LastJobInLstFam][this->m_Machines - 1];
        //cout << "����" << Fac << "��Span��" << FacSpan[Fac] << endl;
    }

    int AfterMakespan = *max_element(FacSpan.begin(), FacSpan.end());
    //cout << "�ع����makespan��" << AfterMakespan << endl;

    //�ж��Ƿ��滻��͵�������
    if (AfterMakespan > nadirpointMS)
        nadirpointMS = AfterMakespan;
    if (AfterMakespan < idealpointMS)
        idealpointMS = AfterMakespan;

    float AfterTEC = -1;
    vector<float> FacEC(FacFamSeq.size(), 0);

    //�����ܺ�
    for (int fac = 0; fac < FacFamSeq.size(); fac++)
    {
        //�����ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
            //int CurFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
            int Fam = FacFamSeq[fac][0];
            for (int g = 0; g < FacFamSeq[fac].size(); g++)
            {
                Fam = FacFamSeq[fac][g];

                if (preFam == -1)
                {
                    TSetupTimeEC += m_SetupTime[m][Fam][Fam] * UnitSetupEC;
                    Tsetuptime += m_SetupTime[m][Fam][Fam];
                }
                else
                {
                    TSetupTimeEC += m_SetupTime[m][preFam][Fam] * UnitSetupEC;
                    Tsetuptime += m_SetupTime[m][preFam][Fam];
                }
                preFam = Fam;
                for (int j = 0; j < JobSeqInFam[Fam].size(); j++)
                {
                    int Job = JobSeqInFam[Fam][j];
                    TpTimeEC += m_TureJobOpertime[Job][m] * UnitPEC[Job][m];
                    Tptime += m_TureJobOpertime[Job][m];
                }

            }
            Tidletime += (FacSpan[fac] - Tptime - Tsetuptime);
            TIdleTimeEC += Tidletime * UnitIdleEC;
        }
        FacEC[fac] = TSetupTimeEC + TIdleTimeEC + TpTimeEC;
        //cout << "����" << fac << "���ܺģ�" << FacEC[fac] << endl;
    }

    AfterTEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);
    //cout << "�ع�������ܺģ�" << AfterTEC << endl;

    //�ж��Ƿ��滻��͵�������
    if (AfterTEC > nadirpointTEC)
        nadirpointTEC = AfterTEC;
    if (AfterTEC < idealpointTEC)
        idealpointTEC = AfterTEC;

    //�ж��Ƿ�Ȳ���ǰ�Ľ������Ľ������ο���
    bool flag = true;
    if (((AfterMakespan < OrgMS) && (AfterTEC < OrgTEC)) || ((AfterMakespan < OrgMS) && (AfterTEC == OrgTEC)) || ((AfterMakespan == OrgMS) && (AfterTEC < OrgTEC)) || (AfterMakespan < idealpointMS) || (AfterTEC < idealpointTEC))
    {
        for (int i = 0; i < CMOEAPopulation.size(); i++)
        {
            if ((AfterMakespan == CMOEAPopulation[i].MS) && (AfterTEC == CMOEAPopulation[i].TEC))
            {
                flag = false;
                break;
            }
        }
        if (flag)
        {
            //cout << "*********************" << endl;
            Individual tempIndi;
            tempIndi.m_FacFamSeqArray = FacFamSeq;
            tempIndi.m_JobSeqInFamArray = JobSeqInFam;
            tempIndi.MS = AfterMakespan;
            tempIndi.TEC = AfterTEC;
            tempIndi.m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���
            CMOEAPopulation.push_back(tempIndi);
        }
    }

    //��һ��
    float AfternormalMS = -1;
    float AfternormalTEC = -1;
    //cout << endl << "normalize" << endl;

    AfternormalMS = (static_cast<float>(AfterMakespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    AfternormalTEC = ((AfterTEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
    //cout << "�����һ�����MS��" << AfternormalMS << "\tTEC��" << AfternormalTEC << endl;

    //����ָ��Ind
    float Afterconvergence_ind = 0;
    Afterconvergence_ind += (AfternormalMS - 1.0) * (AfternormalMS - 1.0);
    Afterconvergence_ind += (AfternormalTEC - 1.0) * (AfternormalTEC - 1.0);
    Afterconvergence_ind = sqrt(Afterconvergence_ind);
    Afterconvergence_ind = 1 / (Afterconvergence_ind + 1);
    //cout << "�ع���ĸ��������ָ�꣺" << Afterconvergence_ind << endl;
    //cout << endl;

    if (Afterconvergence_ind < Orgconvergence_ind)
    {
        OrgMS = AfterMakespan;
        OrgTEC = AfterTEC;
        Orgconvergence_ind = Afterconvergence_ind;
        ObjectMS = AfterMakespan;
        ObjectTEC = AfterTEC;
    }
    else
    {
        FacFamSeq = TempFacFamSeq;
        JobSeqInFam = TempJobSeqInFam;
        CopyJFDJBDTime(m_tempJFDTime1, m_tempJBDTime1, JFDTime, JBDTime);
    }

    //���
    //CheckSol(FacFamSeq, JobSeqInFam, ObjectMS);
    //CheckSolTEC(FacFamSeq, JobSeqInFam, ObjectTEC);
}

void CMOEA::BasedindDestruction_FamsAndJobs(vector<vector<int>>& Sol, vector<vector<int>>& JobSeqInFam, vector<vector<int>>& JFDTime, vector<vector<int>>& JBDTime,
                                            vector<int>& FamsExtracted, unordered_map<int, vector<int>>& JobsExtracted, int FamsExtractedMethods)
{
    unordered_map<int, pair<int, int>> FamPosErasedFromFac;
    auto UpdateFamErasedFromFac = [&FamPosErasedFromFac](int Pos, int fac)
    {
        if (FamPosErasedFromFac.find(fac) == end(FamPosErasedFromFac))
        {
            FamPosErasedFromFac[fac] = { Pos, Pos - 1 };
        }
        else
        {
            if (Pos < FamPosErasedFromFac[fac].first)
            {
                FamPosErasedFromFac[fac].first = Pos;
            }
            if (Pos >= FamPosErasedFromFac[fac].second)
            {
                FamPosErasedFromFac[fac].second = Pos - 1;
            }
            else
            {
                FamPosErasedFromFac[fac].second = FamPosErasedFromFac[fac].second - 1;
            }
        }
    };

    FamsExtracted.clear();
    JobsExtracted.clear();

    //this->m_FamD = JobSeqInFam.size() / 5;
    //this->m_FamD = 2 + rand() % (7 - 2 + 1);
    this->m_FamD = 1;
    //this->m_dmin + rand() % (this->m_dmax - this->m_dmin + 1);
    //this->m_FamD = 5 + rand() % (7 - 5 + 1);
    //switch (FamsExtractedMethods)
    //{
    //case 0: //method_1�������ȡ�����͹�����
    //{

    //	for (int fac = 0; fac < m_Factories; fac++)
    //	{
    //		int Pos = rand() % Sol[fac].size();
    //		int FamExt = Sol[fac][Pos];
    //		Sol[fac].erase(Sol[fac].begin() + Pos);
    //		FamsExtracted.push_back(FamExt);
    //		UpdateFamErasedFromFac(Pos, fac);
    //	}
    //	while (FamsExtracted.size() < this->m_FamD)
    //	{
    //		int fac;
    //		do
    //		{
    //			fac = rand() % this->m_Factories;
    //		} while (Sol[fac].size() <= 1);

    //		int Pos = rand() % Sol[fac].size();
    //		int FamExt = Sol[fac][Pos];
    //		Sol[fac].erase(Sol[fac].begin() + Pos);
    //		FamsExtracted.push_back(FamExt);
    //		UpdateFamErasedFromFac(Pos, fac);
    //	}
    //	break;
    //	//while (FamsExtracted.size() < this->m_FamD)
    //	//{
    //	//	int fac;
    //	//	do
    //	//	{
    //	//		fac = rand() % this->m_Factories;
    //	//	} while (Sol[fac].size() == 0);
    //	//	int Pos = rand() % Sol[fac].size();
    //	//	int FamExt = Sol[fac][Pos];
    //	//	Sol[fac].erase(Sol[fac].begin() + Pos);
    //	//	FamsExtracted.push_back(FamExt);
    //	//	UpdateFamErasedFromFac(Pos, fac);
    //	//}
    //	//break;
    //}
    //}

    /*for (int fac = 0; fac < m_Factories; fac++)
    {
        int Pos = rand() % Sol[fac].size();
        int FamExt = Sol[fac][Pos];
        Sol[fac].erase(Sol[fac].begin() + Pos);
        FamsExtracted.push_back(FamExt);
        UpdateFamErasedFromFac(Pos, fac);
    }*/
    while (FamsExtracted.size() < this->m_FamD)
    {
        int fac;
        do
        {
            fac = rand() % this->m_Factories;
        } while (Sol[fac].size() <= 1);

        int Pos = rand() % Sol[fac].size();
        int FamExt = Sol[fac][Pos];
        Sol[fac].erase(Sol[fac].begin() + Pos);
        FamsExtracted.push_back(FamExt);
        UpdateFamErasedFromFac(Pos, fac);
    }

    //�ӹ�����������ȡ����

    for (int i = 0; i < FamsExtracted.size(); i++) //����������Ĺ�����ȡһ��
    {
        int Fam = FamsExtracted[i];
        if (JobSeqInFam[Fam].size() >= 3)
        {
            for (int j = 0; j <= JobSeqInFam[Fam].size() / 2; j++)
            {
                int JobPos = rand() % JobSeqInFam[Fam].size();
                int Job = JobSeqInFam[Fam][JobPos];
                JobSeqInFam[Fam].erase(JobSeqInFam[Fam].begin() + JobPos);
                //JobsExtracted[i].push_back(Job);
                JobsExtracted[Fam].push_back(Job);
            }
        }
    }
    for (int i = 0; i < Sol.size(); i++)
    {
        int ForwardFamPos = FamPosErasedFromFac[i].first;
        int BackwardFamPos = FamPosErasedFromFac[i].second;

        int FowardJobPos = 0;
        int BackwardJobPos;
        if (BackwardFamPos == -1)
        {
            BackwardJobPos = -1;
        }
        else
        {
            int BackwardFam = Sol[i][BackwardFamPos];
            BackwardJobPos = JobSeqInFam[BackwardFam].size() - 1;
        }

        RefreshJFDJBDTime_InFactory(Sol[i], JobSeqInFam, JFDTime, JBDTime, ForwardFamPos, BackwardFamPos, FowardJobPos, BackwardJobPos);

    }
}

void CMOEA::Basedind_Construction_FamsAndJobs(vector<vector<int>>& Sol, vector<vector<int>>& JobSeqInFam, vector<vector<int>>& JFDTime, vector<vector<int>>& JBDTime,
                                              vector<int>& FamsExtracted, unordered_map<int, vector<int>>& JobsExtracted, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC)
{
    //AllFacTT.resize(this->m_Factories, 0);
    //vector<int> tempFamsExtracted = FamsExtracted;
    while (FamsExtracted.size() > 0)
    {
        int BestFac = -1;
        int BestPos = -1;
        int Pos = rand() % FamsExtracted.size();
        int CurFam = FamsExtracted[Pos];
        FamsExtracted.erase(FamsExtracted.begin() + Pos);
        float Ind = this->FindBestPosToInsertFamForAllFac_Ind_DR(Sol, JobSeqInFam, JFDTime, JBDTime, CurFam,
                                                                 BestFac, BestPos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
        Sol[BestFac].insert(Sol[BestFac].begin() + BestPos, CurFam);
        RefreshJFDJBDTime_InFactory(Sol[BestFac], JobSeqInFam, JFDTime, JBDTime, BestPos, BestPos, 0, JobSeqInFam[CurFam].size() - 1);

        //cout << "�ڹ���" << BestFac << "��" << BestPos << "λ�ò���, ָ��Ϊ��" << Ind << endl;

        int FamPos = find(begin(Sol[BestFac]), end(Sol[BestFac]), CurFam) - begin(Sol[BestFac]);

        for (int i = 0; i < JobsExtracted[CurFam].size(); i++)
        {
            int BestJobPos = -1;
            float minFacInd = INT_MAX;
            //vector<int> TempFamSeqInFac = NewFamSeqInFac;
            for (int Pos = 0; Pos <= JobSeqInFam[CurFam].size(); Pos++)
            {
                float FacInd = GetIndForPerFacAfterInsertJob_DR(BestFac, Sol, Sol[BestFac], JobSeqInFam, JFDTime, JBDTime, CurFam, FamPos, JobsExtracted[CurFam][i], Pos,
                                                                nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);

                if (FacInd < minFacInd)
                {
                    minFacInd = FacInd;
                    BestPos = Pos;
                }

            }
            JobSeqInFam[CurFam].insert(JobSeqInFam[CurFam].begin() + BestPos, JobsExtracted[CurFam][i]);
            RefreshJFDJBDTime_InFactory(Sol[BestFac], JobSeqInFam, JFDTime, JBDTime, FamPos, FamPos, BestPos, BestPos);

            //cout << "�ڹ���" << BestFac << "��" << CurFam << "���" << BestPos << "λ�ò���, ָ��Ϊ��" << minFacInd << endl;
        }


    }
}

void CMOEA::UpdateArchiveGroupJobSet(int Muu)
{

    //��һ��
    //Normalize(CMOEAPopulation, m_NadirPointMS, m_NadirPointTEC, m_IdealPointMS, m_IdealPointTEC);
    //cout << "nadirpointMS��" << m_NadirPointMS << "\tnadirpointTEC��" << m_NadirPointTEC << endl;

    //��������ָ��
    //calc_convergence_ind(CMOEAPopulation);

    //---------------------------------------------------------------------------------------------------------------
    vector<Individual> Temp;
    Temp.clear();
    //---------------------------------------------------------------------------------------------------------------
    //jia
    //for (int i = CMOEAPopulation.size() -1; i > CMOEAPopulation.size() - 1 - m_RefSize; i--)
    //{
    //	Temp.push_back(CMOEAPopulation[i]);
    //}

    //ռ�Ź�ϵ
    Pareto_relation(CMOEAPopulation);

    //����֧��� de-emphasize dominated solutions
    for (int j = 0; j < CMOEAPopulation.size(); j++)
    {
        for (int i = 0; i < CMOEAPopulation.size(); i++)
        {
            if (CMOEAPopulation[i].flag == 0)
            {
                if (CMOEAPopulation[i].pareto_rel[j] == 1)
                {
                    CMOEAPopulation[i].flag = 999;
                }
            }
        }
    }


    for (int i = 0; i < CMOEAPopulation.size(); i++)
    {
        if (CMOEAPopulation[i].flag == 0)
            Temp.push_back(CMOEAPopulation[i]);
    }

    //������mu(Archive��Ⱥ����)
    int mu = Muu;
    if (Temp.size() < mu * m_RefSize)
    {
        while (true)
        {
            for (int i = 0; i < CMOEAPopulation.size(); i++)
            {
                if (CMOEAPopulation[i].flag == 999)
                {
                    Temp.push_back(CMOEAPopulation[i]);
                    if (Temp.size() == 10 * m_RefSize)
                        break;
                }
            }
            break;
        }
    }

    CMOEAPopulation.clear();
    for (int i = 0; i < Temp.size(); i++)
        CMOEAPopulation.push_back(Temp[i]);


    //��һ��
    Normalize(CMOEAPopulation, m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC);
    //cout << "nadirpointMS��" << m_NadirPointMS << "\tnadirpointTEC��" << m_NadirPointTEC << endl;

    //��������ָ��
    calc_convergence_ind(CMOEAPopulation);

    //�����ɢָ��
    calc_distribution_ind(CMOEAPopulation);

    int n, n1, n2, n3, nrank;
    int the_one;
    n = 0;   //Qs.size  Ԥѡ��ĸ��弯��
    n1 = 0;  //Q.size   �ο����ĸ�������
    n2 = 0;  //Qth.size ��distribution threshold�����ļ���
    n3 = 0;  //Qd.size  ��֧�����ļ���

    vector<Individual> TempCMOEAPopulation;
    TempCMOEAPopulation.clear();

    //ռ�Ź�ϵ
    Pareto_relation(CMOEAPopulation);

    n1 = CMOEAPopulation.size();

    while (n1 > 0)
    {
        for (int i = 0; i < CMOEAPopulation.size(); i++)
        {
            if (CMOEAPopulation[i].flag == 0)
            {
                the_one = i;
                break;
            }
        }

        CMOEAPopulation[the_one].flag = 1;
        TempCMOEAPopulation.push_back(CMOEAPopulation[the_one]);
        n++;


        //���������� de-emphasize neighbors
        for (int i = 0; i < CMOEAPopulation.size(); i++)
        {
            if (CMOEAPopulation[i].flag == 0)
            {
                if (CMOEAPopulation[i].distribution_ind[the_one] < thr_zeta)
                {
                    CMOEAPopulation[i].flag = 999;
                    n2++;
                }
            }
        }

        //����֧��� de-emphasize dominated solutions
        for (int i = 0; i < CMOEAPopulation.size(); i++)
        {
            if (CMOEAPopulation[i].flag == 0)
            {
                if (CMOEAPopulation[i].pareto_rel[the_one] == 1)
                {
                    CMOEAPopulation[i].flag = 999;
                    n3++;
                }
            }
        }

        //ʣ������ number of the rest
        n1 = 0;
        for (int i = 0; i < CMOEAPopulation.size(); i++)
        {
            if (CMOEAPopulation[i].flag == 0)
            {
                n1++;
            }
        }
    }

    //���� threshold
    float oldthr_zeta = thr_zeta;
    float ratio = n * 1.0 / (m_RefSize * 1.0);
    if (n3 < (mu - 1) * m_RefSize)
        thr_zeta = thr_zeta * exp((ratio - 1.0) / (2 * 1.0));
    else
        thr_zeta = oldthr_zeta;


    while (n < m_RefSize)
    {
        //ʣ������ number of the rest
        n1 = 0;
        for (int i = 0; i < CMOEAPopulation.size(); i++)
        {
            if (CMOEAPopulation[i].flag == 0)
            {
                n1++;
            }
        }

        if (n1 == 0)
        {
            for (int i = 0; i < CMOEAPopulation.size(); i++)
            {
                if (CMOEAPopulation[i].flag == 999)
                {
                    CMOEAPopulation[i].flag = 0;
                }
            }
        }

        //ѡ��һ��
        for (int i = 0; i < CMOEAPopulation.size(); i++)
        {
            if (CMOEAPopulation[i].flag == 0)
            {
                the_one = i;
                break;
            }
        }
        CMOEAPopulation[the_one].flag = 1;
        TempCMOEAPopulation.push_back(CMOEAPopulation[the_one]);
        n++;

        //���������� de-emphasize neighbors
        for (int i = 0; i < CMOEAPopulation.size(); i++)
        {
            if (CMOEAPopulation[i].flag == 0)
            {
                if (CMOEAPopulation[i].distribution_ind[the_one] < oldthr_zeta)
                {
                    CMOEAPopulation[i].flag = 999;
                    n2++;
                }
            }
        }

        //����֧��� de-emphasize dominated solutions
        for (int i = 0; i < CMOEAPopulation.size(); i++)
        {
            if (CMOEAPopulation[i].flag == 0)
            {
                if (CMOEAPopulation[i].pareto_rel[the_one] == 1)
                {
                    CMOEAPopulation[i].flag = 999;
                    n3++;
                }
            }
        }
    }

    vector<Individual> Temp1;
    Temp1.clear();
    // m_RefSize
    for (int i = 0; i < TempCMOEAPopulation.size(); i++)
    {
        Temp1.push_back(TempCMOEAPopulation[i]);
    }

    //----------------------------------------------------------------------------------------------------------------

    //���
    //cout << "����ǰ��1111111111" << endl;
    //for (int PS = 0; PS < CMOEAPopulation.size(); PS++)
    //{
    //	//m_SpeedMatrix = m_ArchiveSpeedVector[PS];
    //	m_SpeedMatrix = CMOEAPopulation[PS].m_SpeedVector;
    //	//�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
    //	for (int j = 0; j < m_Jobs; j++)
    //	{
    //		for (int i = 0; i < m_Machines; i++)
    //		{
    //			m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
    //			UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
    //		}
    //	}

    //	CheckSol(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].MS);
    //	CheckSolTEC(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].TEC);
    //}

    //����
    Speed_mutation(TempCMOEAPopulation, Temp1, CMOEAPopulation);

    //���
    //cout << "���ٺ�1111111111" << endl;
    //for (int PS = 0; PS < CMOEAPopulation.size(); PS++)
    //{
    //	//m_SpeedMatrix = m_ArchiveSpeedVector[PS];
    //	m_SpeedMatrix = CMOEAPopulation[PS].m_SpeedVector;
    //	//�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
    //	for (int j = 0; j < m_Jobs; j++)
    //	{
    //		for (int i = 0; i < m_Machines; i++)
    //		{
    //			m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
    //			UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
    //		}
    //	}
    //	CheckSol(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].MS);
    //	CheckSolTEC(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].TEC);
    //}

    //---------------------------------------------------------------------------------------------------------------


    for (int PS = 0; PS < m_RefSize; PS++)
    {
        // ���µ���������ֵ��
        m_RefJobSeqinFamArray[PS] = TempCMOEAPopulation[PS].m_JobSeqInFamArray;  //����������
        m_RefFacFamSeqArray[PS] = TempCMOEAPopulation[PS].m_FacFamSeqArray;   //����������
        m_RefSpanArray[PS] = TempCMOEAPopulation[PS].MS;  // ����깤ʱ��
        m_RefTECArray[PS] = TempCMOEAPopulation[PS].TEC;  //���ܺ�
        m_RefSpeedVector[PS] = TempCMOEAPopulation[PS].m_SpeedVector;

        //��������Ⱥ
        m_FacFamSeqArray[PS] = m_RefFacFamSeqArray[PS];

        //���¹�����Ⱥ
        m_JobSeqinFamArray[PS] = m_RefJobSeqinFamArray[PS];

        m_SpeedMatrix = m_RefSpeedVector[PS];
        //�õ������Ĵ���ʱ�� ����λʱ��ӹ��ܺ�
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }

        vector<int> FacSpan(m_Factories, 0);

        GetJFDTime_Forward(m_RefFacFamSeqArray[PS], m_RefJobSeqinFamArray[PS], m_TempJFDTimePop[PS], FacSpan);
        GetJBDTime_Backward(m_RefFacFamSeqArray[PS], m_RefJobSeqinFamArray[PS], m_TempJBDTimePop[PS], FacSpan);


        //���
        //CheckSol(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], m_ArchiveSpanArray[PS]);
        //CheckSolTEC(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], m_ArchiveTECArray[PS]);
    }
}

void CMOEA::SaveTECandDeMS(vector<Individual>& CMOEAPopulation)
{
    vector<Individual> tempCMOEAPopulation;
    tempCMOEAPopulation.clear();
    vector<Individual> tempCMOEAPopulation2;
    tempCMOEAPopulation2.clear();
    vector<Individual> tempCMOEAPopulation3;
    tempCMOEAPopulation3.clear();
    vector<Individual> tempCMOEAPopulation4;
    tempCMOEAPopulation4.clear();

    int orgSize = CMOEAPopulation.size();
    for (int PS = 0; PS < orgSize; PS++)
    {
        tempCMOEAPopulation.push_back(CMOEAPopulation[PS]);
        bool flag2 = false;

        //���� �� �ӹ��ܺ�PEC�ı�
        m_SpeedMatrix = tempCMOEAPopulation[PS].m_SpeedVector;
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }

        vector<int> FacSpan(m_Factories);//�飬�����깤ʱ��

        //���ܲ���1
        vector<vector<int>> DelayTime;
        DelayTime.clear();
        DelayTime.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
        GetDelayTime_Forward(tempCMOEAPopulation[PS].m_FacFamSeqArray, tempCMOEAPopulation[PS].m_JobSeqInFamArray, m_TempJFDTime, FacSpan, DelayTime);
        bool sign = false;

        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 1; i < m_Machines; i++)
            {
                if ((DelayTime[j][i] > 0) && (i < m_Machines - 1))
                {
                    m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[tempCMOEAPopulation[PS].m_SpeedVector[j][i]]));

                    m_TureJobOpertime[j][i + 1] = static_cast<int>(100 * (m_JobOperPTime[j][i + 1] / m_Speed[tempCMOEAPopulation[PS].m_SpeedVector[j][i + 1]]));

                    if ((m_TempJFDTime[j][i] < (m_TempJFDTime[j][i + 1] - m_TureJobOpertime[j][i + 1])))
                    {
                        int Speedlevel = tempCMOEAPopulation[PS].m_SpeedVector[j][i] - 1;

                        for (int level = Speedlevel; level > 0; level--)
                        {
                            if ((static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[level])) - m_TureJobOpertime[j][i]) < DelayTime[j][i])
                            {
                                if ((m_TempJFDTime[j][i] + (static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[level])) - m_TureJobOpertime[j][i])) < (m_TempJFDTime[j][i + 1] - m_TureJobOpertime[j][i + 1]))
                                {
                                    tempCMOEAPopulation[PS].m_SpeedVector[j][i] = level;
                                    sign = true;
                                }
                                else
                                    break;
                            }
                            else
                                break;
                        }
                    }

                }

                else if ((DelayTime[j][i] > 0) && (i == m_Machines - 1))
                {
                    int Speedlevel = tempCMOEAPopulation[PS].m_SpeedVector[j][i] - 1;
                    m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[tempCMOEAPopulation[PS].m_SpeedVector[j][i]]));
                    for (int level = Speedlevel; level > 0; level--)
                    {
                        if ((static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[level])) - m_TureJobOpertime[j][i]) < DelayTime[j][i])
                        {
                            tempCMOEAPopulation[PS].m_SpeedVector[j][i] = level;
                            sign = true;
                        }
                        else
                            break;
                    }
                }
            }//end machine
        }//end job

        if (sign == true)
        {
            int Makespan1 = 0;
            float TotalEC1 = 0.0;

            m_SpeedMatrix = tempCMOEAPopulation[PS].m_SpeedVector;

            //�õ������Ĵ���ʱ�� ����λ�ӹ�ʱ���ܺ�
            for (int j = 0; j < m_Jobs; j++)
            {
                for (int i = 0; i < m_Machines; i++)
                {
                    m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }


            Makespan1 = GetJFDTime_Forward(tempCMOEAPopulation[PS].m_FacFamSeqArray, tempCMOEAPopulation[PS].m_JobSeqInFamArray, m_TempJFDTime, FacSpan);

            TotalEC1 = GetTECForAllFacByJFD(tempCMOEAPopulation[PS].m_FacFamSeqArray, tempCMOEAPopulation[PS].m_JobSeqInFamArray, m_TempJFDTime);

            //���
            //CheckSol(tempCMOEAPopulation[PS].m_FacFamSeqArray, tempCMOEAPopulation[PS].m_JobSeqInFamArray, Makespan1);
            //CheckSolTEC(tempCMOEAPopulation[PS].m_FacFamSeqArray, tempCMOEAPopulation[PS].m_JobSeqInFamArray, TotalEC1);
            tempCMOEAPopulation[PS].MS = Makespan1;
            tempCMOEAPopulation[PS].TEC = TotalEC1;

            //�ж�
            if ((CMOEAPopulation[PS].MS >= tempCMOEAPopulation[PS].MS) && (CMOEAPopulation[PS].TEC > tempCMOEAPopulation[PS].TEC))
            {
                CMOEAPopulation.push_back(tempCMOEAPopulation[PS]);
                flag2 = true;
            }
        }


        //���ܲ���3����1�Ļ����ϣ��ڲ�����Makespan������£����ͷǹؼ�������v��������TEC
        vector<int> FacSpan3(m_Factories);
        if (flag2)
        {
            //float gapTEC = CMOEAPopulation[PS].TEC - tempCMOEAPopulation[PS].TEC;
            tempCMOEAPopulation3.push_back(tempCMOEAPopulation[PS]);
            vector<int> notcirfac;
            notcirfac.clear();
            for (int fac = 0; fac < m_Factories; fac++)
            {
                if (tempCMOEAPopulation3.back().MS > FacSpan[fac])
                {
                    notcirfac.push_back(fac);
                }
                FacSpan3[fac] = FacSpan[fac];
            }
            int diedai = 0;
            float tempTEC = 0.0;

            do {
                diedai++;
                if (diedai == 1)
                    tempTEC = tempCMOEAPopulation[PS].TEC;

                if (diedai == 2)
                    tempTEC = tempCMOEAPopulation3.back().TEC;
                //tempCMOEAPopulation[PS].TEC = tempCMOEAPopulation3[PS].TEC;
                if (diedai == 3)
                    break;

                bool judgeV = true;
                vector<vector<int>> tempSpeedMatrix;
                tempSpeedMatrix.clear();
                tempSpeedMatrix.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
                tempSpeedMatrix = tempCMOEAPopulation3.back().m_SpeedVector;

                //int countfac = 0;
                //int faccount = 0;
                //�Էǹؼ��������б���
                for (int f = 0; f < notcirfac.size(); f++)
                {
                    judgeV = true;
                    //int speedLevel;
                    for (int fam = 0; fam < tempCMOEAPopulation3.back().m_FacFamSeqArray[notcirfac[f]].size(); fam++)
                    {
                        int CurFam = tempCMOEAPopulation3.back().m_FacFamSeqArray[notcirfac[f]][fam];
                        for (int job = 0; job < tempCMOEAPopulation3.back().m_JobSeqInFamArray[CurFam].size(); job++)
                        {
                            int CurJob = tempCMOEAPopulation3.back().m_JobSeqInFamArray[CurFam][job];
                            for (int m = 0; m < m_Machines; m++)
                            {
                                if (tempCMOEAPopulation3.back().m_SpeedVector[CurJob][m] > 0)
                                {
                                    tempCMOEAPopulation3.back().m_SpeedVector[CurJob][m] = tempCMOEAPopulation3.back().m_SpeedVector[CurJob][m] - 1;
                                    judgeV = false;

                                    //��ʵ�Ĵ���ʱ��
                                    //m_TureJobOpertime[CurJob][m] = static_cast<int>(100 * (m_JobOperPTime[CurJob][m] / m_Speed[speedLevel]));
                                }

                            }
                        }
                    }
                    //�����ٶ��б仯
                    if (judgeV == false)
                    {
                        //faccount++;
                        m_SpeedMatrix = tempCMOEAPopulation3.back().m_SpeedVector;

                        //�õ������Ĵ���ʱ�� �� ��λ�ӹ��ܺ�PECϵ��
                        for (int j = 0; j < m_Jobs; j++)
                        {
                            for (int i = 0; i < m_Machines; i++)
                            {
                                m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
                                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                            }
                        }

                        //�õ�����f���ٺ��Span
                        int tempnotcirSpan = GetJFDTime_Forward_InFactory(tempCMOEAPopulation3.back().m_FacFamSeqArray[notcirfac[f]], tempCMOEAPopulation3.back().m_JobSeqInFamArray, m_TempJFDTime);//�õ�ÿ���������깤ʱ��
                        //����ı��ٶȺ�δ����Makspan������������
                        if (tempnotcirSpan > tempCMOEAPopulation3.back().MS)
                        {
                            tempCMOEAPopulation3.back().m_SpeedVector = tempSpeedMatrix;
                            //countfac++;
                            //judgeV = true;
                        }
                        else
                        {
                            tempSpeedMatrix = tempCMOEAPopulation3.back().m_SpeedVector;
                        }

                    }

                }

                /*if (judgeV && (countfac == m_Factories))
                    break;*/

                int Makespan3 = 0;
                float TotalEC3 = 0.0;

                m_SpeedMatrix = tempCMOEAPopulation3.back().m_SpeedVector;

                //�õ������Ĵ���ʱ�� �� ��λ�ӹ��ܺ�PECϵ��
                for (int j = 0; j < m_Jobs; j++)
                {
                    for (int i = 0; i < m_Machines; i++)
                    {
                        m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
                        UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                    }
                }

                Makespan3 = GetJFDTime_Forward(tempCMOEAPopulation3.back().m_FacFamSeqArray, tempCMOEAPopulation3.back().m_JobSeqInFamArray, m_TempJFDTime, FacSpan3);

                TotalEC3 = GetTECForAllFacByJFD(tempCMOEAPopulation3.back().m_FacFamSeqArray, tempCMOEAPopulation3.back().m_JobSeqInFamArray, m_TempJFDTime);

                //���
                //CheckSol(tempCMOEAPopulation3.back().m_FacFamSeqArray, tempCMOEAPopulation3.back().m_JobSeqInFamArray, Makespan3);
                //CheckSolTEC(tempCMOEAPopulation3.back().m_FacFamSeqArray, tempCMOEAPopulation3.back().m_JobSeqInFamArray, TotalEC3);
                tempCMOEAPopulation3.back().MS = Makespan3;
                tempCMOEAPopulation3.back().TEC = TotalEC3;

                //�ж�
                if ((CMOEAPopulation[PS].MS >= tempCMOEAPopulation3.back().MS) && (tempCMOEAPopulation[PS].TEC > tempCMOEAPopulation3.back().TEC))
                {
                    CMOEAPopulation.push_back(tempCMOEAPopulation3.back());
                }

            } while (tempTEC > tempCMOEAPopulation3.back().TEC);

        }

        //���ܲ���4����3�Ļ����ϣ��ڲ�����ԭTEC������£���߹ؼ�������v��������MS
        vector<int> FacSpan4(m_Factories);
        if (flag2)
        {
            //float gapTEC = CMOEAPopulation[PS].TEC - tempCMOEAPopulation[PS].TEC;
            tempCMOEAPopulation4.push_back(tempCMOEAPopulation3.back());
            vector<int> cirfac;
            cirfac.clear();
            for (int fac = 0; fac < m_Factories; fac++)
            {
                if (tempCMOEAPopulation4.back().MS == FacSpan[fac])
                {
                    cirfac.push_back(fac);
                }
                FacSpan4[fac] = FacSpan[fac];
            }
            int diedai = 0;
            int tempMS = 0;

            do {
                diedai++;
                if (diedai == 1)
                    tempMS = CMOEAPopulation[PS].MS;

                if (diedai == 2)
                {
                    tempMS = tempCMOEAPopulation4.back().MS;
                    //cirfac.clear();
                    vector<int> tempcirfac;
                    tempcirfac.clear();
                    //�жϹؼ������Ƿ�ı�
                    for (int fac = 0; fac < m_Factories; fac++)
                    {
                        if (tempCMOEAPopulation4.back().MS == FacSpan4[fac])
                        {
                            tempcirfac.push_back(fac);
                        }
                    }
                    bool judgecirfacifchange = true;
                    for (int i = 0; i < tempcirfac.size(); i++)
                    {
                        if (tempcirfac[i] != cirfac[i])
                        {
                            judgecirfacifchange = false;
                            break;
                        }
                    }
                    if (!judgecirfacifchange)
                    {
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
                tempSpeedMatrix = tempCMOEAPopulation4.back().m_SpeedVector;

                //�Էǹؼ��������б���
                for (int f = 0; f < cirfac.size(); f++)
                {
                    judgeV = true;
                    //int speedLevel;
                    for (int fam = 0; fam < tempCMOEAPopulation4.back().m_FacFamSeqArray[cirfac[f]].size(); fam++)
                    {
                        int CurFam = tempCMOEAPopulation4.back().m_FacFamSeqArray[cirfac[f]][fam];
                        for (int job = 0; job < tempCMOEAPopulation4.back().m_JobSeqInFamArray[CurFam].size(); job++)
                        {
                            int CurJob = tempCMOEAPopulation4.back().m_JobSeqInFamArray[CurFam][job];
                            for (int m = 0; m < m_Machines; m++)
                            {
                                if (tempCMOEAPopulation4.back().m_SpeedVector[CurJob][m] < 2)
                                {
                                    tempCMOEAPopulation4.back().m_SpeedVector[CurJob][m] = tempCMOEAPopulation4.back().m_SpeedVector[CurJob][m] + 1;
                                    judgeV = false;
                                }

                            }
                        }
                    }
                    //�����ٶ��б仯
                    if (judgeV == false)
                    {
                        //faccount++;
                        m_SpeedMatrix = tempCMOEAPopulation4.back().m_SpeedVector;

                        //�õ������Ĵ���ʱ�� �� ��λ�ӹ��ܺ�PECϵ��
                        for (int j = 0; j < m_Jobs; j++)
                        {
                            for (int i = 0; i < m_Machines; i++)
                            {
                                m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
                                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                            }
                        }

                        //�õ�����f���ٺ��TEC
                        GetJFDTime_Forward(tempCMOEAPopulation4.back().m_FacFamSeqArray, tempCMOEAPopulation4.back().m_JobSeqInFamArray, m_TempJFDTime, FacSpan4);
                        float tempTotalEC = GetTECForAllFacByJFD(tempCMOEAPopulation4.back().m_FacFamSeqArray, tempCMOEAPopulation4.back().m_JobSeqInFamArray, m_TempJFDTime);


                        //����ı��ٶȺ�δ����ԭTEC������������
                        if (tempTotalEC > CMOEAPopulation[PS].TEC)
                        {
                            tempCMOEAPopulation4.back().m_SpeedVector = tempSpeedMatrix;
                            break;
                        }

                    }

                }

                int Makespan4 = 0;
                float TotalEC4 = 0.0;

                m_SpeedMatrix = tempCMOEAPopulation4.back().m_SpeedVector;

                //�õ������Ĵ���ʱ�� �� ��λ�ӹ��ܺ�PECϵ��
                for (int j = 0; j < m_Jobs; j++)
                {
                    for (int i = 0; i < m_Machines; i++)
                    {
                        m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
                        UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                    }
                }

                Makespan4 = GetJFDTime_Forward(tempCMOEAPopulation4.back().m_FacFamSeqArray, tempCMOEAPopulation4.back().m_JobSeqInFamArray, m_TempJFDTime, FacSpan4);

                TotalEC4 = GetTECForAllFacByJFD(tempCMOEAPopulation4.back().m_FacFamSeqArray, tempCMOEAPopulation4.back().m_JobSeqInFamArray, m_TempJFDTime);

                //���
                //CheckSol(tempCMOEAPopulation4.back().m_FacFamSeqArray, tempCMOEAPopulation4.back().m_JobSeqInFamArray, Makespan4);
                //CheckSolTEC(tempCMOEAPopulation4.back().m_FacFamSeqArray, tempCMOEAPopulation4.back().m_JobSeqInFamArray, TotalEC4);
                tempCMOEAPopulation4.back().MS = Makespan4;
                tempCMOEAPopulation4.back().TEC = TotalEC4;

                //�ж�
                if ((CMOEAPopulation[PS].MS > tempCMOEAPopulation4.back().MS) && (CMOEAPopulation[PS].TEC >= tempCMOEAPopulation4.back().TEC))
                {
                    CMOEAPopulation.push_back(tempCMOEAPopulation4.back());
                }

            } while (tempMS > tempCMOEAPopulation4.back().MS);

        }
    }

}

void CMOEA::InitialPop(vector<vector<vector<int>>>& JFDTimePop, vector<vector<vector<int>>>& JBDTimePop)
{
    // cout<<"InitialPop start"<<endl;
    //��ʼ�������� Initialize Reference set and Set flags
    for (int PS = 0; PS < m_RefSize; PS++)
    {
        //�õ������Ĵ���ʱ�� ����λʱ��ӹ��ܺ�
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                //m_SpeedMatrix[j][i] = rand() % 3;
                m_SpeedMatrix[j][i] = 0;
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }


        vector<int> FamPrmu(m_Families), FacSpan(m_Factories);//�飬�����깤ʱ��
        vector<float> FacEC(m_Factories);
        vector<vector<int>> JobSeqInFam, FacFamSeq, newFacFamSeq(m_Factories); //������������У����ڹ��������У������ڹ���������
        if (PS == 0 || PS == 1) //��һ���ο�����LPT������������� Generate job sequence in each family and family sequence using LPT
        {
//            GetJobTotalPTime();  //�����ܼӹ�ʱ��
//            GetFamTotalPTime();  //���ܼӹ�ʱ��

            this->SortJobsInFam(0, JobSeqInFam); //LPT
            this->SortFam(0, FamPrmu); //LPT

        }

        else if (PS == 2 || PS == 3)
        {
            this->SortJobsInFam(1, JobSeqInFam); //SPT
            this->SortFam(1, FamPrmu); //SPT
        }

        else // ������ɹ������к������� Generate job sequence in each family and family sequence randomly
        {
            JobSeqInFam = this->m_JobsInEachFamily;
            for (int fam = 0; fam < JobSeqInFam.size(); fam++)
                random_shuffle(JobSeqInFam[fam].begin(), JobSeqInFam[fam].end());//�������ڹ���˳��
            for (int fam = 0; fam < FamPrmu.size(); fam++)
                FamPrmu[fam] = fam;
            random_shuffle(FamPrmu.begin(), FamPrmu.end());//������˳��
        }

        FacFamSeq.clear();
        FacFamSeq.resize(this->m_Factories);

        int CurFam = -1;
        int BestFac = -1, BestPos = -1;

        int Fam = 0;
        for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
        {
            CurFam = FamPrmu[Fam];
            FacFamSeq[Fac].insert(FacFamSeq[Fac].begin() + 0, CurFam);
            RefreshJFDJBDTime_InFactory(FacFamSeq[Fac], JobSeqInFam, JFDTimePop[PS], JBDTimePop[PS], 0, 0, 0, JobSeqInFam[CurFam].size() - 1);
            Fam++;
        }

        if (PS == 0 || PS == 2)
        {
            for (; Fam < FamPrmu.size(); Fam++)
            {
                CurFam = FamPrmu[Fam];
                BestFac = -1;
                BestPos = -1;
                this->FindBestPosToInsertFamForAllFac_Makespan(FacFamSeq, JobSeqInFam, JFDTimePop[PS], JBDTimePop[PS], CurFam, BestFac, BestPos);
                FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
                RefreshJFDJBDTime_InFactory(FacFamSeq[BestFac], JobSeqInFam, JFDTimePop[PS], JBDTimePop[PS], BestPos, BestPos, 0, JobSeqInFam[CurFam].size() - 1);

            }
        }
        else if (PS == 1 || PS == 3)
        {
            for (; Fam < FamPrmu.size(); Fam++)
            {
                CurFam = FamPrmu[Fam];
                BestFac = -1;
                BestPos = -1;
                this->FindBestPosToInsertFamForAllFac_TEC(FacFamSeq, JobSeqInFam, JFDTimePop[PS], JBDTimePop[PS], CurFam, BestFac, BestPos);
                FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
                RefreshJFDJBDTime_InFactory(FacFamSeq[BestFac], JobSeqInFam, JFDTimePop[PS], JBDTimePop[PS], BestPos, BestPos, 0, JobSeqInFam[CurFam].size() - 1);

            }
        }

        else
        {
            int r = rand() % 2;
            if (r == 0)
            {
                for (; Fam < FamPrmu.size(); Fam++)
                {
                    CurFam = FamPrmu[Fam];
                    BestFac = -1;
                    BestPos = -1;
                    this->FindBestPosToInsertFamForAllFac_Makespan(FacFamSeq, JobSeqInFam, JFDTimePop[PS], JBDTimePop[PS], CurFam, BestFac, BestPos);
                    FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
                    RefreshJFDJBDTime_InFactory(FacFamSeq[BestFac], JobSeqInFam, JFDTimePop[PS], JBDTimePop[PS], BestPos, BestPos, 0, JobSeqInFam[CurFam].size() - 1);

                }
            }
            else
            {
                for (; Fam < FamPrmu.size(); Fam++)
                {
                    CurFam = FamPrmu[Fam];
                    BestFac = -1;
                    BestPos = -1;
                    this->FindBestPosToInsertFamForAllFac_TEC(FacFamSeq, JobSeqInFam, JFDTimePop[PS], JBDTimePop[PS], CurFam, BestFac, BestPos);
                    FacFamSeq[BestFac].insert(FacFamSeq[BestFac].begin() + BestPos, CurFam);
                    RefreshJFDJBDTime_InFactory(FacFamSeq[BestFac], JobSeqInFam, JFDTimePop[PS], JBDTimePop[PS], BestPos, BestPos, 0, JobSeqInFam[CurFam].size() - 1);

                }
            }

        }


        int MS = 0; float TEC = 0;
        // cout << endl << "��Ⱥ�е�" << PS << "������" << endl;
        GetMSandTECForPerandToalFacByJFD(FacFamSeq, JobSeqInFam, JFDTimePop[PS], FacSpan, FacEC, MS, TEC);
        // cout << "Makespan��" << MS << "\t" << "TEC��" << TEC << endl;

        // ��ʼ������������ֵ�� Initialize a reference
        m_RefJobSeqinFamArray[PS] = JobSeqInFam;  //����������
        m_RefFacFamSeqArray[PS] = FacFamSeq;   //����������
        m_RefFacSpanArray[PS] = FacSpan;  //���й������깤ʱ��
        m_RefFacECArray[PS] = FacEC;  //���й������ܺ�
        m_RefSpanArray[PS] = MS;  // ����깤ʱ��
        m_RefTECArray[PS] = TEC;  //���ܺ�
        m_RefSpeedVector[PS] = m_SpeedMatrix;
        m_bFlag1[PS] = false;  //���1
        m_bFlag2[PS] = false;  //���2

        //���makespan��TEC
        this->CheckSol(m_RefFacFamSeqArray[PS], m_RefJobSeqinFamArray[PS], m_RefSpanArray[PS]);
        this->CheckSolTEC(m_RefFacFamSeqArray[PS], m_RefJobSeqinFamArray[PS], m_RefTECArray[PS]);

        //���²ο���
        CMOEAPopulation[PS].m_FacFamSeqArray = FacFamSeq;
        CMOEAPopulation[PS].m_JobSeqInFamArray = JobSeqInFam;
        CMOEAPopulation[PS].MS = MS;
        CMOEAPopulation[PS].TEC = TEC;
        CMOEAPopulation[PS].m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���
    }

    // ��ʼ������Ⱥ Initilize Familiy-sequence population, i.e., PS1
    for (int PS = 0; PS < m_PS1; PS++)
    {
        //���������Ľ⸳��������
        m_FacFamSeqArray[PS] = m_RefFacFamSeqArray[PS];
        m_SpanArray1[PS] = m_RefSpanArray[PS];
        m_TECArray1[PS] = m_RefTECArray[PS];
        m_Map1[PS] = PS;  //��ǵ������ж�Ӧ��jobpop���
    }

    // ��ʼ��������Ⱥ Initialize Job-Sequence population, i.e., PS2
    for (int PS = 0; PS < m_PS2; PS++)
    {
        //ǰAS����
        m_JobSeqinFamArray[PS] = m_RefJobSeqinFamArray[PS];
        m_SpanArray2[PS] = m_RefSpanArray[PS];
        m_TECArray2[PS] = m_RefTECArray[PS];
        m_Map2[PS] = PS;  //��Ƕ�Ӧ�ĵ�������fampop���
        m_nCriFacArray2[PS] = m_nRefCriFacArray[PS];
    }

    //��һ��
    Normalize(CMOEAPopulation, m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC);
    //cout << "nadirpointMS��" << m_NadirPointMS << "\tnadirpointTEC��" << m_NadirPointTEC << endl;

    //��������ָ��
    //calc_convergence_ind(CMOEAPopulation);

    //�����ɢָ��
    //calc_distribution_ind(CMOEAPopulation);
    // cout<<"InitialPop end"<<endl;
}

void CMOEA::EvolutionProcess(int mu)
{
   // cout<<"������ʼ"<<endl;
    vector<vector<int>> FacFamSeq(this->m_Factories);
    vector<vector<int>> JobSeqInFam;
    thr_zeta = 1.0;

//    //���
//    cout << "1111111111" << endl;
//    for (int PS = 0; PS < CMOEAPopulation.size(); PS++)
//    {
//    	//m_SpeedMatrix = m_ArchiveSpeedVector[PS];
//    	m_SpeedMatrix = CMOEAPopulation[PS].m_SpeedVector;
//    	//�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
//    	for (int j = 0; j < m_Jobs; j++)
//    	{
//    		for (int i = 0; i < m_Machines; i++)
//    		{
//    			m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
//    			UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
//    		}
//    	}
//    	CheckSol(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].MS);
//    	CheckSolTEC(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].TEC);
//    }

    long InitTime = Base::GetElapsedProcessTime();
    m_InitTime = InitTime;
    int count = 0;	//����
    //����
    while (::GetElapsedProcessTime() - InitTime < m_TimeLimit)
    {
        //Эͬ����--��
        for (int PS = 0; PS < m_PS1; PS++)
        {
            m_SpeedMatrix = m_RefSpeedVector[PS];
            //�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
            for (int j = 0; j < m_Jobs; j++)
            {
                for (int i = 0; i < m_Machines; i++)
                {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }

            vector<vector<int>> FacFamSeq = m_FacFamSeqArray[PS];
            int Map1 = rand() % m_RefSize;// �����ѡһ�������ߣ��������У��ӵ������� form a solution by randomly selecting a RefJobSeq;
            vector<int> FacSpan(m_Factories, 0);
            vector<float> FacTEC(m_Factories, 0);

            int ObjectMS = -1;
            float ObjectTEC = -1;

            //��makespan���Ĺ������������й�������������õ�λ�ã���������ָ����룩
            BasedindRandFamInFacTobestPos(FacFamSeq, m_RefJobSeqinFamArray[Map1], m_TempJFDTime, m_TempJBDTime,
                                          m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC, CMOEAPopulation, ObjectMS, ObjectTEC);

//            //���
//            cout << "after group method 1��1111111111" << endl;
//            for (int PS = 0; PS < CMOEAPopulation.size(); PS++)
//            {
//            	m_SpeedMatrix = CMOEAPopulation[PS].m_SpeedVector;
//            	//�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
//            	for (int j = 0; j < m_Jobs; j++)
//            	{
//            		for (int i = 0; i < m_Machines; i++)
//            		{
//            			m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
//            			UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
//            		}
//            	}
//            	CheckSol(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].MS);
//            	CheckSolTEC(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].TEC);
//            }

            m_SpeedMatrix = m_RefSpeedVector[PS];
            //�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
            for (int j = 0; j < m_Jobs; j++)
            {
            	for (int i = 0; i < m_Machines; i++)
            	{
            		m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
            		UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            	}
            }

            //�ӹؼ����������Ź�����������������飨��������ָ�꽻����
            BasedindSwapFam(FacFamSeq, m_RefJobSeqinFamArray[Map1], m_TempJFDTime, m_TempJBDTime,
                            m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC, CMOEAPopulation, ObjectMS, ObjectTEC);

//            //���
//            cout << "after group method 2��1111111111" << endl;
//            for (int PS = 0; PS < CMOEAPopulation.size(); PS++)
//            {
//            	//m_SpeedMatrix = m_ArchiveSpeedVector[PS];
//            	m_SpeedMatrix = CMOEAPopulation[PS].m_SpeedVector;
//            	//�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
//            	for (int j = 0; j < m_Jobs; j++)
//            	{
//            		for (int i = 0; i < m_Machines; i++)
//            		{
//            			m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
//            			UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
//            		}
//            	}
//            	CheckSol(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].MS);
//            	CheckSolTEC(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].TEC);
//            }

        }

//        //���
//        cout << "after group evolution��1111111111" << endl;
//        for (int PS = 0; PS < CMOEAPopulation.size(); PS++)
//        {
//        	//m_SpeedMatrix = m_ArchiveSpeedVector[PS];
//        	m_SpeedMatrix = CMOEAPopulation[PS].m_SpeedVector;
//        	//�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
//        	for (int j = 0; j < m_Jobs; j++)
//        	{
//        		for (int i = 0; i < m_Machines; i++)
//        		{
//        			m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
//        			UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
//        		}
//        	}
//        	CheckSol(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].MS);
//        	CheckSolTEC(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].TEC);
//        }

        //Эͬ����--����
        for (int PS = 0; PS < m_PS2; PS++)
        {
            m_SpeedMatrix = m_RefSpeedVector[PS];
            //�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
            for (int j = 0; j < m_Jobs; j++)
            {
                for (int i = 0; i < m_Machines; i++)
                {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }

            vector<vector<int>> JobSeqInFam = m_JobSeqinFamArray[PS];
            int Map2 = rand() % m_RefSize;// �����ѡһ�������ߣ������У� form a solution by randomly selecting a ReFacFamseq
            //int Map2 = PS;
            vector<int> FacSpan(m_Factories, 0);
            vector<float> FacTEC(m_Factories, 0);

            int ObjectMS = -1;
            float ObjectTEC = -1;

            //�Թؼ������е����ڹ��������λ��
            BasedindCirJobInFamTobestPos(m_RefFacFamSeqArray[Map2], JobSeqInFam, m_TempJFDTime, m_TempJBDTime,
                                         m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC, CMOEAPopulation, ObjectMS, ObjectTEC);

            //�Թؼ������е����ڹ����������
            BasedindSwapJob(m_RefFacFamSeqArray[Map2], JobSeqInFam, m_TempJFDTime, m_TempJBDTime,
                            m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC, CMOEAPopulation, ObjectMS, ObjectTEC);

        }

//        //���
//        cout << "after job evolution��1111111111" << endl;
//        for (int PS = 0; PS < CMOEAPopulation.size(); PS++)
//        {
//        	//m_SpeedMatrix = m_ArchiveSpeedVector[PS];
//        	m_SpeedMatrix = CMOEAPopulation[PS].m_SpeedVector;
//        	//�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
//        	for (int j = 0; j < m_Jobs; j++)
//        	{
//        		for (int i = 0; i < m_Machines; i++)
//        		{
//        			m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
//        			UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
//        		}
//        	}
//        	CheckSol(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].MS);
//        	CheckSolTEC(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].TEC);
//        }

        //Эͬ����--������
        for (int PS = 0; PS < m_RefSize; PS++)
        {

            m_SpeedMatrix = m_RefSpeedVector[PS];
            //�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
            for (int j = 0; j < m_Jobs; j++)
            {
                for (int i = 0; i < m_Machines; i++)
                {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }

            //����͹���ִ���ƻ��ع�����
            BasedindRefInsert(m_RefFacFamSeqArray[PS], m_RefJobSeqinFamArray[PS], m_TempJFDTimePop[PS], m_TempJBDTimePop[PS],
                              m_nadirpointMS, m_nadirpointTEC, m_idealpointMS, m_idealpointTEC, CMOEAPopulation, m_RefSpanArray[PS], m_RefTECArray[PS]);

            //����͹���ִ�н���
            /*BasedindRefSwap(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], m_TempJFDTimePop[PS], m_TempJBDTimePop[PS],
                m_NadirPointMS, m_NadirPointTEC, m_IdealPointMS, m_IdealPointTEC, CMOEAPopulation, m_ArchiveSpanArray[PS], m_ArchiveTECArray[PS]);*/

        }

//        //���
//        cout << " behind archiveGroupJobSet evolution��1111111111" << endl;
//        for (int PS = 0; PS < CMOEAPopulation.size(); PS++)
//        {
//        	m_SpeedMatrix = CMOEAPopulation[PS].m_SpeedVector;
//        	//�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
//        	for (int j = 0; j < m_Jobs; j++)
//        	{
//        		for (int i = 0; i < m_Machines; i++)
//        		{
//        			m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
//        			UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
//        		}
//        	}
//        	CheckSol(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].MS);
//        	CheckSolTEC(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].TEC);
//        }

        //ͨ���ο������µ��������飬���������ұ���
        UpdateArchiveGroupJobSet(mu);
        count++;

//        //���
//        cout <<"after 1 iter��1111111111" << endl;
//        for (int PS = 0; PS < CMOEAPopulation.size(); PS++)
//        {
//        	m_SpeedMatrix = CMOEAPopulation[PS].m_SpeedVector;
//        	//�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
//        	for (int j = 0; j < m_Jobs; j++)
//        	{
//        		for (int i = 0; i < m_Machines; i++)
//        		{
//        			m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
//        			UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
//        		}
//        	}
//        	CheckSol(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].MS);
//        	CheckSolTEC(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].TEC);
//        }

    }//END While

    //���
//    cout << "last" << endl;
//    for (int PS = 0; PS < CMOEAPopulation.size(); PS++)
//    {
//    	//m_SpeedMatrix = m_ArchiveSpeedVector[PS];
//    	m_SpeedMatrix = CMOEAPopulation[PS].m_SpeedVector;
//    	//�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
//    	for (int j = 0; j < m_Jobs; j++)
//    	{
//    		for (int i = 0; i < m_Machines; i++)
//    		{
//    			m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
//    			UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
//    		}
//    	}
//    	CheckSol(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].MS);
//    	CheckSolTEC(CMOEAPopulation[PS].m_FacFamSeqArray, CMOEAPopulation[PS].m_JobSeqInFamArray, CMOEAPopulation[PS].TEC);
//    }

    //���ܺͽ���MS
    SaveTECandDeMS(CMOEAPopulation);

//    cout << "����������" << count << endl;
//    cout << "zeta��" << thr_zeta << endl;
  //  cout<<"��������"<<endl;
}

int CMOEA::RunEvolution(int CPUTime, vector<vector<Individual>>& CMOEAFinalAfterRepParetoSet, int AN, int mu)
{
    ReadInstanceFileNameList("..\\Benchmark\\");
    int FamD = 6;
    int NoImproveNumber = 40;
    int BlockLengthChoice = 1;
    int Instances = 405;

    int Reps = 10; //���ÿ�������ظ����еĴ���

    vector<string> NameArray;
    vector<int> FacArray, FamArray, JobArray, MacArray, RunTimeArray, SetupTypeArray;
    vector<int> MakeSpanArray, OrgMakespanArray;

    string FileDirectory = "..\\Result\\";
    ostringstream str;
    str << "CMOEA_" << CPUTime << "_experiment" << ".txt"; //��ͬ���㷨
    ofstream ofile;
    ofile.open(FileDirectory + str.str());

    for (int ins = 0; ins < Instances; ins++)
    {

        srand((unsigned int)time(NULL));
        this->ReadInstance(ins);
        this->GetJobTotalPTime();
        this->GetFamTotalPTime();
        this->GetFamAvgSetupTime();
        this->GetFamTotalPTimeOnLastMachine();
        this->GetFamTotalPTimeOnFirstMachine();

        //��
        vector<Individual> FinalCMOEAParetoSet;
        vector<Individual> AfterRepParetoSet;

        vector<Individual> TempAfterRepParetoSet;
        TempAfterRepParetoSet.clear();

        for (int r = 0; r < Reps; r++)
        {
            long TimeLimit = CPUTime * m_Machines * m_Families; //original: 20 * m_Jobs * m_Machines

            this->SetParameters(AN, NoImproveNumber, FamD, TimeLimit, m_AllJobTotalPTime);

            m_JFDTime.clear();
            m_JFDTime.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
            m_TempJFDTime.clear();
            m_TempJFDTime.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
            m_BestJFDTime.clear();
            m_BestJFDTime.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
            m_tempJFDTime1.clear();
            m_tempJFDTime1.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
            m_JFDTimePop.clear();
            m_JFDTimePop.resize(this->m_Popsize, vector<vector<int>>(this->m_Jobs, vector<int>(this->m_Machines, 0)));
            m_TempJFDTimePop.clear();
            m_TempJFDTimePop.resize(this->m_Popsize, vector<vector<int>>(this->m_Jobs, vector<int>(this->m_Machines, 0)));
            m_BestJFDTimePop.clear();
            m_BestJFDTimePop.resize(this->m_Popsize, vector<vector<int>>(this->m_Jobs, vector<int>(this->m_Machines, 0)));
            m_tempJFDTime1Pop.clear();
            m_tempJFDTime1Pop.resize(this->m_Popsize, vector<vector<int>>(this->m_Jobs, vector<int>(this->m_Machines, 0)));

            m_JBDTime.clear();
            m_JBDTime.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
            m_TempJBDTime.clear();
            m_TempJBDTime.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
            m_BestJBDTime.clear();
            m_BestJBDTime.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
            m_tempJBDTime1.clear();
            m_tempJBDTime1.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
            m_JBDTimePop.clear();
            m_JBDTimePop.resize(this->m_Popsize, vector<vector<int>>(this->m_Jobs, vector<int>(this->m_Machines, 0)));
            m_TempJBDTimePop.clear();
            m_TempJBDTimePop.resize(this->m_Popsize, vector<vector<int>>(this->m_Jobs, vector<int>(this->m_Machines, 0)));
            m_BestJBDTimePop.clear();
            m_BestJBDTimePop.resize(this->m_Popsize, vector<vector<int>>(this->m_Jobs, vector<int>(this->m_Machines, 0)));
            m_tempJBDTime1Pop.clear();
            m_tempJBDTime1Pop.resize(this->m_Popsize, vector<vector<int>>(this->m_Jobs, vector<int>(this->m_Machines, 0)));

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

            long StartTime_IG = ::GetElapsedProcessTime();

            //��ʼ����Ⱥ
            InitialPop(m_TempJFDTimePop, m_TempJBDTimePop);

            //��������
//            EvolutionProcess(mu);

            //��֧���
            CMOEAParetoSet.clear();

            Pareto_relation(CMOEAPopulation);

            //����֧��� de-emphasize dominated solutions
            for (int j = 0; j < CMOEAPopulation.size(); j++)
            {
                for (int i = 0; i < CMOEAPopulation.size(); i++)
                {
                    if (CMOEAPopulation[i].flag == 0)
                    {
                        if (CMOEAPopulation[i].pareto_rel[j] == 1)
                        {
                            CMOEAPopulation[i].flag = 999;
                        }
                    }
                }
            }

            for (int i = 0; i < CMOEAPopulation.size(); i++)
            {
                if (CMOEAPopulation[i].flag == 0)
                    CMOEAParetoSet.push_back(CMOEAPopulation[i]);
            }

            //ȥ���ظ�
            FinalCMOEAParetoSet.clear();
            for (int i = 0; i < CMOEAParetoSet.size(); i++)
            {
                bool fg = true;
                for (int j = 0; j < FinalCMOEAParetoSet.size(); j++)
                {
                    if ((FinalCMOEAParetoSet[j].MS == CMOEAParetoSet[i].MS) && (FinalCMOEAParetoSet[j].TEC == CMOEAParetoSet[i].TEC))
                    {
                        fg = false;
                        break;
                    }
                }
                if (fg)
                {
                    FinalCMOEAParetoSet.push_back(CMOEAParetoSet[i]);
                }

            }

            long EndTime_IG = ::GetElapsedProcessTime();

            for (int PS = 0; PS < FinalCMOEAParetoSet.size(); PS++)
            {
                TempAfterRepParetoSet.push_back(FinalCMOEAParetoSet[PS]);
            }
        }//end rep

        //��֧���
        AfterRepParetoSet.clear();

        Pareto_relation(TempAfterRepParetoSet);

        //����֧��� de-emphasize dominated solutions
        for (int j = 0; j < TempAfterRepParetoSet.size(); j++)
        {
            for (int i = 0; i < TempAfterRepParetoSet.size(); i++)
            {
                if (TempAfterRepParetoSet[i].flag == 0)
                {
                    if (TempAfterRepParetoSet[i].pareto_rel[j] == 1)
                    {
                        TempAfterRepParetoSet[i].flag = 999;
                    }
                }
            }
        }

        for (int i = 0; i < TempAfterRepParetoSet.size(); i++)
        {
            if (TempAfterRepParetoSet[i].flag == 0)
                AfterRepParetoSet.push_back(TempAfterRepParetoSet[i]);
        }

        //ȥ���ظ�
        vector<Individual> FinalAfterRepParetoSet;
        FinalAfterRepParetoSet.clear();
        for (int i = 0; i < AfterRepParetoSet.size(); i++)
        {
            bool fg = true;
            for (int j = 0; j < FinalAfterRepParetoSet.size(); j++)
            {
                if ((FinalAfterRepParetoSet[j].MS == AfterRepParetoSet[i].MS) && (FinalAfterRepParetoSet[j].TEC == AfterRepParetoSet[i].TEC))
                {
                    fg = false;
                    break;
                }
            }
            if (fg)
            {
                FinalAfterRepParetoSet.push_back(AfterRepParetoSet[i]);
            }
        }

        CMOEAFinalAfterRepParetoSet[ins].clear();
        cout << ins + 1 << "\t" << "Factories:" << this->m_Factories<< "\t" << "Machines :" << this->m_Machines<< "\t" << "Families:" << this->m_Families<< "\t" << "Jobs :" << this->m_Jobs << "\t";

        for (int PS = 0; PS < FinalAfterRepParetoSet.size(); PS++)
        {
            // ��� MS �� TEC�����нⶼ����ͬһ��
            cout << "MS:" << FinalAfterRepParetoSet[PS].MS
                 << " " << "TEC:" << FinalAfterRepParetoSet[PS].TEC << "\t";
            ofile << FinalAfterRepParetoSet[PS].MS << "  "<< FinalAfterRepParetoSet[PS].TEC << ","<< "\t";
            CMOEAFinalAfterRepParetoSet[ins].push_back(FinalAfterRepParetoSet[PS]);
        }
        ofile << endl;  // ÿ��ʵ�������������
        cout << endl;   // ����̨������У���������ʵ������ʱ����

    }//end ins
    ofile.close();
    return 0;
}

