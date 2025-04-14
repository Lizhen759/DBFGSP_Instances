#include "MOSA.h"

MOSA::MOSA()
{
}

MOSA::~MOSA()
{
}

void MOSA::SetParameters(int RefSize, int PS1, int PS2, int AgeLimit, int DesLen, long TimeLimit) //PS1>RefSize; PS2>RefSize;
{
    //monituihuo
    T0 = 1000.00;
    TF = 10.00;
    iteration_number = m_Families;
    cool_ratio = 0.96;

    //-------------------------
    m_RefSize = RefSize;
    m_PS1 = PS1;
    m_PS2 = PS2;
    m_AgeLimit = AgeLimit;
    m_DesLen = DesLen;
    m_TimeLimit = TimeLimit;

    MOSAPopulation.clear();
    //MOSAPopulation.resize(m_RefSize);

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


    m_RefFacFamSeqArray.clear();//������������
    m_RefFacFamSeqArray.resize(m_RefSize);

    m_RefJobSeqinFamArray.clear();//�������鹤��
    m_RefJobSeqinFamArray.resize(m_RefSize);

    m_bFlag1.clear();//
    m_bFlag1.resize(m_RefSize);

    m_bFlag2.clear();//
    m_bFlag2.resize(m_RefSize);

    //��
    m_TECArray1.clear(); //�ܺ�
    m_TECArray1.resize(m_PS1);
    m_SpanArray1.clear();
    m_SpanArray1.resize(m_PS1);
    m_Map1.clear();
    m_Map1.resize(m_PS1);
    m_Age1.clear();
    m_Age1.resize(m_PS1);
    m_FacFamSeqArray.clear();//������
    m_FacFamSeqArray.resize(m_PS1);

    //����
    m_TECArray2.clear(); //�ܺ�
    m_TECArray2.resize(m_PS1);
    m_SpanArray2.clear();
    m_SpanArray2.resize(m_PS2);
    m_Map2.clear();
    m_Map2.resize(m_PS2);
    m_Age2.clear();
    m_Age2.resize(m_PS2);
    m_JobSeqinFamArray.clear();
    m_JobSeqinFamArray.resize(m_PS2);
}

void MOSA::InitPop()
{
    //ģ���˻�
    //�ٶȳ�ʼ��
    for (int j = 0; j < m_Jobs; j++)
    {
        for (int i = 0; i < m_Machines; i++)
        {
            m_SpeedMatrix[j][i] = rand() % 3;
            m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
            UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
        }
    }

    //���г�ʼ��
    vector<int> FamPrmu(m_Families), FacSpan(m_Factories);//�飬�����깤ʱ��
    vector<float> FacEC(m_Factories);

    vector<vector<int>> JobSeqinFam, FacFamSeq; //������������У����ڹ��������У������ڹ���������

    JobSeqinFam = m_JobsInEachFamily;
    for (int fam = 0; fam < JobSeqinFam.size(); fam++)
        random_shuffle(JobSeqinFam[fam].begin(), JobSeqinFam[fam].end());//�������ڹ���˳��
    for (int fam = 0; fam < FamPrmu.size(); fam++)
        FamPrmu[fam] = fam;
    random_shuffle(FamPrmu.begin(), FamPrmu.end());//������˳��

    // ���������ɵ�����䵽��������һ��������� Generte a solution using assignment rule--include
    GetSol_Include(FamPrmu, JobSeqinFam, FacFamSeq, FacSpan);

    int MS = 0; float TEC = 0;
    GetMSandTECForPerandToalFac(FacFamSeq, JobSeqinFam, FacSpan, FacEC, MS, TEC);
   // cout << "Makespan��" << MS << "\t" << "TEC��" << TEC << endl;

    //����
    MOSASol.m_FacFamSeqArray = FacFamSeq;
    MOSASol.m_JobSeqInFamArray = JobSeqinFam;
    MOSASol.MS = MS;
    MOSASol.TEC = TEC;
    MOSASol.m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���

    Individual tempIndi;
    tempIndi.m_FacFamSeqArray = FacFamSeq;
    tempIndi.m_JobSeqInFamArray = JobSeqinFam;
    tempIndi.MS = MS;
    tempIndi.TEC = TEC;
    tempIndi.m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���
    MOSAPopulation.push_back(tempIndi);

    //���
    this->CheckSol(MOSASol.m_FacFamSeqArray, MOSASol.m_JobSeqInFamArray, MOSASol.MS);
    this->CheckSolTEC(MOSASol.m_FacFamSeqArray, MOSASol.m_JobSeqInFamArray, MOSASol.TEC);

}

int MOSA::Destruction_Construction(int Len, vector<vector<int>>& FacFamSeq, vector<vector<int>> JobSeqinFam, vector<int>& FacSpan)
{
    // �ҵ��ؼ����� Find critical factories
    int mSpan = *max_element(FacSpan.begin(), FacSpan.end());
    vector<int> CriFac;
    for (int Fac = 0; Fac < m_Factories; Fac++)
    {
        if (mSpan == FacSpan[Fac])
            CriFac.push_back(Fac);
    }

    // Extract Len/2 families from cirtical factory
    vector<int> ExtractFamSeq;
    do {
        int Fac;
        if (ExtractFamSeq.size() < Len / 2) //�ӹؼ�������ѡd/2����
            Fac = CriFac[rand() % CriFac.size()];
        else  //������������ѡʣ�µ�d/2����
            Fac = rand() % m_Factories;
        if (FacFamSeq[Fac].size())
        {
            int Pos = rand() % FacFamSeq[Fac].size();
            ExtractFamSeq.push_back(FacFamSeq[Fac][Pos]);  //����ѡ��������ӵ�ExtractFamSeq
            FacFamSeq[Fac].erase(FacFamSeq[Fac].begin() + Pos);  //��ԭ������ɾ��posλ�õ���
        }
    } while (ExtractFamSeq.size() < Len);

    GetSpan(FacFamSeq, JobSeqinFam, FacSpan);
    // ����ѡ������d������뵽��õ�λ�� Insert the extracted Families into the best Positions
    return NEHFam(ExtractFamSeq, JobSeqinFam, FacFamSeq, FacSpan);
}

void MOSA::NGP1(vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqinFam, vector<vector<int>>& m_SpeedMatrix)
{
    //srand((unsigned int)time(NULL));

    int randomFac = -1;
    do {
        randomFac = rand() % FacFamSeq.size();
    } while (FacFamSeq[randomFac].size() < (m_Families / m_Factories));

    int randomPos = -1;

    do {
        randomPos = rand() % FacFamSeq[randomFac].size();

    } while ((randomPos == 0) || (randomPos == FacFamSeq[randomFac].size() - 1));


    vector<int> tempSeq1((randomPos + 1), 0);

    vector<int> tempSeq2((FacFamSeq[randomFac].size() - randomPos - 1), 0);

    for (int i = 0; i < randomPos + 1; i++)
    {
        tempSeq1[i] = FacFamSeq[randomFac][i];
    }

    int num = 0;
    for (int i = randomPos + 1; i < FacFamSeq[randomFac].size(); i++)
    {
        tempSeq2[num] = FacFamSeq[randomFac][i];
        num++;
    }

    //�ı�ԭ����
    for (int i = 0; i < (FacFamSeq[randomFac].size()- randomPos - 1); i++)
    {
        FacFamSeq[randomFac][i] = tempSeq2[i];
    }

    num = 0;
    for (int i = (FacFamSeq[randomFac].size() - randomPos - 1); i < FacFamSeq[randomFac].size(); i++)
    {
        FacFamSeq[randomFac][i] = tempSeq1[num];
        num++;
    }

    if (tempSeq1.size() < tempSeq2.size())
    {
        for (int i = 0; i < tempSeq1.size(); i++)
        {
            int CurFam = tempSeq1[i];
            if (JobSeqinFam[CurFam].size() > 2)
            {
                int randomJobPos = rand() % JobSeqinFam[CurFam].size();
                int CurJob = JobSeqinFam[CurFam][randomJobPos];
                JobSeqinFam[CurFam].erase(JobSeqinFam[CurFam].begin() + randomJobPos);

                randomJobPos = rand() % JobSeqinFam[CurFam].size();
                JobSeqinFam[CurFam].insert(JobSeqinFam[CurFam].begin() + randomJobPos, CurJob);

                //����
                for (int i = 0; i < m_Machines; i++)
                {
                    if (m_SpeedMatrix[CurJob][i] == 0)
                    {
                        m_SpeedMatrix[CurJob][i]++;
                    }
                    else if (m_SpeedMatrix[CurJob][i] == 1)
                    {
                        double p1 = rand() / (double)RAND_MAX;
                        if (p1 > 0.5)
                            m_SpeedMatrix[CurJob][i]++;
                        else
                            m_SpeedMatrix[CurJob][i]--;
                    }
                    else
                        m_SpeedMatrix[CurJob][i]--;
                }
            }

        }
    }
    else
    {
        for (int i = 0; i < tempSeq2.size(); i++)
        {
            int CurFam = tempSeq2[i];
            if (JobSeqinFam[CurFam].size() > 2)
            {
                int randomJobPos = rand() % JobSeqinFam[CurFam].size();
                int CurJob = JobSeqinFam[CurFam][randomJobPos];
                JobSeqinFam[CurFam].erase(JobSeqinFam[CurFam].begin() + randomJobPos);

                randomJobPos = rand() % JobSeqinFam[CurFam].size();
                JobSeqinFam[CurFam].insert(JobSeqinFam[CurFam].begin() + randomJobPos, CurJob);

                //����
                for (int i = 0; i < m_Machines; i++)
                {
                    if (m_SpeedMatrix[CurJob][i] == 0)
                    {
                        m_SpeedMatrix[CurJob][i]++;
                    }
                    else if (m_SpeedMatrix[CurJob][i] == 1)
                    {
                        double p1 = rand() / (double)RAND_MAX;
                        if (p1 > 0.5)
                            m_SpeedMatrix[CurJob][i]++;
                        else
                            m_SpeedMatrix[CurJob][i]--;
                    }
                    else
                        m_SpeedMatrix[CurJob][i]--;
                }
            }

        }
    }

}

void MOSA::NGP2(vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqinFam, vector<vector<int>>& m_SpeedMatrix)
{
    //srand((unsigned int)time(NULL));
    //������
    int randomFac1 = -1;
    do {
        randomFac1 = rand() % FacFamSeq.size();
    } while (FacFamSeq[randomFac1].size() < (m_Families/m_Factories));

    int randomPos1 = -1;
    int CurFam1 = -1;
    int A = 0;
    do {
        randomPos1 = rand() % FacFamSeq[randomFac1].size();
        CurFam1 = FacFamSeq[randomFac1][randomPos1];
        A++;
        if (A == 1000)
            return;
    } while (JobSeqinFam[CurFam1].size() < 3);

    //FacFamSeq[randomFac1].erase(FacFamSeq[randomFac1].begin() + randomPos1);

    int randomFac2 = -1;
    do {
        randomFac2 = rand() % FacFamSeq.size();
    } while (FacFamSeq[randomFac2].size() < (m_Families / m_Factories));

    int randomPos2 = -1;
    int CurFam2 = -1;
    int B = 0;
    do {
        randomPos2 = rand() % FacFamSeq[randomFac2].size();
        CurFam2 = FacFamSeq[randomFac2][randomPos2];
        B++;
        if (B == 1000)
            return;
    } while ((JobSeqinFam[CurFam2].size() < 2) || (CurFam2 == CurFam1));

    //FacFamSeq[randomFac2].erase(FacFamSeq[randomFac2].begin() + randomPos2);

    FacFamSeq[randomFac1][randomPos1] = CurFam2;
    FacFamSeq[randomFac2][randomPos2] = CurFam1;

    //������1����
    int randomJobPos1 = rand() % JobSeqinFam[CurFam1].size();
    int CurJob1 = JobSeqinFam[CurFam1][randomJobPos1];

    int randomJobPos2 = -1;
    int CurJob2 = -1;
    do {
        randomJobPos2 = rand() % JobSeqinFam[CurFam1].size();
        CurJob2 = JobSeqinFam[CurFam1][randomJobPos2];
    } while (CurJob1 == CurJob2);

    JobSeqinFam[CurFam1][randomJobPos1] = CurJob2;
    JobSeqinFam[CurFam1][randomJobPos2] = CurJob1;

    //����
    for (int i = 0; i < m_Machines; i++)
    {
        if (m_SpeedMatrix[CurJob1][i] == 0)
        {
            m_SpeedMatrix[CurJob1][i]++;
        }
        else if (m_SpeedMatrix[CurJob1][i] == 1)
        {
            double p1 = rand() / (double)RAND_MAX;
            if (p1 > 0.5)
                m_SpeedMatrix[CurJob1][i]++;
            else
                m_SpeedMatrix[CurJob1][i]--;
        }
        else
            m_SpeedMatrix[CurJob1][i]--;
    }

    //����
    for (int i = 0; i < m_Machines; i++)
    {
        if (m_SpeedMatrix[CurJob2][i] == 0)
        {
            m_SpeedMatrix[CurJob2][i]++;
        }
        else if (m_SpeedMatrix[CurJob2][i] == 1)
        {
            double p1 = rand() / (double)RAND_MAX;
            if (p1 > 0.5)
                m_SpeedMatrix[CurJob2][i]++;
            else
                m_SpeedMatrix[CurJob2][i]--;
        }
        else
            m_SpeedMatrix[CurJob2][i]--;
    }

    //������2����
    randomJobPos1 = rand() % JobSeqinFam[CurFam2].size();
    CurJob1 = JobSeqinFam[CurFam2][randomJobPos1];

    randomJobPos2 = -1;
    CurJob2 = -1;
    do {
        randomJobPos2 = rand() % JobSeqinFam[CurFam2].size();
        CurJob2 = JobSeqinFam[CurFam2][randomJobPos2];
    } while (CurJob1 == CurJob2);

    JobSeqinFam[CurFam2][randomJobPos1] = CurJob2;
    JobSeqinFam[CurFam2][randomJobPos2] = CurJob1;

    //����
    for (int i = 0; i < m_Machines; i++)
    {
        if (m_SpeedMatrix[CurJob1][i] == 0)
        {
            m_SpeedMatrix[CurJob1][i]++;
        }
        else if (m_SpeedMatrix[CurJob1][i] == 1)
        {
            double p1 = rand() / (double)RAND_MAX;
            if (p1 > 0.5)
                m_SpeedMatrix[CurJob1][i]++;
            else
                m_SpeedMatrix[CurJob1][i]--;
        }
        else
            m_SpeedMatrix[CurJob1][i]--;
    }

    //����
    for (int i = 0; i < m_Machines; i++)
    {
        if (m_SpeedMatrix[CurJob2][i] == 0)
        {
            m_SpeedMatrix[CurJob2][i]++;
        }
        else if (m_SpeedMatrix[CurJob2][i] == 1)
        {
            double p1 = rand() / (double)RAND_MAX;
            if (p1 > 0.5)
                m_SpeedMatrix[CurJob2][i]++;
            else
                m_SpeedMatrix[CurJob2][i]--;
        }
        else
            m_SpeedMatrix[CurJob2][i]--;
    }

}

void MOSA::NGP3(vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqinFam, vector<vector<int>>& m_SpeedMatrix)
{
    //srand((unsigned int)time(NULL));
    //������
    int randomFac = -1;
    do {
        randomFac = rand() % FacFamSeq.size();
    } while (FacFamSeq[randomFac].size() < (m_Families / m_Factories));

    int randomPos = -1;
    int CurFam = -1;
    int A = 0;
    do {
        //srand((unsigned int)time(NULL));
        randomPos = (rand() % FacFamSeq[randomFac].size());
        CurFam = FacFamSeq[randomFac][randomPos];
        A++;
        if (A == 1000)
            return;
    } while (JobSeqinFam[CurFam].size() < 3);

    FacFamSeq[randomFac].erase(FacFamSeq[randomFac].begin() + randomPos);

    randomFac = rand() % FacFamSeq.size();
    randomPos = rand() % FacFamSeq[randomFac].size();
    FacFamSeq[randomFac].insert(FacFamSeq[randomFac].begin() + randomPos, CurFam);

    //���빤��
    int randomJobPos = rand() % JobSeqinFam[CurFam].size();
    int CurJob = JobSeqinFam[CurFam][randomJobPos];
    JobSeqinFam[CurFam].erase(JobSeqinFam[CurFam].begin() + randomJobPos);

    randomJobPos = rand() % JobSeqinFam[CurFam].size();
    JobSeqinFam[CurFam].insert(JobSeqinFam[CurFam].begin() + randomJobPos, CurJob);

    //����
    for (int i = 0; i < m_Machines; i++)
    {
        if (m_SpeedMatrix[CurJob][i] == 0)
        {
            m_SpeedMatrix[CurJob][i]++;
        }
        else if (m_SpeedMatrix[CurJob][i] == 1)
        {
            double p1 = rand() / (double)RAND_MAX;
            if (p1 > 0.5)
                m_SpeedMatrix[CurJob][i]++;
            else
                m_SpeedMatrix[CurJob][i]--;
        }
        else
            m_SpeedMatrix[CurJob][i]--;
    }

}

void MOSA::LS_strategy(vector<vector<int>>& FacFamSeq, vector<vector<int>>& JobSeqinFam, vector<vector<int>>& m_SpeedMatrix)
{
    srand((unsigned int)time(NULL));
    int orgMS = 0;
    float orgTEC = 0.0;
    vector<int> orgFacSpan(m_Factories, 0);
    vector<float> orgFacEC(m_Factories, 0);
    //�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
    for (int j = 0; j < m_Jobs; j++)
    {
        for (int i = 0; i < m_Machines; i++)
        {
            m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
            UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
        }
    }
    GetMSandTECForPerandToalFac(FacFamSeq, JobSeqinFam, orgFacSpan, orgFacEC, orgMS, orgTEC);

    //�ҵ�ĳ�������ĵ�һ����
    int randomFac = -1;
    do {
        randomFac = rand() % FacFamSeq.size();
    } while (FacFamSeq[randomFac].size() < 2);

    int FirstFam = FacFamSeq[randomFac][0];

    if (JobSeqinFam[FirstFam].size() > 2)
    {
        //�������ڹ���
        int randomJobPos1 = rand() % JobSeqinFam[FirstFam].size();
        int CurJob1 = JobSeqinFam[FirstFam][randomJobPos1];

        int randomJobPos2 = -1;
        int CurJob2 = -1;
        do {
            randomJobPos2 = rand() % JobSeqinFam[FirstFam].size();
            CurJob2 = JobSeqinFam[FirstFam][randomJobPos2];
        } while (CurJob1 == CurJob2);

        JobSeqinFam[FirstFam][randomJobPos1] = CurJob2;
        JobSeqinFam[FirstFam][randomJobPos2] = CurJob1;

        //�����������ٶ�
        for (int i = 0; i < m_Machines; i++)
        {
            int tempSpeed;
            tempSpeed = m_SpeedMatrix[CurJob1][i];
            m_SpeedMatrix[CurJob1][i] = m_SpeedMatrix[CurJob2][i];
            m_SpeedMatrix[CurJob2][i] = tempSpeed;
        }

        int newMS = 0;
        float newTEC = 0.0;
        vector<int> newFacSpan(m_Factories, 0);
        vector<float> newFacEC(m_Factories, 0);
        //�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }
        GetMSandTECForPerandToalFac(FacFamSeq, JobSeqinFam, newFacSpan, newFacEC, newMS, newTEC);

        if ((newMS < orgMS) || (newTEC < orgTEC))
        {
            orgMS = newMS;
            orgTEC = newTEC;
            vector<vector<int>> tempFacFamSeq = FacFamSeq;
            vector<vector<int>> tempJobSeqinFam = JobSeqinFam;
            vector<vector<int>> tempm_SpeedMatrix = m_SpeedMatrix;

            for (int i = 1; i < tempFacFamSeq[randomFac].size(); i++)
            {
                int CurFam = FacFamSeq[randomFac][0];

                //�������ڹ���
                randomJobPos1 = rand() % JobSeqinFam[CurFam].size();
                CurJob1 = JobSeqinFam[CurFam][randomJobPos1];

                randomJobPos2 = -1;
                CurJob2 = -1;
                do {
                    randomJobPos2 = rand() % JobSeqinFam[CurFam].size();
                    CurJob2 = JobSeqinFam[CurFam][randomJobPos2];
                } while (CurJob1 == CurJob2);

                JobSeqinFam[CurFam][randomJobPos1] = CurJob2;
                JobSeqinFam[CurFam][randomJobPos2] = CurJob1;

                //�����������ٶ�
                for (int i = 0; i < m_Machines; i++)
                {
                    int tempSpeed;
                    tempSpeed = tempm_SpeedMatrix[CurJob1][i];
                    tempm_SpeedMatrix[CurJob1][i] = tempm_SpeedMatrix[CurJob2][i];
                    tempm_SpeedMatrix[CurJob2][i] = tempSpeed;
                }
            }
            newMS = 0;
            newTEC = 0.0;
            //�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
            for (int j = 0; j < m_Jobs; j++)
            {
                for (int i = 0; i < m_Machines; i++)
                {
                    m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[tempm_SpeedMatrix[j][i]]));
                    UnitPEC[j][i] = 4 * m_Speed[tempm_SpeedMatrix[j][i]] * m_Speed[tempm_SpeedMatrix[j][i]];
                }
            }
            GetMSandTECForPerandToalFac(FacFamSeq, JobSeqinFam, newFacSpan, newFacEC, newMS, newTEC);
            if ((newMS < orgMS) && (newTEC < orgTEC))
            {
                FacFamSeq = tempFacFamSeq;
                JobSeqinFam = tempJobSeqinFam;
                m_SpeedMatrix = tempm_SpeedMatrix;
            }
        }

    }

    else
        return;

}

void MOSA::MOSAEvolution()//randomly setup the mapping relationship, compare with the original relationship, considering span and the number of critical pathes
{
    m_InitTime = Base::GetElapsedProcessTime();
    int count = 0;
    while (true)
    {
        //--------------------------------------------------
        double Tempture = T0;
        count++;
        vector<vector<int>> FacFamSeq = MOSASol.m_FacFamSeqArray;
        vector<vector<int>> JobSeqinFam = MOSASol.m_JobSeqInFamArray;
        int OrgMS = MOSASol.MS;
        float OrgTEC = MOSASol.TEC;
        m_SpeedMatrix = MOSASol.m_SpeedVector;
        //�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }

        vector<vector<int>> bestFacFamSeq = FacFamSeq;
        vector<vector<int>> bestJobSeqinFam = JobSeqinFam;
        vector<vector<int>> bestm_SpeedMatrix = m_SpeedMatrix;
        for (int n = 0; n < iteration_number; n++)
        {
            int random = rand() % 3;

            if (random == 0)
                NGP1(FacFamSeq, JobSeqinFam, m_SpeedMatrix);	//��������1
            else if (random == 1)
                NGP2(FacFamSeq, JobSeqinFam, m_SpeedMatrix);	//��������2
            else
                NGP3(FacFamSeq, JobSeqinFam, m_SpeedMatrix);	//��������3

            int newMS = 0;
            float newTEC = 0.0;
            vector<int> newFacSpan(m_Factories, 0);
            vector<float> newFacEC(m_Factories, 0);
            //�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
            for (int j = 0; j < m_Jobs; j++)
            {
                for (int i = 0; i < m_Machines; i++)
                {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }
            GetMSandTECForPerandToalFac(FacFamSeq, JobSeqinFam, newFacSpan, newFacEC, newMS, newTEC);
            //���
            this->CheckSol(FacFamSeq, JobSeqinFam, newMS);
            this->CheckSolTEC(FacFamSeq, JobSeqinFam, newTEC);

            bool flag = true;
            //����ο���
            if ((newMS < OrgMS) || (newTEC < OrgTEC))
            {
                for (int i = 0; i < MOSAPopulation.size(); i++)
                {
                    if ((newMS == MOSAPopulation[i].MS) && (newTEC == MOSAPopulation[i].TEC))
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
                    tempIndi.m_JobSeqInFamArray = JobSeqinFam;
                    tempIndi.MS = newMS;
                    tempIndi.TEC = newTEC;
                    tempIndi.m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���
                    MOSAPopulation.push_back(tempIndi);
                }
            }

            //�滻
            if ((newMS < OrgMS) && (newTEC < OrgTEC))
            {
                OrgMS = newMS;
                OrgTEC = newTEC;
                bestFacFamSeq = FacFamSeq;
                bestJobSeqinFam = JobSeqinFam;
                bestm_SpeedMatrix = m_SpeedMatrix;

                MOSASol.MS = newMS;
                MOSASol.TEC = newTEC;
                MOSASol.m_FacFamSeqArray = FacFamSeq;
                MOSASol.m_JobSeqInFamArray = JobSeqinFam;
                MOSASol.m_SpeedVector = m_SpeedMatrix;
            }
            else
            {
                double r1 = rand() / (double)RAND_MAX;
                double temp = 0.0;
                double temp1 = OrgMS - newMS;
                double temp2 = OrgTEC - newTEC;
                temp = temp1 > temp2 ? temp1 : temp2;
                if (r1 < exp((temp) / Tempture))
                {
                    OrgMS = newMS;
                    OrgTEC = newTEC;
                    bestFacFamSeq = FacFamSeq;
                    bestJobSeqinFam = JobSeqinFam;
                    bestm_SpeedMatrix = m_SpeedMatrix;
                }
            }
        }

        LS_strategy(bestFacFamSeq, bestJobSeqinFam, bestm_SpeedMatrix);

        int afterMS = 0;
        float afterTEC = 0.0;
        vector<int> afterFacSpan(m_Factories, 0);
        vector<float> afterFacEC(m_Factories, 0);
        //�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[bestm_SpeedMatrix[j][i]] * m_Speed[bestm_SpeedMatrix[j][i]];
            }
        }
        GetMSandTECForPerandToalFac(bestFacFamSeq, bestJobSeqinFam, afterFacSpan, afterFacEC, afterMS, afterTEC);
        //���
        this->CheckSol(bestFacFamSeq, bestJobSeqinFam, afterMS);
        this->CheckSolTEC(bestFacFamSeq, bestJobSeqinFam, afterTEC);

        bool flag = true;
        //����ο���
        if ((afterMS < OrgMS) || (afterTEC < OrgTEC))
        {
            for (int i = 0; i < MOSAPopulation.size(); i++)
            {
                if ((afterMS == MOSAPopulation[i].MS) && (afterTEC == MOSAPopulation[i].TEC))
                {
                    flag = false;
                    break;
                }
            }
            if (flag)
            {
                //cout << "*********************" << endl;
                Individual tempIndi;
                tempIndi.m_FacFamSeqArray = bestFacFamSeq;
                tempIndi.m_JobSeqInFamArray = bestJobSeqinFam;
                tempIndi.MS = afterMS;
                tempIndi.TEC = afterTEC;
                tempIndi.m_SpeedVector = bestm_SpeedMatrix;  //�ٶȾ���
                MOSAPopulation.push_back(tempIndi);
            }
        }

        MOSASol.m_FacFamSeqArray = bestFacFamSeq;
        MOSASol.m_JobSeqInFamArray = bestJobSeqinFam;
        MOSASol.m_SpeedVector = bestm_SpeedMatrix;
        MOSASol.MS = afterMS;
        MOSASol.TEC = afterTEC;

        if (Tempture > TF)
            Tempture = cool_ratio * Tempture;
        else
            Tempture = Tempture;

        //--------------------------------------------------

        long CurTime = Base::GetElapsedProcessTime();
        long ElapsedTime = CurTime - m_InitTime;
        if (ElapsedTime >= m_TimeLimit) break;
    }
   // cout << "��������: " << count << endl;

}

void MOSA::RunEvolution(int CPUFactor, int Rep, vector<vector<Individual>>& MOSAFinalAfterRepParetoSet)
{
    ReadInstanceFileNameList("..\\Benchmark\\");
    int RefSize = 5; //Five paramters are Calibrated on 20190724
    int PS1 = 10;
    int PS2 = 10;
    int AgeLimit = 5;
    int DesLen = 2;
    double Rate = 1.0;
    int Instances = 405;
    vector<int> result;

    string FileDirectory = "..\\Result\\";
    ostringstream str;
    str << "MOSA_" << CPUFactor << "_experiment" << ".txt"; //��ͬ���㷨
    ofstream ofile;
    ofile.open(FileDirectory + str.str());

    for (int ins = 0; ins < Instances; ins++)
    {
        srand((unsigned int)time(NULL));
        ReadInstance(ins);

        //XIN
        vector<Individual> FinalMOSAParetoSet;
        vector<Individual> AfterRepParetoSet;

        vector<Individual> TempAfterRepParetoSet;
        TempAfterRepParetoSet.clear();

        for (int r = 0; r < Rep; r++) {

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

            long StartTime_Ins = Base::GetElapsedProcessTime();
            SetParameters(RefSize, PS1, PS2, AgeLimit, DesLen, CPUFactor * m_Families * m_Machines);

            InitPop();

            MOSAEvolution();

            //��֧���
            MOSAParetoSet.clear();

            Pareto_relation(MOSAPopulation);

            //����֧��� de-emphasize dominated solutions
            for (int j = 0; j < MOSAPopulation.size(); j++)
            {
                for (int i = 0; i < MOSAPopulation.size(); i++)
                {
                    if (MOSAPopulation[i].flag == 0)
                    {
                        if (MOSAPopulation[i].pareto_rel[j] == 1)
                        {
                            MOSAPopulation[i].flag = 999;
                        }
                    }
                }
            }

            for (int i = 0; i < MOSAPopulation.size(); i++)
            {
                if (MOSAPopulation[i].flag == 0)
                    MOSAParetoSet.push_back(MOSAPopulation[i]);
            }


            //����
            vector<Individual> Temp1;
            Temp1.clear();
            for (int i = 0; i < MOSAParetoSet.size(); i++)
            {
                Temp1.push_back(MOSAParetoSet[i]);
            }

            Speed_mutation(MOSAParetoSet, Temp1);

            //ȥ���ظ�
            FinalMOSAParetoSet.clear();
            for (int i = 0; i < MOSAParetoSet.size(); i++)
            {
                bool fg = true;
                for (int j = 0; j < FinalMOSAParetoSet.size(); j++)
                {
                    if ((FinalMOSAParetoSet[j].MS == MOSAParetoSet[i].MS) && (FinalMOSAParetoSet[j].TEC == MOSAParetoSet[i].TEC))
                    {
                        fg = false;
                        break;
                    }
                }
                if (fg)
                {
                    FinalMOSAParetoSet.push_back(MOSAParetoSet[i]);
                }

            }

            long EndTime_Ins = Base::GetElapsedProcessTime();

            //���
            for (int PS = 0; PS < FinalMOSAParetoSet.size(); PS++)
            {
                //m_SpeedMatrix = m_ArchiveSpeedVector[PS];
                m_SpeedMatrix = FinalMOSAParetoSet[PS].m_SpeedVector;
                //�õ������Ĵ���ʱ�� ����λ�ӹ��ܺ�
                for (int j = 0; j < m_Jobs; j++)
                {
                    for (int i = 0; i < m_Machines; i++)
                    {
                        m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
                        UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                    }
                }

                CheckSol(FinalMOSAParetoSet[PS].m_FacFamSeqArray, FinalMOSAParetoSet[PS].m_JobSeqInFamArray, FinalMOSAParetoSet[PS].MS);
                CheckSolTEC(FinalMOSAParetoSet[PS].m_FacFamSeqArray, FinalMOSAParetoSet[PS].m_JobSeqInFamArray, FinalMOSAParetoSet[PS].TEC);
            }


            for (int PS = 0; PS < FinalMOSAParetoSet.size(); PS++)
            {
                TempAfterRepParetoSet.push_back(FinalMOSAParetoSet[PS]);
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


        MOSAFinalAfterRepParetoSet[ins].clear();
        cout << ins + 1 << "\t" << "Factories:" << this->m_Factories<< "\t" << "Machines :" << this->m_Machines<< "\t" << "Families:" << this->m_Families<< "\t" << "Jobs :" << this->m_Jobs << "\t";

        for (int PS = 0; PS < FinalAfterRepParetoSet.size(); PS++)
        {
            // ��� MS �� TEC�����нⶼ����ͬһ��
            cout << "MS:" << FinalAfterRepParetoSet[PS].MS
                 << " " << "TEC:" << FinalAfterRepParetoSet[PS].TEC << "\t";
            ofile << FinalAfterRepParetoSet[PS].MS << "  "<< FinalAfterRepParetoSet[PS].TEC << ","<< "\t";
            MOSAFinalAfterRepParetoSet[ins].push_back(FinalAfterRepParetoSet[PS]);
        }
        ofile << endl;  // ÿ��ʵ�������������
        cout << endl;   // ����̨������У���������ʵ������ʱ����

    }//end ins

    ofile.close();
}