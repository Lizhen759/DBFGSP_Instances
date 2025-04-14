#include <iostream>
#include <numeric>
#include <valarray>
#include "NOperator.h"
#include "Individual.h"
#include <algorithm>
#include <iomanip>

using namespace std;

NOperator::NOperator()
{
}

NOperator::~NOperator()
{
}

void NOperator::CheckSol(vector <int> FamSeq, vector <vector <int>> JobSeqinFam, int Span)
{
    if (!FamSeq.size())
    {
        cout << "Error in CheckSol: FamSeq.size()=0" << endl;
        getchar();
        exit(0);
    }

    //Check Families
    vector <bool> bFamArray(this->m_Families, false);

    for (int Fam = 0; Fam < FamSeq.size(); Fam++)
        bFamArray[FamSeq[Fam]] = true;
    for (int i = 0; i < this->m_Families; i++)
    {
        if (!bFamArray[i])
        {
            cout << "Family " << i << " is lost!" << endl;
            getchar();
            exit(0);
        }
    }

    //Check Jobs
    vector <bool> bJobArray(this->m_Jobs, false);
    for (int Fam = 0; Fam < JobSeqinFam.size(); Fam++)
    {
        for (int j = 0; j < JobSeqinFam[Fam].size(); j++)
        {
            int Job = JobSeqinFam[Fam][j];
            bJobArray[Job] = true;
        }
    }

    for (int i = 0; i < this->m_Jobs; i++)
    {
        if (!bJobArray[i])
        {
            cout << "Job " << i << " is lost!" << endl;
            getchar();
            exit(0);
        }
    }

    //Check Jobs in each Family
    for (int Fam = 0; Fam < this->m_JobSeqInFam.size(); Fam++)
    {
        for (int j = 0; j < this->m_JobSeqInFam[Fam].size(); j++)
        {
            int Job = this->m_JobSeqInFam[Fam][j];
            if (JobSeqinFam[Fam].end() == find(JobSeqinFam[Fam].begin(), JobSeqinFam[Fam].end(), Job))
            {
                cout << "Job is not found in Family:\t" << Fam << endl;
                getchar();
                exit(0);
            }
        }
    }

    //----- Check makespan---------
    int TSpan = this->GetSpan(FamSeq, JobSeqinFam);
    if (TSpan != Span)
    {
        cout << "Span is Erro!" << Span << "\t" << TSpan << endl;
        getchar();
        exit(0);
    }
    cout << "��ȷ��" << endl;
}
//��
void NOperator::CheckSol(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam, int Span)
{
    if (!FacFamSeq.size())
    {
        cout << "Error in CheckSol: FamSeq.size()=0" << endl;
        getchar();
        exit(0);
    }
    vector<bool> bFamArray(this->m_Families, false); //Check Families
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
        for (int Fam = 0; Fam < FacFamSeq[Fac].size(); Fam++)
            bFamArray[FacFamSeq[Fac][Fam]] = true;
    for (int i = 0; i < this->m_Families; i++)
    {
        if (!bFamArray[i])
        {
            cout << "Family " << i << " is lost!" << endl;
            getchar();
            exit(0);
        }
    }
    vector<bool> bJobArray(this->m_Jobs, false); //Check Jobs
    for (int Fam = 0; Fam < JobSeqInFam.size(); Fam++)
    {
        for (int j = 0; j < JobSeqInFam[Fam].size(); j++)
        {
            int Job = JobSeqInFam[Fam][j];
            bJobArray[Job] = true;
        }
    }
    for (int i = 0; i < this->m_Jobs; i++)
    {
        if (!bJobArray[i])
        {
            cout << "Job " << i << " is lost!" << endl;
            getchar();
            exit(0);
        }
    }
    for (int Fam = 0; Fam < this->m_JobSeqInFam.size(); Fam++) //Check Jobs in each Family
    {
        for (int j = 0; j < this->m_JobSeqInFam[Fam].size(); j++)
        {
            int Job = this->m_JobSeqInFam[Fam][j];
            if (JobSeqInFam[Fam].end() == find(JobSeqInFam[Fam].begin(), JobSeqInFam[Fam].end(), Job))
            {
                cout << "Job is not found in Family:\t" << Fam << endl;
                getchar();
                exit(0);
            }
        }
    }

    //----- Check makespan---------
    int TSpan = this->GetSpan(FacFamSeq, JobSeqInFam);  //�õ��깤ʱ��
    if (TSpan != Span)
    {
        cout << "Span is Erro!" << Span << "\t" << TSpan << endl;
        getchar();
        exit(0);
    }

}
//xin
int NOperator::GetSpan(vector <vector <int>> FacFamSeq, vector <vector <int>> JobSeqInFam)
{
    vector<vector<int>> JFDTime(this->m_Jobs);
    for (int j = 0; j < JFDTime.size(); j++)
        JFDTime[j].resize(this->m_Machines, 0);
    vector<int> FacSpan(FacFamSeq.size(), 0);
    int mSpan = this->GetJFDTime_Forward(FacFamSeq, JobSeqInFam, JFDTime, FacSpan);   //����GetJFDTime_Forward
    return mSpan;
}
int NOperator::GetSpan(vector <int> FamSeq, vector <vector <int>> JobSeqinFam)
{
    vector<vector<int>> JFDTime(this->m_Jobs);
    for (int j = 0; j < JFDTime.size(); j++)
        JFDTime[j].resize(this->m_Machines, 0);
    int mSpan = GetJFDTime_Forward_InFactory(FamSeq, JobSeqinFam, JFDTime);
    return mSpan;
}

//�õ��ܵ��깤ʱ�䣨�������أ�20190524
int NOperator::GetSpan(vector <vector <int>> FamSeq, vector <vector <int>> JobSeqinFam, vector<int>& FacSpan)
{
    vector<vector<int>> JFDTime(this->m_Jobs);
    for (int j = 0; j < JFDTime.size(); j++)
        JFDTime[j].resize(this->m_Machines, 0);
    int mSpan = this->GetJFDTime_Forward(FamSeq, JobSeqinFam, JFDTime, FacSpan);
    return mSpan;
}

void NOperator::CheckSolTEC(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam, float TEC)
{
    if (!FacFamSeq.size())
    {
        cout << "Error in CheckSol: FamSeq.size()=0" << endl;
        getchar();
        exit(0);
    }
    vector<bool> bFamArray(this->m_Families, false); //Check Families
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
        for (int Fam = 0; Fam < FacFamSeq[Fac].size(); Fam++)
            bFamArray[FacFamSeq[Fac][Fam]] = true;
    for (int i = 0; i < this->m_Families; i++)
    {
        if (!bFamArray[i])
        {
            cout << "Family " << i << " is lost!" << endl;
            getchar();
            exit(0);
        }
    }
    vector<bool> bJobArray(this->m_Jobs, false); //Check Jobs
    for (int Fam = 0; Fam < JobSeqInFam.size(); Fam++)
    {
        for (int j = 0; j < JobSeqInFam[Fam].size(); j++)
        {
            int Job = JobSeqInFam[Fam][j];
            bJobArray[Job] = true;
        }
    }
    for (int i = 0; i < this->m_Jobs; i++)
    {
        if (!bJobArray[i])
        {
            cout << "Job " << i << " is lost!" << endl;
            getchar();
            exit(0);
        }
    }
    for (int Fam = 0; Fam < this->m_JobSeqInFam.size(); Fam++) //Check Jobs in each Family
    {
        for (int j = 0; j < this->m_JobSeqInFam[Fam].size(); j++)
        {
            int Job = this->m_JobSeqInFam[Fam][j];
            if (JobSeqInFam[Fam].end() == find(JobSeqInFam[Fam].begin(), JobSeqInFam[Fam].end(), Job))
            {
                cout << "Job is not found in Family:\t" << Fam << endl;
                getchar();
                exit(0);
            }
        }
    }

    //----- Check TEC---------
    float tempTEC = this->GetTEC(FacFamSeq, JobSeqInFam);  //�õ��깤ʱ��

    if (tempTEC != TEC)
    {
        cout.precision(10);
        cout << "TEC is Erro!" << "TEC:" << TEC << "\t" << "tempTEC:" << tempTEC << endl;
        getchar();
        exit(0);
    }
}
//�� �õ����ܺ�
float NOperator::GetTEC(vector <vector <int>> FacFamSeq, vector <vector <int>> JobSeqInFam)
{
    double TEC = 0;
    vector<vector<int>> JFDTime(this->m_Jobs);
    for (int j = 0; j < JFDTime.size(); j++)
        JFDTime[j].resize(this->m_Machines, 0);
    vector<float> FacEC(FacFamSeq.size(), 0);
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        FacEC[Fac] = GetPerFacEC(FacFamSeq[Fac], JobSeqInFam, JFDTime);//�õ�ÿ��������EC
        TEC += FacEC[Fac];
    }
    return TEC;
}

//��
float NOperator::GetPerFacEC(vector <int> FamSeqInFac, vector <vector <int>> JobSeqInFam, vector<vector<int>>& JFDTime)
{

    int Span = GetJFDTime_Forward_InFactory(FamSeqInFac, JobSeqInFam, JFDTime);//�õ�ÿ���������깤ʱ��

    //�����ܺ�
    double TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

    for (int m = 0; m < m_Machines; m++)
    {
        int preFam = -1;
        //int CurFam = -1;
        int Tptime = 0, Tsetuptime = 0, Tidletime = 0;

        int Fam = FamSeqInFac[0];
        for (int g = 0; g < FamSeqInFac.size(); g++)
        {
            Fam = FamSeqInFac[g];

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
        Tidletime += (Span - Tptime - Tsetuptime);
        TIdleTimeEC += Tidletime * UnitIdleEC;
    }
    double FacEC = TSetupTimeEC + TIdleTimeEC + TpTimeEC;
    return FacEC;
}


/**
 * ����ÿ�������������뿪ʱ��
 * @param FamSeq
 * @param JobSeqinFam
 * @param JFDTime
 * @param Span
 */
int NOperator::GetJFDTime_Forward(const vector<vector<int>> &FamSeq, const vector<vector<int>> &JobSeqInFam,
                                  vector<vector<int>> &JFDTime, vector<int> &FacSpan)
{
    for (int Fac = 0; Fac < FamSeq.size(); Fac++)
        FacSpan[Fac] = GetJFDTime_Forward_InFactory(FamSeq[Fac], JobSeqInFam, JFDTime);
    //�õ�ÿ���������뿪ʱ��
    return *max_element(FacSpan.begin(), FacSpan.end()); //�õ��ؼ�����
}

/**
 * ����ÿ�������ڹ����������뿪ʱ��
 * @param FamSeqInFac
 * @param JobSeqinFam
 * @param JFDTime
 * @return
 */
int
NOperator::GetJFDTime_Forward_InFactory(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,vector<vector<int>> &JFDTime)
{

    for (int Fam = 0; Fam < FamSeqInFac.size(); Fam++)
    {
        int CurFam = FamSeqInFac[Fam];//��ǰ��
        if (Fam == 0)//the first group of jobs ��һ����
        {
            for (int m = 0; m < this->m_Machines; m++)
                MachReadyTime[m] = this->m_SetupTime[m][CurFam][CurFam];
        }
        else //�ڶ�������� from the second group of jobs to the end;
        {
            int PreFam = FamSeqInFac[Fam - 1];
            for (int m = 0; m < this->m_Machines; m++)
                MachReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
        }

        for (int j = 0; j < JobSeqInFam[CurFam].size(); j++)//��ǰ��Ĺ������� Scheduling Jobs in CurFam
        {
            int CurJob = JobSeqInFam[CurFam][j];//��ǰ����
            JFDTime[CurJob][0] = max(MachReadyTime[0] + this-> m_TureJobOpertime[CurJob][0], MachReadyTime[1]);

            for (int m = 1; m < this->m_Machines; m++)
            {
                if (m == this->m_Machines - 1)
                {
                    JFDTime[CurJob][m] = JFDTime[CurJob][m - 1] + this->m_TureJobOpertime[CurJob][m];
                }
                else
                {
                    JFDTime[CurJob][m] = max(JFDTime[CurJob][m - 1] + m_TureJobOpertime[CurJob][m],
                                             MachReadyTime[m + 1]);
//                    cout<<JFDTime[CurJob][m]<<endl;
                }
            }
            MachReadyTime = JFDTime[CurJob];//����
        }
    }
    return MachReadyTime[this->m_Machines - 1];
}

// �������õ��ؼ������깤ʱ��
int NOperator::GetJBDTime_Backward(const vector<vector<int>> &FamSeq, const vector<vector<int>> &JobSeqInFam,
                                   vector<vector<int>> &JBDTime, vector<int> &FacSpan) //Backward pass calculation
{
    for (int Fac = 0; Fac < FamSeq.size(); Fac++)
        FacSpan[Fac] = this->GetJBDTime_Backward_InFactory(FamSeq[Fac], JobSeqInFam, JBDTime);//�õ�ÿ���������깤ʱ��
    return *max_element(FacSpan.begin(), FacSpan.end());
}
/**
 * �������õ������깤ʱ��
 * @param FamSeqInFac
 * @param JobSeqinFam
 * @param JBDTime
 * @return
 */
int NOperator::GetJBDTime_Backward_InFactory(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                             vector<vector<int>> &JBDTime)
{
    for (int Fam = FamSeqInFac.size() - 1; Fam >= 0; Fam--)  //������һ���鿪ʼ
    {
        int CurFam = FamSeqInFac[Fam];  //��ǰ��
        if (Fam < FamSeqInFac.size() - 1)//�������һ���� Not the last group of jobs
        {
            int NextFam = FamSeqInFac[Fam + 1];  //��һ����
            for (int m = this->m_Machines - 1; m >= 0; m--)
                MachReadyTime[m] += this->m_SetupTime[m][CurFam][NextFam];  //׼��ʱ��
        }

        for (int j = JobSeqInFam[CurFam].size() - 1; j >= 0; j--)//���������һ��������ʼ Scheduling Jobs in CurFam
        {
            int CurJob = JobSeqInFam[CurFam][j];
            JBDTime[CurJob][this->m_Machines - 1] = max(
                    MachReadyTime[this->m_Machines - 1] +  m_TureJobOpertime[CurJob][this->m_Machines - 1], MachReadyTime[this->m_Machines - 2]);//���һ̨����on the last machine
            for (int m = this->m_Machines - 2; m >= 0; m--) // ʣ�µĻ��� on the rest machine
            {
                if (m == 0)
                {
                    JBDTime[CurJob][m] = JBDTime[CurJob][m + 1] + m_TureJobOpertime[CurJob][m];
                } else
                {
                    JBDTime[CurJob][m] = max(JBDTime[CurJob][m + 1] + m_TureJobOpertime[CurJob][m], MachReadyTime[m - 1]);
                }
            }
            MachReadyTime = JBDTime[CurJob];
        }

        if (Fam == 0) // ��һ���� the first Family
        {
            for (int m = this->m_Machines - 1; m >= 0; m--)
                MachReadyTime[m] += this->m_SetupTime[m][CurFam][CurFam];
        }
    }
    return MachReadyTime[0];;
}

int NOperator::GetSpanForPerFacAfterInsertFam(const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam,
                                              const vector<vector<int>>& JFDTime, const vector<vector<int>>& JBDTime, int InsFam, int Pos)
{
    if (Pos == 0)
    {
        for (int m = 0; m < this->m_Machines; m++)
        {
            MachReadyTime[m] = this->m_SetupTime[m][InsFam][InsFam]; //MachReadyTime ��Ϊ׼��ʱ��
        }
    }
    else
    {
        int PreFam = FamSeqInFac[Pos - 1]; //ǰһ����
        int LastJobinPreFam = JobSeqInFam[PreFam][JobSeqInFam[PreFam].size() - 1]; //ǰһ��������һ������
        for (int m = 0; m < this->m_Machines; m++)
        {
            MachReadyTime[m] = JFDTime[LastJobinPreFam][m] + this->m_SetupTime[m][PreFam][InsFam];
        }
    }

    // ��������﹤������ Scheduling jobs in InsFam
    for (int j = 0; j < JobSeqInFam[InsFam].size(); j++)
    {
        int CurJob = JobSeqInFam[InsFam][j]; //��ǰ����
        //�ڵ�һ��������
        MachReadyTime[0] = max(MachReadyTime[0] + this->m_TureJobOpertime[CurJob][0], MachReadyTime[1]);

        for (int m = 1; m < this->m_Machines; m++)
        {
            if (m == this->m_Machines - 1)
            {
                MachReadyTime[m] = MachReadyTime[m - 1] + this->m_TureJobOpertime[CurJob][m];
            }
            else
            {
                MachReadyTime[m] = max(MachReadyTime[m - 1] + m_TureJobOpertime[CurJob][m],MachReadyTime[m + 1]);
            }
        }
    }

    int Span = 0;
    if (Pos < FamSeqInFac.size())
    {
        int NextFam = FamSeqInFac[Pos]; //��һ����

        int FirstJobinNextFam = JobSeqInFam[NextFam][0]; //��һ����ĵ�һ������

        for (int m = 0; m < this->m_Machines; m++)
        {
            if (MachReadyTime[m] + this->m_SetupTime[m][InsFam][NextFam] + JFDTime[FirstJobinNextFam][m] >= Span)
            {
                Span = MachReadyTime[m] + this->m_SetupTime[m][InsFam][NextFam] + JBDTime[FirstJobinNextFam][m];
            }
        }
    }
    else
    {
        Span = MachReadyTime[this->m_Machines - 1];
    }
    return Span;
}

/**
 * ����������TEC
 * @param FamSeqInFac
 * @param JobSeqinFam
 * @param JFDTime
 * @param JBDTime
 * @param InsFam
 * @param Pos
 * @return
 */
float NOperator::GetECForPerFacAfterInsertFam(const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam,
                                              const vector<vector<int>>& JFDTime, const vector<vector<int>>& JBDTime, int InsFam, int Pos)
{
    if (Pos == 0)
    {
        for (int m = 0; m < this->m_Machines; m++)
        {
            MachReadyTime[m] = this->m_SetupTime[m][InsFam][InsFam]; //MachReadyTime ��Ϊ׼��ʱ��
        }
    }
    else
    {
        int PreFam = FamSeqInFac[Pos - 1]; //ǰһ����
        int LastJobinPreFam = JobSeqInFam[PreFam][JobSeqInFam[PreFam].size() - 1]; //ǰһ��������һ������
        for (int m = 0; m < this->m_Machines; m++)
        {
            MachReadyTime[m] = JFDTime[LastJobinPreFam][m] + this->m_SetupTime[m][PreFam][InsFam];
        }
    }

    // ��������﹤������ Scheduling jobs in InsFam
    for (int j = 0; j < JobSeqInFam[InsFam].size(); j++)
    {
        int CurJob = JobSeqInFam[InsFam][j]; //��ǰ����
        //�ڵ�һ��������
        MachReadyTime[0] = max(MachReadyTime[0] + this->m_TureJobOpertime[CurJob][0], MachReadyTime[1]);

        for (int m = 1; m < this->m_Machines; m++)
        {
            if (m == this->m_Machines - 1)
            {
                MachReadyTime[m] = MachReadyTime[m - 1] + this->m_TureJobOpertime[CurJob][m];
            }
            else
            {
                MachReadyTime[m] = max(MachReadyTime[m - 1] + m_TureJobOpertime[CurJob][m],MachReadyTime[m + 1]);
            }
        }
    }

    int Span = 0;
    if (Pos < FamSeqInFac.size())
    {
        int NextFam = FamSeqInFac[Pos]; //��һ����

        int FirstJobinNextFam = JobSeqInFam[NextFam][0]; //��һ����ĵ�һ������

        for (int m = 0; m < this->m_Machines; m++)
        {
            if (MachReadyTime[m] + this->m_SetupTime[m][InsFam][NextFam] + JFDTime[FirstJobinNextFam][m] >= Span)
            {
                Span = MachReadyTime[m] + this->m_SetupTime[m][InsFam][NextFam] + JBDTime[FirstJobinNextFam][m];
            }
        }
    }
    else
    {
        Span = MachReadyTime[this->m_Machines - 1];
    }

    //�����ܺ�
    float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

    vector<int> TempFamSeqInFac;
    TempFamSeqInFac = FamSeqInFac;
    TempFamSeqInFac.insert(TempFamSeqInFac.begin() + Pos, InsFam);

    for (int m = 0; m < m_Machines; m++)
    {
        int preFam = -1;
        int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
        int Fam = TempFamSeqInFac[0];
        for (int g = 0; g < TempFamSeqInFac.size(); g++)
        {
            Fam = TempFamSeqInFac[g];

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
        Tidletime += (Span - Tptime - Tsetuptime);
        TIdleTimeEC += Tidletime * UnitIdleEC;
    }
    float FacEC = TSetupTimeEC + TIdleTimeEC + TpTimeEC;
    return FacEC;
}


/**
 * ���ڼ��ټ�������ǰ���򹤼��뿪ʱ��
 * @param FamSeqinFac
 * @param JobSeqinFam
 * @param JFDTime
 * @param JBDTime
 * @param ForwardFamPos
 * @param BackwardFamPos
 * @param ForwardJobPos
 * @param BackwardJobPos
 * @return
 */
void NOperator::RefreshJFDJBDTime_InFactory(const vector<int> &FamSeqinFac, const vector<vector<int>> &JobSeqinFam,
                                            vector<vector<int>> &JFDTime, vector<vector<int>> &JBDTime,
                                            int ForwardFamPos, int BackwardFamPos, int ForwardJobPos,
                                            int BackwardJobPos)
{
    if (FamSeqinFac.empty())
    {
        return;
    }

    vector<int> ReadyTime(this->m_Machines, 0);//����׼��ʱ���ʼ��Ϊ0

    if (ForwardFamPos <= FamSeqinFac.size() - 1)
    {
        int ForwardFam = FamSeqinFac[ForwardFamPos];
        if (ForwardFamPos == 0)
        {
            //��һ����
            if (ForwardJobPos == 0)
            {
                //��һ������
                for (int m = 0; m < this->m_Machines; ++m)
                {
                    ReadyTime[m] = this->m_SetupTime[m][ForwardFam][ForwardFam];
                }
            } else
            {
                //���ǵ�һ������
                int PreJob = JobSeqinFam[ForwardFam][ForwardJobPos - 1];
                //����׼��ʱ��=ǰһ���������뿪ʱ��
                copy(begin(JFDTime[PreJob]), end(JFDTime[PreJob]), begin(ReadyTime));
            }
        } else
        {
            //���ǵ�һ����
            if (ForwardJobPos == 0)
            {
                //��һ������
                int PreFam = FamSeqinFac[ForwardFamPos - 1];
                int LastJobInPreFam = *JobSeqinFam[PreFam].rbegin();
                for (int m = 0; m < this->m_Machines; ++m)
                {
                    ReadyTime[m] = JFDTime[LastJobInPreFam][m] + this->m_SetupTime[m][PreFam][ForwardFam];
                }
            } else
            {
                //���ǵ�һ������
                int PreJob = JobSeqinFam[ForwardFam][ForwardJobPos - 1];
                copy(begin(JFDTime[PreJob]), end(JFDTime[PreJob]), begin(ReadyTime));
            }
        }

        for (int j = ForwardJobPos; j < JobSeqinFam[ForwardFam].size(); ++j)
        {
            //��ǰ��Ĺ�������
            int CurJob = JobSeqinFam[ForwardFam][j];//��ǰ����
            JFDTime[CurJob][0] = max(ReadyTime[0] + this->m_TureJobOpertime[CurJob][0], ReadyTime[1]);//�ڵ�һ̨�����ϵ��뿪ʱ��
            for (int m = 1; m < this->m_Machines; m++)
            {
                //�����������ϵ��뿪ʱ��
                if (m == this->m_Machines - 1)
                {
                    //���һ̨����
                    JFDTime[CurJob][m] = JFDTime[CurJob][m - 1] + this->m_TureJobOpertime[CurJob][m];
                } else
                {
                    JFDTime[CurJob][m] = max(JFDTime[CurJob][m - 1] + this->m_TureJobOpertime[CurJob][m],
                                             ReadyTime[m + 1]);
                }
            }
            copy(begin(JFDTime[CurJob]), end(JFDTime[CurJob]), begin(ReadyTime));
        }

        int PreFam = FamSeqinFac[ForwardFamPos];
        for (int Fam = ForwardFamPos + 1; Fam < FamSeqinFac.size(); Fam++)
        {
            int CurFam = FamSeqinFac[Fam];//��ǰ��
            for (int m = 0; m < this->m_Machines; ++m)
            {
                ReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
            }

            for (int j = 0; j < JobSeqinFam[CurFam].size(); ++j)
            {
                int CurJob = JobSeqinFam[CurFam][j];//��ǰ����
                JFDTime[CurJob][0] = max(ReadyTime[0] + this->m_TureJobOpertime[CurJob][0], ReadyTime[1]);//�ڵ�һ̨�����ϵ��뿪ʱ��
                for (int m = 1; m < this->m_Machines; m++)
                {
                    //�����������ϵ��뿪ʱ��
                    if (m == this->m_Machines - 1)
                    {
                        //���һ̨����
                        JFDTime[CurJob][m] = JFDTime[CurJob][m - 1] + this->m_TureJobOpertime[CurJob][m];
                    } else
                    {
                        JFDTime[CurJob][m] = max(JFDTime[CurJob][m - 1] + this->m_TureJobOpertime[CurJob][m],
                                                 ReadyTime[m + 1]);
                    }
                }
                copy(begin(JFDTime[CurJob]), end(JFDTime[CurJob]), begin(ReadyTime));
            }

            PreFam = CurFam;
        }
    }

    if (BackwardFamPos >= 0)
    {
        fill(begin(ReadyTime), end(ReadyTime), 0);//��ʼ������׼��ʱ��
        int BackwardFam = FamSeqinFac[BackwardFamPos]; //backward family
        if (BackwardFamPos == FamSeqinFac.size() - 1)
        {
            //���һ����
            if (BackwardJobPos != JobSeqinFam[BackwardFam].size() - 1)
            {
                //�������һ������
                int NextJob = JobSeqinFam[BackwardFam][BackwardJobPos + 1];
                for (int m = this->m_Machines - 1; m >= 0; m--)
                {
                    ReadyTime[m] = JBDTime[NextJob][m];
                }
            }
        } else
        {
            //�������һ����
            if (BackwardJobPos != JobSeqinFam[BackwardFam].size() - 1)
            {
                //�������һ������
                int NextJob = JobSeqinFam[BackwardFam][BackwardJobPos + 1];
                for (int m = this->m_Machines - 1; m >= 0; m--)
                {
                    ReadyTime[m] = JBDTime[NextJob][m];
                }
            } else
            {
                //���һ������
                int NextFam = FamSeqinFac[BackwardFamPos + 1]; //��һ����
                int FitstJobInNextFam = JobSeqinFam[NextFam][0];
                for (int m = this->m_Machines - 1; m >= 0; m--)
                {
                    ReadyTime[m] = JBDTime[FitstJobInNextFam][m] + this->m_SetupTime[m][BackwardFam][NextFam]; //׼��ʱ��
                }
            }
        }

        for (int j = BackwardJobPos; j >= 0; j--)
        {
            //�����һ��������ʼ
            int CurJob = JobSeqinFam[BackwardFam][j];
            //���һ̨����
            JBDTime[CurJob][this->m_Machines - 1] = max(
                    ReadyTime[this->m_Machines - 1] + this->m_TureJobOpertime[CurJob][this->m_Machines - 1],
                    ReadyTime[this->m_Machines - 2]);
            for (int m = this->m_Machines - 2; m >= 0; m--)
            {
                //ʣ��Ļ���
                if (m == 0)
                {
                    //��һ̨����
                    JBDTime[CurJob][m] = JBDTime[CurJob][m + 1] + this->m_TureJobOpertime[CurJob][m];
                } else
                {
                    JBDTime[CurJob][m] = max(JBDTime[CurJob][m + 1] + this->m_TureJobOpertime[CurJob][m],
                                             ReadyTime[m - 1]);
                }
            }
            copy(begin(JBDTime[CurJob]), end(JBDTime[CurJob]), begin(ReadyTime));
        }

        int NextFam = FamSeqinFac[BackwardFamPos];
        for (int Fam = BackwardFamPos - 1; Fam >= 0; Fam--)
        {
            int CurFam = FamSeqinFac[Fam]; //��ǰ��
            for (int m = this->m_Machines - 1; m >= 0; m--)
            {
                ReadyTime[m] += this->m_SetupTime[m][CurFam][NextFam];
            }

            for (int j = JobSeqinFam[CurFam].size() - 1; j >= 0; j--)
            {
                //���������һ��������ʼ
                int CurJob = JobSeqinFam[CurFam][j];
                //���һ̨����
                JBDTime[CurJob][this->m_Machines - 1] = max(
                        ReadyTime[this->m_Machines - 1] + this->m_TureJobOpertime[CurJob][this->m_Machines - 1],
                        ReadyTime[this->m_Machines - 2]);
                for (int m = this->m_Machines - 2; m >= 0; m--)
                {
                    //ʣ��Ļ���
                    if (m == 0)
                    {
                        //��һ̨����
                        JBDTime[CurJob][m] = JBDTime[CurJob][m + 1] + this->m_TureJobOpertime[CurJob][m];
                    } else
                    {
                        JBDTime[CurJob][m] = max(JBDTime[CurJob][m + 1] + this->m_TureJobOpertime[CurJob][m],
                                                 ReadyTime[m - 1]);
                    }
                }
                copy(begin(JBDTime[CurJob]), end(JBDTime[CurJob]), begin(ReadyTime));
            }
            NextFam = CurFam;
        }
    }

}


void NOperator::CopyJFDJBDTime(const vector<vector<int>>& srcJFDTime, const vector<vector<int>>& srcJBDTime,
                               vector<vector<int>>& dstJFDTime, vector<vector<int>>& dstJBDTime)
{
    for (int i = 0; i < m_Jobs; ++i)
    {
        copy(begin(srcJFDTime[i]), end(srcJFDTime[i]), begin(dstJFDTime[i]));
        copy(begin(srcJBDTime[i]), end(srcJBDTime[i]), begin(dstJBDTime[i]));
    }
}

// Select the factory with lowest makespan after the factory includes the appended Fam
int NOperator::GetSol_Include(vector<int> FamPrmu, vector<vector<int>>JobSeqInFam, vector<vector<int>>& FacFamSeq, vector<int>& FacSpan)
{
    FacFamSeq.clear();
    FacFamSeq.resize(this->m_Factories);
    FacSpan.clear();
    FacSpan.resize(this->m_Factories, 0);

    vector<vector<int>> FacMachReadyTime(this->m_Factories);
    for (int Fac = 0; Fac < FacMachReadyTime.size(); Fac++)
        FacMachReadyTime[Fac].resize(this->m_Machines, 0);

    vector<vector<int>> JFDTime(this->m_Jobs);
    for (int j = 0; j < this->m_Jobs; j++)
        JFDTime[j].resize(this->m_Machines, 0);

    // Assign one job to each Factory
    int MaxFac = min(this->m_Factories,static_cast<int>(FamPrmu.size()));
    for (int Fac = 0; Fac < MaxFac; Fac++)
    {
        int CurFam = FamPrmu[Fac];

        // Family setup: initialize the ready time for each machine in the factory
        for (int m = 0; m < this->m_Machines; m++)
            FacMachReadyTime[Fac][m] = this->m_SetupTime[m][CurFam][CurFam];  // Setup time for the first group

        // Process jobs one by one in the Family
        for (int j = 0; j < JobSeqInFam[CurFam].size(); j++)
        {
            int CurJob = JobSeqInFam[CurFam][j];

            // The first machine's completion time
            JFDTime[CurJob][0] = max(FacMachReadyTime[Fac][0] + this->m_TureJobOpertime[CurJob][0], FacMachReadyTime[Fac][1]);

            // Iterate through the remaining machines
            for (int m = 1; m < this->m_Machines; m++)
            {
                // For the last machine, just add operation time to the previous machine's completion time
                if (m == this->m_Machines - 1)
                {
                    JFDTime[CurJob][m] =  JFDTime[CurJob][m - 1] + this->m_TureJobOpertime[CurJob][m];
                }
                else
                {
                    // For the other machines, check the max between the ready time and previous job's completion time
                    JFDTime[CurJob][m] = max(JFDTime[CurJob][m - 1] + this->m_TureJobOpertime[CurJob][m], FacMachReadyTime[Fac][m + 1]);
                }
            }
            // After processing the job, update the ready time for each machine in the factory
            FacMachReadyTime[Fac] = JFDTime[CurJob];
        }

        // Update the factory's sequence and span
        FacFamSeq[Fac].push_back(CurFam);
        FacSpan[Fac] = FacMachReadyTime[Fac][this->m_Machines - 1];
    }

    //Assign the remaining Jobs
    vector<vector<int>> TempFacMachReadyTime = FacMachReadyTime;
    for (int Fam = this->m_Factories; Fam < FamPrmu.size(); Fam++)
    {
        GetFamSumSetupTime();
        vector<vector<int>>D = this->m_FamSumSetupTime;

        int CurFam = FamPrmu[Fam];
        int SelFac = -1, minFacSpan = INT_MAX;
        for (int Fac = 0; Fac < this->m_Factories; Fac++) //find the Factory that can complete the Family earliest
        {
            int PreFam = FacFamSeq[Fac][FacFamSeq[Fac].size() - 1];
            int CurFam = FamPrmu[Fam];
            for (int m = 0; m < this->m_Machines; m++)//Family setup time
            {
                TempFacMachReadyTime[Fac][m] = FacMachReadyTime[Fac][m] + this->m_SetupTime[m][PreFam][CurFam];
            }

            for (int j = 0; j < JobSeqInFam[CurFam].size(); j++) //Process Jobs one by one in the Family
            {
                int CurJob = JobSeqInFam[CurFam][j];
                JFDTime[CurJob][0] = max(FacMachReadyTime[Fac][0] + this->m_TureJobOpertime[CurJob][0], FacMachReadyTime[Fac][1]);

                for (int m = 1; m < this->m_Machines; m++)
                {
                    // For the last machine, just add operation time to the previous machine's completion time
                    if (m == this->m_Machines - 1)
                    {
                        JFDTime[CurJob][m] =  JFDTime[CurJob][m - 1] + this->m_TureJobOpertime[CurJob][m];
                    }
                    else
                    {
                        // For the other machines, check the max between the ready time and previous job's completion time
                        JFDTime[CurJob][m] = max(JFDTime[CurJob][m - 1] + this->m_TureJobOpertime[CurJob][m], FacMachReadyTime[Fac][m + 1]);
                    }
                }
                // After processing the job, update the ready time for each machine in the factory
                FacMachReadyTime[Fac] = JFDTime[CurJob];
            }

            if (TempFacMachReadyTime[Fac][this->m_Machines - 1] < minFacSpan) //Record the best Factory
            {
                minFacSpan = TempFacMachReadyTime[Fac][this->m_Machines - 1];
                SelFac = Fac;
            }
        }

        //Assign CurFam to the selected factory
        FacMachReadyTime[SelFac] = TempFacMachReadyTime[SelFac];
        FacFamSeq[SelFac].push_back(CurFam);
        FacSpan[SelFac] = minFacSpan;
    }
    return *max_element(FacSpan.begin(), FacSpan.end());
}

int NOperator::GetDelayTime_Forward(vector <vector <int>> FacFamSeq, vector <vector <int>> JobSeqInFam, vector<vector<int>>& JFDTime, vector<int>& FacSpan, vector<vector<int>>& DelayTime)//Forward pass calculation
{
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
        FacSpan[Fac] = GetDelayTime_Forward_InFactory(FacFamSeq[Fac], JobSeqInFam, JFDTime, DelayTime);//�õ�ÿ���������깤ʱ��
    return *max_element(FacSpan.begin(), FacSpan.end()); //�õ��ؼ�����
}
int NOperator::GetDelayTime_Forward_InFactory(vector <int> FamSeqInFac, vector <vector <int>> JobSeqinFam, vector<vector<int>>& JFDTime, vector<vector<int>>& DelayTime)
{
    for (int Fam = 0; Fam < FamSeqInFac.size(); Fam++)
    {
        int CurFam = FamSeqInFac[Fam];  // ��ǰ��
        if (Fam == 0)  // ��һ����
        {
            for (int m = 0; m < this->m_Machines; m++)
                MachReadyTime[m] = this->m_SetupTime[m][CurFam][CurFam];
        }
        else  // �ӵڶ��鵽�����
        {
            int PreFam = FamSeqInFac[Fam - 1];
            for (int m = 0; m < this->m_Machines; m++)
                MachReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
        }

        for (int j = 0; j < JobSeqinFam[CurFam].size(); j++)  // ��ǰ��Ĺ�������
        {
            int CurJob = JobSeqinFam[CurFam][j];  // ��ǰ����
            // ���ڵ�һ�����������㵱ǰ�������깤ʱ��
            JFDTime[CurJob][0] = max(MachReadyTime[0] + this->m_TureJobOpertime[CurJob][0], MachReadyTime[1]);

            for (int m = 1; m < this->m_Machines; m++)
            {
                // �������һ̨������ֱ�Ӽ����깤ʱ��
                if (m == this->m_Machines - 1)
                {
                    JFDTime[CurJob][m] = JFDTime[CurJob][m - 1] + this->m_TureJobOpertime[CurJob][m];
                }
                else
                {
                    // �����������Ͽ������������
                    JFDTime[CurJob][m] = max(JFDTime[CurJob][m - 1] + m_TureJobOpertime[CurJob][m], MachReadyTime[m + 1]);
                }

                // ��¼ DelayTime[PreJob][m] �ڴ˴�
                if (j > 0)  // ������ǵ�һ������
                {
                    int PreJob = JobSeqinFam[CurFam][j - 1];  // ǰһ������
                    // ���ǰһ���������깤ʱ�����ڻ�����׼��ʱ�䣬���¼�ӳ�ʱ��
                    if (JFDTime[CurJob][m - 1] > MachReadyTime[m])
                    {
                        DelayTime[PreJob][m] = JFDTime[CurJob][m - 1] - MachReadyTime[m];  // ��¼�ӳ�ʱ��
                    }
                }
            }
            MachReadyTime = JFDTime[CurJob];  // ���»���׼��ʱ��
        }
    }

    return MachReadyTime[this->m_Machines - 1];  // �������һ̨�������깤ʱ��
}


int NOperator::GetDelayTime_Forward_NoAC(vector <vector <int>> FacFamSeq, vector <vector <int>> JobSeqInFam,
                                         vector<int>& FacSpan, vector<vector<int>>& DelayTime)//Forward pass calculation
{
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
        FacSpan[Fac] = GetDelayTime_Forward_InFactory_NoAC(FacFamSeq[Fac], JobSeqInFam, DelayTime);//�õ�ÿ���������깤ʱ��
    return *max_element(FacSpan.begin(), FacSpan.end()); //�õ��ؼ�����
}
int NOperator::GetDelayTime_Forward_InFactory_NoAC(vector <int> FamSeqInFac, vector <vector <int>> JobSeqInFam,  vector<vector<int>>& DelayTime)
{
    vector<vector<int>> JFDTime(this->m_Jobs);
    for (int j = 0; j < JFDTime.size(); j++)
        JFDTime[j].resize(this->m_Machines, 0);

    for (int Fam = 0; Fam < FamSeqInFac.size(); Fam++)
    {
        int CurFam = FamSeqInFac[Fam];  // ��ǰ��
        if (Fam == 0)  // ��һ����
        {
            for (int m = 0; m < this->m_Machines; m++)
                MachReadyTime[m] = this->m_SetupTime[m][CurFam][CurFam];
        }
        else  // �ӵڶ��鵽�����
        {
            int PreFam = FamSeqInFac[Fam - 1];
            for (int m = 0; m < this->m_Machines; m++)
                MachReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
        }

        for (int j = 0; j < JobSeqInFam[CurFam].size(); j++)  // ��ǰ��Ĺ�������
        {
            int CurJob = JobSeqInFam[CurFam][j];  // ��ǰ����
            // ���ڵ�һ�����������㵱ǰ�������깤ʱ��
            JFDTime[CurJob][0] = max(MachReadyTime[0] + this->m_TureJobOpertime[CurJob][0], MachReadyTime[1]);

            for (int m = 1; m < this->m_Machines; m++)
            {
                // �������һ̨������ֱ�Ӽ����깤ʱ��
                if (m == this->m_Machines - 1)
                {
                    JFDTime[CurJob][m] = JFDTime[CurJob][m - 1] + this->m_TureJobOpertime[CurJob][m];
                }
                else
                {
                    // �����������Ͽ������������
                    JFDTime[CurJob][m] = max(JFDTime[CurJob][m - 1] + m_TureJobOpertime[CurJob][m], MachReadyTime[m + 1]);
                }

                // ��¼ DelayTime[PreJob][m] �ڴ˴�
                if (j > 0)  // ������ǵ�һ������
                {
                    int PreJob = JobSeqInFam[CurFam][j - 1];  // ǰһ������
                    // ���ǰһ���������깤ʱ�����ڻ�����׼��ʱ�䣬���¼�ӳ�ʱ��
                    if (JFDTime[CurJob][m - 1] > MachReadyTime[m])
                    {
                        DelayTime[PreJob][m] = JFDTime[CurJob][m - 1] - MachReadyTime[m];  // ��¼�ӳ�ʱ��
                    }
                }
            }
            MachReadyTime = JFDTime[CurJob];  // ���»���׼��ʱ��
        }
    }

    return MachReadyTime[this->m_Machines - 1];  // �������һ̨�������깤ʱ��
}


/**
 *
 * @param FamSeq
 * @param JobSeqinFam
 * @param JFDTime
 * @param FacTT ÿ��������makespan
 * @param scenario
 * @return ���й���makespan���ܺ�
 */
int NOperator::GetSpan_Forward(const vector<vector<int>> &FamSeq, const vector<vector<int>> &JobSeqInFam, vector<int> &FacSpan)
{
    for (int Fac = 0; Fac < FamSeq.size(); Fac++)
        FacSpan[Fac] = GetSpan_Forward_InFactory(FamSeq[Fac], JobSeqInFam);
    return accumulate(FacSpan.begin(), FacSpan.end(), 0);
}

/**
 * ����ÿ��������TT
 * @param FamSeqInFac
 * @param JobSeqInFam
 * @param scenario
 * @return
 */
int NOperator::GetSpan_Forward_InFactory(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam)
{
    for (int Fam = 0; Fam < FamSeqInFac.size(); Fam++)
    {
        int CurFam = FamSeqInFac[Fam];//��ǰ��
        int LastJobinCurFam = JobSeqInFam[CurFam][JobSeqInFam[CurFam].size() - 1];
        if (Fam == 0)//the first group of jobs ��һ����
        {
            for (int m = 0; m < this->m_Machines; m++)
                MachReadyTime[m] = this->m_SetupTime[m][CurFam][CurFam];
        }
        else //�ڶ�������� from the second group of jobs to the end;
        {
            int PreFam = FamSeqInFac[Fam - 1];
            for (int m = 0; m < this->m_Machines; m++)
                MachReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
        }

        for (int j = 0; j < JobSeqInFam[CurFam].size(); j++)//��ǰ��Ĺ������� Scheduling Jobs in CurFam
        {
            int CurJob = JobSeqInFam[CurFam][j];//��ǰ����
            MachReadyTime[0] = max(MachReadyTime[0] + this->m_JobOperPTime[CurJob][0], MachReadyTime[1]);

            for (int m = 1; m < this->m_Machines; m++)
            {
                if (m == this->m_Machines - 1)
                {
                    MachReadyTime[m] = MachReadyTime[m - 1] + this->m_JobOperPTime[CurJob][m];
                }
                else
                {
                    MachReadyTime[m] = max(MachReadyTime[m - 1] + m_JobOperPTime[CurJob][m],
                                           MachReadyTime[m + 1]);
                }
            }
        }
    }
    return MachReadyTime[this->m_Machines - 1];
}

/**
 * �ҵ�makespan��С��λ�ò����飨��ʹ�ÿ�����
 * @param FamSeq
 * @param JobSeqinFam
 * @param InsFam
 * @param bestFac
 * @param bestPos
 * @return
 */
int NOperator::FindBestPosToInsertFam_NoAC(const vector<vector<int>> &FamSeq, const vector<vector<int>> &JobSeqinFam,
                                           int InsFam, int &bestFac, int &bestPos)
{
    int minSpan = INT_MAX;
    for (int Fac = 0; Fac < FamSeq.size(); Fac++)
    {
        int Pos;
        int Span = FindBestPosToInsertFam_NoAC_InFactory(FamSeq[Fac], JobSeqinFam, InsFam, Pos);
        if (Span < minSpan)
        {
            minSpan = Span;
            bestFac = Fac;
            bestPos = Pos;
        }
    }
    return minSpan;
}

/**
 * �ڹ������ҵ�makespan��С�ĵ�λ�ò����飨��ʹ�ÿ�����
 * @param FamSeqInFac
 * @param JobSeqinFam
 * @param InsFam
 * @param bestPos
 * @return
 */
int
NOperator::FindBestPosToInsertFam_NoAC_InFactory(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqinFam,
                                                 int InsFam, int &bestPos)
{
    int minSpan = INT_MAX;
    //vector<int> MachSpan(this->m_Machines);
    //vector<int> MachReadyTime(this->m_Machines, 0);  //��ʼ����Ϊ0

    vector<vector<int>> JFDTime(this->m_Jobs);
    for (int j = 0; j < JFDTime.size(); j++)
    {
        JFDTime[j].resize(this->m_Machines, 0);
    }
    vector<int> tempFamSeqInFac = FamSeqInFac;
    for (int Pos = 0; Pos <= FamSeqInFac.size(); Pos++)  //���ڹ����е�����λ��
    {
        tempFamSeqInFac.insert(begin(tempFamSeqInFac) + Pos, InsFam);

        int Span = GetJFDTime_Forward_InFactory(tempFamSeqInFac, JobSeqinFam, JFDTime);
        if (Span < minSpan)
        {
            minSpan = Span;  //ȡ����λ����span��С��
            bestPos = Pos;   //��¼λ��
        }
    }
    return minSpan;  //������С��span
}


/**
 * ���������뵽makespan��С��λ�ã���ʹ�ÿ�����
 * @param FamSeqInFac
 * @param JobSeqinFam
 * @param Fam
 * @param InsJob
 * @param bestPos
 * @return
 */
int
NOperator::FindBestPosToInsertJob_NoAC(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqinFam, int Fam,
                                       int InsJob, int &bestPos)
{
    vector<vector<int>> JFDTime(m_Jobs, vector<int>(m_Machines, 0));

    vector<vector<int>> tempJobSeqinFam = JobSeqinFam;
    int minSpan = INT_MAX;
    for (int Pos = 0; Pos <= JobSeqinFam[Fam].size(); Pos++)
    {
        tempJobSeqinFam[Fam].insert(begin(tempJobSeqinFam[Fam]) + Pos, InsJob);

        //int Span = *max_element(MachSpan.begin(), MachSpan.end());

        int Span = GetJFDTime_Forward_InFactory(FamSeqInFac, tempJobSeqinFam, JFDTime);
        if (Span < minSpan)
        {
            minSpan = Span;
            bestPos = Pos;   //��¼��õ�λ��
        }
        tempJobSeqinFam[Fam].erase(begin(tempJobSeqinFam[Fam]) + Pos);
    }
    return minSpan;
}

/**
 * ���������й������ҵ�makespan��С��λ�ã������룩
 * ��ʹ�ÿ���
 * @param FamSeq
 * @param JobSeqinFam
 * @param FacTT_AllScenarios
 * @param JFDTime
 * @param InsFam
 * @param bestFac
 * @param bestPos
 * @return
 */
int NOperator::FindBestPosToInsertFam_NoAC_Span(vector<vector<int>> &FamSeq, const vector<vector<int>> &JobSeqinFam,int InsFam, int &bestFac, int &bestPos)

{
    int minSpan = INT_MAX;
    for (int Fac = 0; Fac < FamSeq.size(); Fac++)
    {
        int Pos;
        int Span = FindBestPosToInsertFam_NoAC_InFactory_Span(FamSeq, JobSeqinFam, InsFam, Pos, Fac);
        if (Span < minSpan)
        {
            minSpan = Span;
            bestFac = Fac;
            bestPos = Pos;
        }
    }
    return minSpan;
}

/**
 * ���������й������ҵ�TT��С��λ�ã������룩
 * ��ʹ�ÿ���
 * @param FamSeq
 * @param JobSeqinFam
 * @param FacTT_AllScenarios
 * @param JFDTime
 * @param InsFam
 * @param bestFac
 * @param bestPos
 * @return
 */
int NOperator::FindBestPosToInsertFam_NoAC_InFactory_Span(vector<vector<int>> &FamSeq,
                                                          const vector<vector<int>> &JobSeqInFam, int InsFam,
                                                          int &bestPos, int factory)

{
    int minSpan = INT_MAX;
    int Span;

    for (int Pos = 0; Pos <= FamSeq[factory].size(); ++Pos)
    {
        // Insert the group at position Pos
        FamSeq[factory].insert(FamSeq[factory].begin() + Pos, InsFam);

        // Calculate the total TT after insertion
        Span = GetSpan_Forward_InFactory(FamSeq[factory], JobSeqInFam);

        // Remove the group after calculation to restore the sequence
        FamSeq[factory].erase(FamSeq[factory].begin() + Pos);

        if (Span < minSpan)
        {
            minSpan = Span;
            bestPos = Pos;
        }
    }
    return minSpan;
}

/**
 * �������ҵ�makespan��С��λ�ã������룩
 * ��ʹ�ÿ���
 * @param FamSeq
 * @param JobSeqinFam
 * @param JFDTime
 * @param Fam
 * @param InsJob
 * @param bestPos
 * @return
 */
int NOperator::FindBestPosToInsertJob_NoAC_Span(const vector<vector<int>> &FamSeq, vector<vector<int>> &JobSeqinFam, vector<int> &FacTT, int Fam, int InsJob,
                                                int &bestPos)
{
    //�ҵ�Fam���ڵĹ���
    int Fac = 0;
    for (int i = 0; i < FamSeq.size(); ++i)
    {
        auto it = find(FamSeq[i].begin(), FamSeq[i].end(), Fam);
        if (it != FamSeq[i].end())
        {
            Fac = i;
            break;
        }
    }
    int minSpan = INT_MAX;
    vector<int> temp_FacTT= FacTT;
    for (int Pos = 0; Pos <= JobSeqinFam[Fam].size(); ++Pos)
    {
        JobSeqinFam[Fam].insert(begin(JobSeqinFam[Fam]) + Pos, InsJob);

        int Span = GetSpan_Forward_InFactory(FamSeq[Fac], JobSeqinFam);

        JobSeqinFam[Fam].erase(begin(JobSeqinFam[Fam]) + Pos);


        // If the new TT is smaller, update bestPos
        if (Span < minSpan)
        {
            minSpan =Span;
            bestPos = Pos;
        }
    }

    return minSpan;
}



float NOperator::GetTECForPerFac(const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam)
{
    double TEC = 0;
    vector<vector<int>> JFDTime(this->m_Jobs);
    for (int j = 0; j < JFDTime.size(); j++)
        JFDTime[j].resize(this->m_Machines, 0);
    vector<float> FacEC(m_Factories, 0);
    for (int Fac = 0; Fac < m_Factories; Fac++)
    {
        FacEC[Fac] = GetPerFacEC(FamSeqInFac, JobSeqInFam, JFDTime);//�õ�ÿ��������EC
        TEC += FacEC[Fac];
    }
    return TEC;
}

float NOperator::GetTECForAllFac(const vector<vector<int>> &FamSeq, const vector<vector<int>>& JobSeqInFam)
{
    double TEC = 0;
    vector<vector<int>> JFDTime(this->m_Jobs);
    for (int j = 0; j < JFDTime.size(); j++)
        JFDTime[j].resize(this->m_Machines, 0);
    vector<float> FacEC(m_Factories, 0);
    for (int Fac = 0; Fac < m_Factories; Fac++)
    {
        // ����ù���������Ϊ�գ�������
        if (FamSeq[Fac].empty()) {
            continue;
        }
        FacEC[Fac] = GetPerFacEC(FamSeq[Fac], JobSeqInFam, JFDTime);//�õ�ÿ��������EC
        TEC += FacEC[Fac];
    }
    return TEC;
}

float NOperator::FindBestPosToInsertFamForAllFac_TEC(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam,  int InsFam, int& BestFac, int& BestPos)
{
    float minTEC = INT_MAX;
    //����δ����֮ǰ��������EC
    vector <float> FacEC(FacFamSeq.size(), 0.0);
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {   // ����ù���������Ϊ�գ�������
        if (FacFamSeq[Fac].empty()) {
            continue;
        }
        FacEC[Fac] = GetTECForPerFac(FacFamSeq[Fac], JobSeqInFam);
        //cout << "δ����ʱ����" << Fac << "��EC��" << FacEC[Fac] << endl;
    }

    vector <float> TempFacEC(FacFamSeq.size(), 0.0);
    TempFacEC = FacEC;

    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        int Pos = -1;
        TempFacEC[Fac] = this->FindBestPosToInsertFamForPerFac_TEC(FacFamSeq[Fac], JobSeqInFam, InsFam, Pos);

        float TEC = accumulate(TempFacEC.begin(), TempFacEC.end(), 0);

        if (TEC < minTEC)
        {
            minTEC = TEC;
            BestFac = Fac;
            BestPos = Pos;
        }
        TempFacEC[Fac] = FacEC[Fac];
    }
    return minTEC;
}

double NOperator::FindBestPosToInsertFamForAllFacs_MC(const vector<vector<int>> &FacFamSeq,
                                                      const vector<vector<int>> &JobSeqInFam,
                                                      int InsFam, int &BestFac, int &BestPos) {
    int minMC = INT_MAX;
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++) {
        int Pos;
        int MC = FindBestPosToInsertFam_InFactory_MC(Fac, FacFamSeq, FacFamSeq[Fac], JobSeqInFam, InsFam,Pos);
        if (MC < minMC) {
            minMC = MC;
            BestFac = Fac;
            BestPos = Pos;
        }
    }
    return minMC;
}


double NOperator::FindBestPosToInsertFam_InFactory_MC(int Fac , const vector<vector<int>> &FacFamSeq,const vector<int> &NewFamSeqInFac,
                                                      const vector<vector<int>> &JobSeqInFam,int InsFam,int &BestPos) {
    int minMC = INT_MAX;
    int MC = 0;
    for (int Pos = 0; Pos <= NewFamSeqInFac.size(); ++Pos) {
        auto re = GetMCForPerFacAfterInsertFam(Fac,FacFamSeq, NewFamSeqInFac, JobSeqInFam, InsFam, Pos);
        MC = re;
        if (MC < minMC) {
            minMC = MC;
            BestPos = Pos;
        }
    }

    return minMC;
}

double
NOperator::GetMCForPerFacAfterInsertFam(int Fac , const vector<vector<int>> &FacFamSeq,const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,int InsFam, int Pos) {

    //����δ����֮ǰ��������span
    vector <int> FacSpan(FacFamSeq.size(), 0);
    vector<vector<int>> TempFacFamSeq;
    TempFacFamSeq.clear();
    TempFacFamSeq.resize(FacFamSeq.size());
    TempFacFamSeq = FacFamSeq;
    TempFacFamSeq[Fac].insert(TempFacFamSeq[Fac].begin() + Pos, InsFam);
    int Span = GetSpan_Forward_InFactory(TempFacFamSeq[Fac], JobSeqInFam);

    for (int f = 0; f < FacFamSeq.size(); f++) {
        FacSpan[f]= GetSpan_Forward_InFactory(TempFacFamSeq[Fac], JobSeqInFam);
    }
    //�������󹤳�Fac��span
    FacSpan[Fac]= Span;
    int Makespan = *max_element(FacSpan.begin(), FacSpan.end());
    vector <int> FacEC(FacFamSeq.size(), 0);
    for (int fac = 0; fac < TempFacFamSeq.size(); fac++)
    {
        //�����ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
            //int CurFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
            int Fam = TempFacFamSeq[fac][0];
            for (int g = 0; g < TempFacFamSeq[fac].size(); g++)
            {
                Fam = TempFacFamSeq[fac][g];

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
    }

    float TEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);

    double MC = 0.5f * Makespan + 0.5f * TEC ;
    return MC;
}


float NOperator::FindBestPosToInsertFamForPerFac_TEC(const vector<int>& NewFamSeqInFac, const vector<vector<int>>& JobSeqInFam, int InsFam, int& BestPos)
{
    float minFacEC = INT_MAX;
    vector <int> TempFamSeqInFac = NewFamSeqInFac;
    for (int Pos = 0; Pos <= NewFamSeqInFac.size(); Pos++)
    {
        TempFamSeqInFac.insert(TempFamSeqInFac.begin() + Pos, InsFam);
        float FacEC = GetTECForPerFac(TempFamSeqInFac, JobSeqInFam);
        if (FacEC < minFacEC)
        {
            minFacEC = FacEC;
            BestPos = Pos;
        }
        TempFamSeqInFac.erase(TempFamSeqInFac.begin() + Pos);
    }
    return minFacEC;
}


int NOperator::FindBestPosToInsertFam(vector <vector <int>> FamSeq, vector <vector <int>> JobSeqinFam, vector<vector<int>> JFDTime, vector<vector<int>> JBDTime, int InsFam, int& bestFac, int& bestPos)
{
    int minSpan = INT_MAX;
    for (int Fac = 0; Fac < FamSeq.size(); Fac++)
    {
        int Pos;
        int Span = FindBestPosToInsertFam_InFactory(FamSeq[Fac], JobSeqinFam, JFDTime, JBDTime, InsFam, Pos);
        if (Span < minSpan)
        {
            minSpan = Span;
            bestFac = Fac;
            bestPos = Pos;
        }
    }
    return minSpan;
}

/**
 * �ڹ������ҵ���õ�λ�ò�����
 * @param FamSeqInFac
 * @param JobSeqinFam
 * @param JFDTime
 * @param JBDTime
 * @param InsFam
 * @param bestPos
 * @return
 */
//gai
int NOperator::FindBestPosToInsertFam_InFactory(vector <int> FamSeqInFac, vector <vector <int>> JobSeqinFam, vector<vector<int>> JFDTime,vector<vector<int>> JBDTime,  int InsFam, int& bestPos)
{
    int minSpan = INT_MAX;
    for (int Pos = 0; Pos <= FamSeqInFac.size(); Pos++)  //���ڹ����е�����λ��
    {
        if (Pos == 0)//��һ��λ�ò���Insert at Postion 0
        {
            for (int m = 0; m < this->m_Machines; m++)
            {
                MachReadyTime[m] = this->m_SetupTime[m][InsFam][InsFam];  //MachReadyTime ��Ϊ׼��ʱ��
            }
        }
        else // ������λ�ò��� Insert at other Postions
        {
            int PreFam = FamSeqInFac[Pos - 1];  //ǰһ����
            int LastJobinPreFam = JobSeqinFam[PreFam][JobSeqinFam[PreFam].size() - 1];  //ǰһ��������һ������
            for (int m = 0; m < this->m_Machines; m++)
            {
                MachReadyTime[m] = JFDTime[LastJobinPreFam][m] + this->m_SetupTime[m][PreFam][InsFam];  //ǰһ�������һ���������깤ʱ��+׼��ʱ��
            }
        }

        // ���﹤������ Scheduling jobs in InsFam
        for (int CurJob: JobSeqinFam[InsFam])
        {
            //��ǰ����
            MachReadyTime[0] = max(MachReadyTime[0] + this->m_TureJobOpertime[CurJob][0], MachReadyTime[1]);   //�ڵ�һ��������,����׼��ʱ��=max(��ǰ�����뿪ʱ�䣬��һ������׼��ʱ��)
            for (int m = 1; m < this->m_Machines; m++)
            {
                if (m == this->m_Machines - 1)
                {
                    MachReadyTime[m] = MachReadyTime[m - 1] + this->m_TureJobOpertime[CurJob][m];//���һ�������ϵĹ����뿪ʱ������깤ʱ��
                } else
                {
                    MachReadyTime[m] = max(MachReadyTime[m - 1] + this->m_TureJobOpertime[CurJob][m], MachReadyTime[m + 1]);//�뿪ʱ��=max(����׼��ʱ��+����ʱ��,��һ������׼��ʱ��)
                }
            }
        }

        int span = INT_MIN;
        if (Pos < FamSeqInFac.size()) //�������һ��λ��Compute makespan;
        {
            int NextFam = FamSeqInFac[Pos];  //��һ����
            int FirstJobinNextFam = JobSeqinFam[NextFam][0];   //��һ����ĵ�һ������
            for (int m = 0; m < this->m_Machines; m++)
            {
                span = max(span, MachReadyTime[m] + this->m_SetupTime[m][InsFam][NextFam] + JBDTime[FirstJobinNextFam][m]);  //�깤ʱ��=�ӿ�ͷ����ǰ���ϵ��뿪ʱ��+��ǰ������һ�����׼��ʱ��+������㵽��ǰ����깤ʱ��
            }
        } else
        {
            span = *max_element(MachReadyTime.begin(), MachReadyTime.end());   //spanȡ���ֵ
        }
        if (span < minSpan)
        {
            minSpan = span;  //ȡ����λ����span��С��
            bestPos = Pos;   //��¼λ��
        }
    }
    return minSpan;  //������С��span
}
/**
 * ���ݹ����������뿪ʱ����㵥��������MAKESPAN
 * @param FamSeqInFac
 * @param JobSeqInFam
 * @param JFDTime
 * @return
 */
int NOperator::GetSpanPerFacByJFD(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                  const vector<vector<int>> &JFDTime)
{
    int Span = 0;
    for (int f = 0; f < FamSeqInFac.size(); ++f)
    {
        int Fam = FamSeqInFac[f];
        int LastJobInFam = JobSeqInFam[Fam][JobSeqInFam[Fam].size() - 1];//������һ������
        Span = JFDTime[LastJobInFam][this->m_Machines - 1]  ;
    }
    return Span;
}

/**
 * ���ݹ�����ǰ���뿪ʱ��������ܺ�TEC
 * @param FamSeqInFac
 * @param JobSeqInFam
 * @param JFDTime
 * @return
 */
float NOperator::GetTECForAllFacByJFD(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam,
                                      const vector<vector<int>>& JFDTime)
{
    float TEC = 0;
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        float TempTEC = GetTECForPerFacByJFD(FacFamSeq[Fac], JobSeqInFam, JFDTime);
        TEC += TempTEC;
    }

    return TEC;
}
/**
 * ���ݹ����������뿪ʱ��������ܺ�TEC
 * @param FamSeqInFac
 * @param JobSeqInFam
 * @param JFDTime
 * @return
**/
float NOperator::GetTECForPerFacByJFD(const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam,const vector<vector<int>>& JFDTime)
{
    int LastFam = *FamSeqInFac.rbegin();
    int LastJobInLstFam = *JobSeqInFam[LastFam].rbegin(); //������һ������
    int makespan = JFDTime[LastJobInLstFam][this->m_Machines - 1];

    //�����ܺ�
    float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

    for (int m = 0; m < m_Machines; m++)
    {
        int preFam = -1;
        //int CurFam = -1;
        int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
        int Fam = FamSeqInFac[0];
        for (int g = 0; g < FamSeqInFac.size(); g++)
        {
            Fam = FamSeqInFac[g];

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
        Tidletime += (makespan - Tptime - Tsetuptime);
        TIdleTimeEC += Tidletime * UnitIdleEC;
    }
    float FacEC = TSetupTimeEC + TIdleTimeEC + TpTimeEC;
    return FacEC;
}

/**
 * �Դ����ȹ�����InsFam�����й���FacFamSeq��Ѱ�����λ��(Makespan)
 * @param FacFamSeq
 * @param JobSeqInFam
 * @param JFDTime
 * @param InsFam
 * @param BestFac
 * @param BestPos
 * @return
 */
int NOperator::FindBestPosToInsertFamForAllFac_Makespan(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam,
                                                        const vector<vector<int>>& JFDTime, const vector<vector<int>>& JBDTime, int InsFam, int& BestFac, int& BestPos)
{
    int minMakespan = INT_MAX;

    //����δ����֮ǰ��������span
    vector <int> FacSpan(FacFamSeq.size(), 0);
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        int LastFam = *FacFamSeq[Fac].rbegin();
        int LastJobInLstFam = *JobSeqInFam[LastFam].rbegin(); //������һ������
        FacSpan[Fac] = JFDTime[LastJobInLstFam][this->m_Machines - 1];
    }
    vector <int> TempFacSpan(FacFamSeq.size(), 0);
    TempFacSpan = FacSpan;

    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        int Pos = -1;
        TempFacSpan[Fac] = this->FindBestPosToInsertFamForPerFac_Makespan(FacFamSeq[Fac], JobSeqInFam, JFDTime, JBDTime, InsFam, Pos);
        int Makespan = *max_element(TempFacSpan.begin(), TempFacSpan.end());

        if (Makespan < minMakespan)
        {
            minMakespan = Makespan;
            BestFac = Fac;
            BestPos = Pos;
        }
        TempFacSpan[Fac] = FacSpan[Fac];
    }
    return minMakespan;
}

/**
 * �Դ����ȹ�����InsFam�ڵ�ǰ����NewFamSeqInFac��Ѱ�����λ��(Makespan)
 * @param FacFamSeq
 * @param JobSeqInFam
 * @param JFDTime
 * @param InsFam
 * @param BestFac
 * @param BestPos
 * @return
 */

int NOperator::FindBestPosToInsertFamForPerFac_Makespan(const vector<int>& NewFamSeqInFac, const vector<vector<int>>& JobSeqInFam,
                                                        const vector<vector<int>>& JFDTime, const vector<vector<int>>& JBDTime, int InsFam, int& BestPos)
{
    int minFacMakespan = INT_MAX;
    for (int Pos = 0; Pos <= NewFamSeqInFac.size(); Pos++)
    {
        int FacMakespan = GetSpanForPerFacAfterInsertFam(NewFamSeqInFac, JobSeqInFam, JFDTime, JBDTime, InsFam, Pos);
        if (FacMakespan < minFacMakespan)
        {
            minFacMakespan = FacMakespan;
            BestPos = Pos;
        }
    }
    return minFacMakespan;
}

/**
 * �¿���
 * �Դ����ȹ�����InsFam�����й���FacFamSeq��Ѱ�����λ��(TEC)
 * @param FacFamSeq
 * @param JobSeqInFam
 * @param JFDTime
 * @param InsFam
 * @param BestFac
 * @param BestPos
 * @return
 */
float NOperator::FindBestPosToInsertFamForAllFac_TEC(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam,
                                                     const vector<vector<int>>& JFDTime, const vector<vector<int>>& JBDTime, int InsFam, int& BestFac, int& BestPos)
{
    float minTEC = INT_MAX;

    //����δ����֮ǰ��������EC
    vector <float> FacEC(FacFamSeq.size(), 0.0);
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        FacEC[Fac] = GetTECForPerFacByJFD(FacFamSeq[Fac], JobSeqInFam, JFDTime);

    }
    vector <float> TempFacEC(FacFamSeq.size(), 0.0);
    TempFacEC = FacEC;

    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        int Pos = -1;
        TempFacEC[Fac] = this->FindBestPosToInsertFamForPerFac_TEC(FacFamSeq[Fac], JobSeqInFam, JFDTime, JBDTime, InsFam, Pos);

        float TEC = accumulate(TempFacEC.begin(), TempFacEC.end(), 0.0);

        if (TEC < minTEC)
        {
            minTEC = TEC;
            BestFac = Fac;
            BestPos = Pos;
        }
        TempFacEC[Fac] = FacEC[Fac];
    }
    return minTEC;
}

/**
 * �Դ����ȹ�����InsFam�ڵ�ǰ����NewFamSeqInFac��Ѱ�����λ��(TEC)
 * @param FacFamSeq
 * @param JobSeqInFam
 * @param JFDTime
 * @param InsFam
 * @param BestFac
 * @param BestPos
 * @return
 */
float NOperator::FindBestPosToInsertFamForPerFac_TEC(const vector<int>& NewFamSeqInFac, const vector<vector<int>>& JobSeqInFam,
                                                     const vector<vector<int>>& JFDTime, const vector<vector<int>>& JBDTime, int InsFam, int& BestPos)
{
    float minFacEC = INT_MAX;

    for (int Pos = 0; Pos <= NewFamSeqInFac.size(); Pos++)
    {
        float FacEC = GetECForPerFacAfterInsertFam(NewFamSeqInFac, JobSeqInFam, JFDTime, JBDTime, InsFam, Pos);

        if (FacEC < minFacEC)
        {
            minFacEC = FacEC;
            BestPos = Pos;
        }

    }
    return minFacEC;
}

//��
void NOperator::GetMSandTECForPerandToalFac(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam,
                                            vector<int>& FacSpan, vector<float>& FacEC, int& Makespan, float& TotalEC)
{

    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        vector<vector<int>> JFDTime(this->m_Jobs);
        for (int j = 0; j < JFDTime.size(); j++)
            JFDTime[j].resize(this->m_Machines, 0);
        this->GetJFDTime_Forward(FacFamSeq, JobSeqInFam, JFDTime, FacSpan);

//        int LastFam = *FacFamSeq[Fac].rbegin();
//        int LastJobInLstFam = *JobSeqInFam[LastFam].rbegin(); //������һ������
//        FacSpan[Fac] = JFDTime[LastJobInLstFam][this->m_Machines - 1];
        if (!FacFamSeq[Fac].empty())
        {
            int LastFam = *FacFamSeq[Fac].rbegin();
            if (!JobSeqInFam[LastFam].empty())
            {
                int LastJobInLstFam = *JobSeqInFam[LastFam].rbegin(); // ������һ������
                FacSpan[Fac] = JFDTime[LastJobInLstFam][this->m_Machines - 1];
            }
            else
            {
                FacSpan[Fac] = 0; // �������޹������趨�ʵ���Ĭ��ֵ
            }
        }
        else
        {
            FacSpan[Fac] = 0; // ���ù������飬�趨�ʵ���Ĭ��ֵ
        }
        //�����ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
          //  int Fam = FacFamSeq[Fac][0];
            // �жϵ�ǰ�������������Ƿ�Ϊ��
            if (FacFamSeq[Fac].empty()) {
                continue; // ֱ��������ǰ����
            }
            int Fam = FacFamSeq[Fac][0]; // ���ڿ��԰�ȫ����
            for (int g = 0; g < FacFamSeq[Fac].size(); g++)
            {
                Fam = FacFamSeq[Fac][g];

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
            Tidletime += (FacSpan[Fac] - Tptime - Tsetuptime);
            TIdleTimeEC += Tidletime * UnitIdleEC;
        }
        FacEC[Fac] = TSetupTimeEC + TIdleTimeEC + TpTimeEC;
    }
    Makespan = *max_element(FacSpan.begin(), FacSpan.end()); //�õ��ؼ�����
    TotalEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);

}

void NOperator::GetMSandTECForPerandToalFacByJFD(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam,const vector<vector<int>>& JFDTime,
                                                 vector<int>& FacSpan, vector<float>& FacEC, int& Makespan, float& TotalEC)
{

    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        int LastFam = *FacFamSeq[Fac].rbegin();
        int LastJobInLstFam = *JobSeqInFam[LastFam].rbegin(); //������һ������
        FacSpan[Fac] = JFDTime[LastJobInLstFam][this->m_Machines - 1];

        //�����ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
            int Fam = FacFamSeq[Fac][0];
            for (int g = 0; g < FacFamSeq[Fac].size(); g++)
            {
                Fam = FacFamSeq[Fac][g];

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
            Tidletime += (FacSpan[Fac] - Tptime - Tsetuptime);
            TIdleTimeEC += Tidletime * UnitIdleEC;
        }
        FacEC[Fac] = TSetupTimeEC + TIdleTimeEC + TpTimeEC;
    }
    Makespan = *max_element(FacSpan.begin(), FacSpan.end()); //�õ��ؼ�����
    TotalEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);

}

//Pareto_relation
void NOperator::Pareto_relation(vector<Individual>& Population)
{
    int pflag, pflag1;
    for (int i = 0; i < Population.size(); i++)
    {
        Population[i].flag = 0;
        for (int j = 0; j < Population.size(); j++)
        {
            pflag = 1;
            pflag1 = 0;
            if (Population[i].MS > Population[j].MS)
                pflag = 0;
            if (Population[i].MS < Population[j].MS)
                pflag1 = 1;
            if (Population[i].TEC > Population[j].TEC)
                pflag = 0;
            if (Population[i].TEC < Population[j].TEC)
                pflag1 = 1;


            if (pflag == 1 && pflag1 == 1)
                Population[j].pareto_rel[i] = 1;	//iռ��j
            else
                Population[j].pareto_rel[i] = 0;
        }

    }
}

void NOperator::Speed_mutation_new(vector<Individual> &Population, vector<Individual> &Temp) {
    vector<Individual> tempCCEAPopulation;
    tempCCEAPopulation.clear();
    int Orgsize = Population.size();

    for (int PS = 0; PS < Orgsize; PS++)
    {
        //���� ����λʱ��ӹ��ܺ�
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                Temp[PS].m_SpeedVector[j][i] = rand() % 3;
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[Temp[PS].m_SpeedVector[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[Temp[PS].m_SpeedVector[j][i]] * m_Speed[Temp[PS].m_SpeedVector[j][i]];
            }
        }


        vector<int> FacSpan(m_Factories, 0);
        vector<float> FacEC(m_Factories, 0);
        int Makespan = 0;
        float TotalEC = 0;

        GetMSandTECForPerandToalFac(Temp[PS].m_FacFamSeqArray, Temp[PS].m_JobSeqInFamArray, FacSpan, FacEC, Makespan, TotalEC);

        Temp[PS].MS = Makespan;
        Temp[PS].TEC = TotalEC;
        if (((Temp[PS].MS < Population[PS].MS) && (Temp[PS].TEC < Population[PS].TEC)) || ((Temp[PS].MS < Population[PS].MS) && (Temp[PS].TEC == Population[PS].TEC)) || ((Temp[PS].MS == Population[PS].MS) && (Temp[PS].TEC < Population[PS].TEC)))
        {
            Population[PS].m_SpeedVector = Temp[PS].m_SpeedVector;
            Population[PS].MS = Temp[PS].MS;
            Population[PS].TEC = Temp[PS].TEC;

        }

        tempCCEAPopulation.push_back(Population[PS]);

        m_SpeedMatrix = tempCCEAPopulation[PS].m_SpeedVector;
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }

        int tempMakespan = 0;
        float tempTotalEC = 0.0;
        vector<int> tempFacSpan(m_Factories, 0);
        vector<float> tempFacEC(m_Factories, 0);

        GetMSandTECForPerandToalFac(Population[PS].m_FacFamSeqArray, Population[PS].m_JobSeqInFamArray, tempFacSpan, tempFacEC, tempMakespan, tempTotalEC);

        //���ܲ���1
        vector<vector<int>> DelayTime;
        DelayTime.clear();
        DelayTime.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
        vector<vector<int>> JFDTime(this->m_Jobs);
        for (int j = 0; j < JFDTime.size(); j++)
            JFDTime[j].resize(this->m_Machines, 0);
        GetDelayTime_Forward(tempCCEAPopulation[PS].m_FacFamSeqArray, tempCCEAPopulation[PS].m_JobSeqInFamArray, JFDTime, FacSpan, DelayTime);

        vector<vector<int>> JBDTime(this->m_Jobs);
        for (int j = 0; j < JBDTime.size(); j++)
            JBDTime[j].resize(this->m_Machines, 0);

        GetJBDTime_Backward(tempCCEAPopulation[PS].m_FacFamSeqArray, tempCCEAPopulation[PS].m_JobSeqInFamArray, JBDTime, FacSpan);

        bool sign = false;

        // ��ʼ�����в����Ĺؼ�·�����Ϊ false
        vector<vector<bool>> criticalPathFlag(m_Jobs, vector<bool>(m_Machines, false));

        // �����һ�������Ĳ�����ʼ��������һ������Ϊ�ؼ�·��
        for (int j = 0; j < m_Jobs; j++) {
            criticalPathFlag[j][m_Machines - 1] = true; // ���һ���������Ϊ�ؼ�·��
        }

        // ������֤ÿ�������Ƿ��ڹؼ�·����
        for (int j = 0; j < m_Jobs; j++) {
            for (int i = m_Machines - 2; i >= 0; i--) {
                // �жϲ����Ƿ��������ģ�����������ģ���Ϊ�ؼ�·����һ����
                if (JFDTime[j][i] == JBDTime[j][i + 1] || JFDTime[j][i] == JBDTime[j][i - 1]) {
                    criticalPathFlag[j][i] = true; // ��ǰ�����ǹؼ�·����һ����
                }
                // �ж��Ƿ��������Լ��
                if (JFDTime[j][i] - this->m_TureJobOpertime[j][i] == JBDTime[j][i + 1]) {
                    criticalPathFlag[j][i] = true; // ��������Լ������ǰ�����ǹؼ�·��
                }
            }
        }

        for (int j = 0; j < m_Jobs; j++) {
            for (int i = 1; i < m_Machines; i++) {
                // �жϵ�ǰ�����Ƿ��ڹؼ�·���ϣ������Ƿ��п���ʱ��
                if ((DelayTime[j][i] > 0) && (i < m_Machines - 1)) {
                    // ���㵱ǰ������ʵ�ʼӹ�ʱ��
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[tempCCEAPopulation[PS].m_SpeedVector[j][i]]);
                    m_TureJobOpertime[j][i + 1] = static_cast<int>(m_JobOperPTime[j][i + 1] / m_Speed[tempCCEAPopulation[PS].m_SpeedVector[j][i + 1]]);

                    // �����ǰ�������ڹؼ�·�������п���ʱ�䣬�����
                    if (criticalPathFlag[j][i] == 0) {
                        if ((JFDTime[j][i] < (JFDTime[j][i + 1] - m_TureJobOpertime[j][i + 1]))) {
                            int Speedlevel = tempCCEAPopulation[PS].m_SpeedVector[j][i] - 1;
                            for (int level = Speedlevel; level > 0; level--) {
                                if ((static_cast<int>(m_JobOperPTime[j][i] / m_Speed[level]) - m_TureJobOpertime[j][i]) < DelayTime[j][i]) {
                                    tempCCEAPopulation[PS].m_SpeedVector[j][i] = level;
                                    sign = true;
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                        // �����ǰ�����ڹؼ�·������û�п���ʱ�䣬������
                    else {
                        if ((JFDTime[j][i] < (JFDTime[j][i + 1] - m_TureJobOpertime[j][i + 1]))) {
                            int Speedlevel = tempCCEAPopulation[PS].m_SpeedVector[j][i] - 1;
                            for (int level = Speedlevel; level > 0; level--) {
                                if ((static_cast<int>(m_JobOperPTime[j][i] / m_Speed[level]) - m_TureJobOpertime[j][i]) < DelayTime[j][i]) {
                                    if ((JFDTime[j][i] + (static_cast<int>(m_JobOperPTime[j][i] / m_Speed[level]) - m_TureJobOpertime[j][i])) < (JFDTime[j][i + 1] - m_TureJobOpertime[j][i + 1])) {
                                        tempCCEAPopulation[PS].m_SpeedVector[j][i] = level;
                                        sign = true;
                                    } else {
                                        break;
                                    }
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                }
                    // �������һ̨�������������
                else if ((DelayTime[j][i] > 0) && (i == m_Machines - 1)) {
                    int Speedlevel = tempCCEAPopulation[PS].m_SpeedVector[j][i] - 1;
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[tempCCEAPopulation[PS].m_SpeedVector[j][i]]);

                    // �ж��Ƿ��ڹؼ�·����
                    if (criticalPathFlag[j][i] == 1) {
                        // ����ڹؼ�·���ϲ���û�п���ʱ�䣬ִ������
                        for (int level = Speedlevel; level > 0; level--) {
                            if ((static_cast<int>(m_JobOperPTime[j][i] / m_Speed[level]) - m_TureJobOpertime[j][i]) < DelayTime[j][i]) {
                                tempCCEAPopulation[PS].m_SpeedVector[j][i] = level;
                                sign = true;
                            } else {
                                break;
                            }
                        }
                    }
                    else {
                        // ������ڹؼ�·���ϲ����п���ʱ�䣬ִ�м���
                        for (int level = Speedlevel; level > 0; level--) {
                            if ((static_cast<int>(m_JobOperPTime[j][i] / m_Speed[level]) - m_TureJobOpertime[j][i]) < DelayTime[j][i]) {
                                tempCCEAPopulation[PS].m_SpeedVector[j][i] = level;
                                sign = true;
                            } else {
                                break;
                            }
                        }
                    }
                }
            }  // end machine
        }  // end job


        if (sign == true)
        {

            m_SpeedMatrix = tempCCEAPopulation[PS].m_SpeedVector;

            //�õ������Ĵ���ʱ�� ����λʱ��ӹ��ܺ�
            for (int j = 0; j < m_Jobs; j++)
            {
                for (int i = 0; i < m_Machines; i++)
                {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }


            int Makespan1 = 0;
            float TotalEC1 = 0.0;
            vector<int> FacSpan1(m_Factories, 0);
            vector<float> FacEC1(m_Factories, 0);

            GetMSandTECForPerandToalFac(tempCCEAPopulation[PS].m_FacFamSeqArray, tempCCEAPopulation[PS].m_JobSeqInFamArray, FacSpan1, FacEC1, Makespan1, TotalEC1);

            tempCCEAPopulation[PS].MS = Makespan1;
            tempCCEAPopulation[PS].TEC = TotalEC1;

            //���
            //CheckSol(tempCCEAPopulation[PS].m_FacFamSeqArray, tempCCEAPopulation[PS].m_JobSeqInFamArray, tempCCEAPopulation[PS].MS);
            //CheckSolTEC(tempCCEAPopulation[PS].m_FacFamSeqArray, tempCCEAPopulation[PS].m_JobSeqInFamArray, tempCCEAPopulation[PS].TEC);

            //�ж�
            if ((Population[PS].MS == tempCCEAPopulation[PS].MS) && (Population[PS].TEC > tempCCEAPopulation[PS].TEC))
            {
                //Population.push_back(tempCCEAPopulation[PS]);
                Population[PS].m_SpeedVector=tempCCEAPopulation[PS].m_SpeedVector;
                Population[PS].MS = tempCCEAPopulation[PS].MS;
                Population[PS].TEC = tempCCEAPopulation[PS].TEC;
            }
        }

    }

}

//����
void NOperator::Speed_mutation(vector<Individual>& Population, vector<Individual>& Temp)
{
    vector<Individual> tempCCEAPopulation;
    tempCCEAPopulation.clear();
    int Orgsize = Population.size();

    for (int PS = 0; PS < Orgsize; PS++)
    {
        //���� ����λʱ��ӹ��ܺ�
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                Temp[PS].m_SpeedVector[j][i] = rand() % 3;
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[Temp[PS].m_SpeedVector[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[Temp[PS].m_SpeedVector[j][i]] * m_Speed[Temp[PS].m_SpeedVector[j][i]];
            }
        }


        vector<int> FacSpan(m_Factories, 0);
        vector<float> FacEC(m_Factories, 0);
        int Makespan = 0;
        float TotalEC = 0;

        GetMSandTECForPerandToalFac(Temp[PS].m_FacFamSeqArray, Temp[PS].m_JobSeqInFamArray, FacSpan, FacEC, Makespan, TotalEC);

        Temp[PS].MS = Makespan;
        Temp[PS].TEC = TotalEC;
        if (((Temp[PS].MS < Population[PS].MS) && (Temp[PS].TEC < Population[PS].TEC)) || ((Temp[PS].MS < Population[PS].MS) && (Temp[PS].TEC == Population[PS].TEC)) || ((Temp[PS].MS == Population[PS].MS) && (Temp[PS].TEC < Population[PS].TEC)))
        {
            Population[PS].m_SpeedVector = Temp[PS].m_SpeedVector;
            Population[PS].MS = Temp[PS].MS;
            Population[PS].TEC = Temp[PS].TEC;

        }

        tempCCEAPopulation.push_back(Population[PS]);

        m_SpeedMatrix = tempCCEAPopulation[PS].m_SpeedVector;
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }

        int tempMakespan = 0;
        float tempTotalEC = 0.0;
        vector<int> tempFacSpan(m_Factories, 0);
        vector<float> tempFacEC(m_Factories, 0);

        GetMSandTECForPerandToalFac(Population[PS].m_FacFamSeqArray, Population[PS].m_JobSeqInFamArray, tempFacSpan, tempFacEC, tempMakespan, tempTotalEC);

        //���ܲ���1
        vector<vector<int>> DelayTime;
        DelayTime.clear();
        DelayTime.resize(this->m_Jobs, vector<int>(this->m_Machines, 0));
        vector<vector<int>> JFDTime(this->m_Jobs);
        for (int j = 0; j < JFDTime.size(); j++)
            JFDTime[j].resize(this->m_Machines, 0);

        GetDelayTime_Forward(tempCCEAPopulation[PS].m_FacFamSeqArray, tempCCEAPopulation[PS].m_JobSeqInFamArray, JFDTime, FacSpan, DelayTime);

        bool sign = false;

        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 1; i < m_Machines; i++)
            {
                if ((DelayTime[j][i] > 0) && (i < m_Machines - 1))
                {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[tempCCEAPopulation[PS].m_SpeedVector[j][i]]);

                    m_TureJobOpertime[j][i + 1] = static_cast<int>(m_JobOperPTime[j][i + 1] / m_Speed[tempCCEAPopulation[PS].m_SpeedVector[j][i + 1]]);

                    if ((JFDTime[j][i] < (JFDTime[j][i + 1] - m_TureJobOpertime[j][i + 1])))
                    {
                        int Speedlevel = tempCCEAPopulation[PS].m_SpeedVector[j][i] - 1;

                        for (int level = Speedlevel; level > 0; level--)
                        {
                            if ((static_cast<int>(m_JobOperPTime[j][i] / m_Speed[level]) - m_TureJobOpertime[j][i]) < DelayTime[j][i])
                            {
                                if ((JFDTime[j][i] + (static_cast<int>(m_JobOperPTime[j][i] / m_Speed[level])) - m_TureJobOpertime[j][i]) < (JFDTime[j][i + 1] - m_TureJobOpertime[j][i + 1]))
                                {
                                    tempCCEAPopulation[PS].m_SpeedVector[j][i] = level;
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
                    int Speedlevel = tempCCEAPopulation[PS].m_SpeedVector[j][i] - 1;
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[tempCCEAPopulation[PS].m_SpeedVector[j][i]]);
                    for (int level = Speedlevel; level > 0; level--)
                    {
                        if ((static_cast<int>(m_JobOperPTime[j][i] / m_Speed[level]) - m_TureJobOpertime[j][i]) < DelayTime[j][i])
                        {
                            tempCCEAPopulation[PS].m_SpeedVector[j][i] = level;
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

            m_SpeedMatrix = tempCCEAPopulation[PS].m_SpeedVector;

            //�õ������Ĵ���ʱ�� ����λʱ��ӹ��ܺ�
            for (int j = 0; j < m_Jobs; j++)
            {
                for (int i = 0; i < m_Machines; i++)
                {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }


            int Makespan1 = 0;
            float TotalEC1 = 0.0;
            vector<int> FacSpan1(m_Factories, 0);
            vector<float> FacEC1(m_Factories, 0);

            GetMSandTECForPerandToalFac(tempCCEAPopulation[PS].m_FacFamSeqArray, tempCCEAPopulation[PS].m_JobSeqInFamArray, FacSpan1, FacEC1, Makespan1, TotalEC1);

            tempCCEAPopulation[PS].MS = Makespan1;
            tempCCEAPopulation[PS].TEC = TotalEC1;

            //���
            //CheckSol(tempCCEAPopulation[PS].m_FacFamSeqArray, tempCCEAPopulation[PS].m_JobSeqInFamArray, tempCCEAPopulation[PS].MS);
            //CheckSolTEC(tempCCEAPopulation[PS].m_FacFamSeqArray, tempCCEAPopulation[PS].m_JobSeqInFamArray, tempCCEAPopulation[PS].TEC);

            //�ж�
            if ((Population[PS].MS == tempCCEAPopulation[PS].MS) && (Population[PS].TEC > tempCCEAPopulation[PS].TEC))
            {
                Population[PS].m_SpeedVector=tempCCEAPopulation[PS].m_SpeedVector;
                Population[PS].MS = tempCCEAPopulation[PS].MS;
                Population[PS].TEC = tempCCEAPopulation[PS].TEC;
            }
        }

    }

}


void NOperator::Speed_mutation(vector<Individual>& CMOEAPopulation, vector<Individual>& Temp, vector<Individual>& tureCMOEAPopulation)
{

    for (int PS = 0; PS < CMOEAPopulation.size(); PS++)
    {
        //���� ����λ�ӹ�ʱ���ܺ�
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                Temp[PS].m_SpeedVector[j][i] = rand() % 3;
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[Temp[PS].m_SpeedVector[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[Temp[PS].m_SpeedVector[j][i]] * m_Speed[Temp[PS].m_SpeedVector[j][i]];
            }
        }


        vector<int> FacSpan(m_Factories, 0);
        vector<float> FacEC(m_Factories, 0);
        int Makespan = 0;
        float TotalEC = 0;
        Makespan = GetJFDTime_Forward(Temp[PS].m_FacFamSeqArray, Temp[PS].m_JobSeqInFamArray, m_TempJFDTime, FacSpan);
        TotalEC = GetTECForAllFacByJFD(Temp[PS].m_FacFamSeqArray, Temp[PS].m_JobSeqInFamArray, m_TempJFDTime);
        Temp[PS].MS = Makespan;
        Temp[PS].TEC = TotalEC;
        if (((Temp[PS].MS < CMOEAPopulation[PS].MS) && (Temp[PS].TEC < CMOEAPopulation[PS].TEC)) || ((Temp[PS].MS < CMOEAPopulation[PS].MS) && (Temp[PS].TEC == CMOEAPopulation[PS].TEC)) || ((Temp[PS].MS == CMOEAPopulation[PS].MS) && (Temp[PS].TEC < CMOEAPopulation[PS].TEC)))
        {
            CMOEAPopulation[PS].m_SpeedVector = Temp[PS].m_SpeedVector;
            CMOEAPopulation[PS].MS = Temp[PS].MS;
            CMOEAPopulation[PS].TEC = Temp[PS].TEC;
            tureCMOEAPopulation.push_back(Temp[PS]);  //����ο���
        }

    }

}


void NOperator::Speed_mutation_NoAC(vector<Individual>& QCCEAPopulation, vector<Individual>& Temp, vector<Individual>& tureQCCEAPopulation)
{
    //����NR_s�еĸ�������ٶ�ͻ�䣬(��,��,v)to(��,��,v��)?�����½�ռ��ԭ�⣬����v'����v?

    for (int PS = 0; PS < QCCEAPopulation.size(); PS++)
    {
        //���� ����λ�ӹ�ʱ���ܺ�
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int m = 0; m < m_Machines; m++)
            {
                Temp[PS].m_SpeedVector[j][m] = rand() % 3;
            }
        }

        vector<int> FacSpan(m_Factories, 0);
        vector<float> FacEC(m_Factories, 0);
        int Makespan = 0;
        float TotalEC = 0;
        GetMSandTECForPerandToalFac(Temp[PS].m_FacFamSeqArray, Temp[PS].m_JobSeqInFamArray,FacSpan, FacEC, Makespan, TotalEC);
        Temp[PS].MS = Makespan;
        Temp[PS].TEC = TotalEC;
        if (((Temp[PS].MS < QCCEAPopulation[PS].MS) && (Temp[PS].TEC < QCCEAPopulation[PS].TEC)) || ((Temp[PS].MS < QCCEAPopulation[PS].MS) && (Temp[PS].TEC == QCCEAPopulation[PS].TEC)) ||
            ((Temp[PS].MS == QCCEAPopulation[PS].MS) && (Temp[PS].TEC < QCCEAPopulation[PS].TEC)))
        {
            QCCEAPopulation[PS].m_SpeedVector = Temp[PS].m_SpeedVector;
            QCCEAPopulation[PS].MS = Temp[PS].MS;
            QCCEAPopulation[PS].TEC = Temp[PS].TEC;
            tureQCCEAPopulation.push_back(Temp[PS]);  //����ο���
        }
    }

}

void NOperator::SwapFam(vector <int>& FamSeq)
{
    this->SwapFaminFac(FamSeq);
}


//��������
void NOperator::SwapFam(vector<vector <int>>& FamSeq)
{
    if (rand() % 2)
        this->SwapFamBetweenFacs(FamSeq);
    else
        this->SwapFaminFac(FamSeq[rand() % this->m_Factories]);
}

void NOperator::SwapFaminFac(vector <int>& FamSeqInFac)
{
    if (FamSeqInFac.size() < 2) return;
    int pt1 = rand() % FamSeqInFac.size();
    int pt2;
    do {
        pt2 = rand() % FamSeqInFac.size();
    } while (pt1 == pt2);

    int Fam = FamSeqInFac[pt1];
    FamSeqInFac[pt1] = FamSeqInFac[pt2];
    FamSeqInFac[pt2] = Fam;
}
void NOperator::SwapFamBetweenFacs(vector<vector <int>>& FamSeq)
{
    int F1 = rand() % this->m_Factories;  //���ѡ��һ������
    int F2;
    do
    {
        F2 = rand() % this->m_Factories;
    } while (F1 == F2);  //ʹ������������ͬ

    if (FamSeq[F1].size() && FamSeq[F2].size())  //��������������������Ϊ��
    {
        int pt1 = rand() % FamSeq[F1].size();  //���ѡһ��λ��
        int pt2 = rand() % FamSeq[F2].size();
        //����������
        int Fam = FamSeq[F1][pt1];
        FamSeq[F1][pt1] = FamSeq[F2][pt2];
        FamSeq[F2][pt2] = Fam;
    }
    else if (FamSeq[F1].size())  //����F1��������Ϊ0
    {
        int pt1 = rand() % FamSeq[F1].size();  //�ڹ���F1�����ѡһ��λ��
        //����ѡ������빤��F2������F1ɾ����
        int Fam = FamSeq[F1][pt1];
        FamSeq[F2].push_back(Fam);
        FamSeq[F1].erase(FamSeq[F1].begin() + pt1);
    }
    else if (FamSeq[F2].size())
    {
        int pt1 = rand() % FamSeq[F2].size();
        int Fam = FamSeq[F2][pt1];
        FamSeq[F1].push_back(Fam);
        FamSeq[F2].erase(FamSeq[F2].begin() + pt1);
    }
}


void NOperator::JobInsert(vector <int>FamSeq, vector <vector <int>>& JobSeqinFam, int& Span) {
    for (int index = 0; index < FamSeq.size(); index++) {
        int d = FamSeq[index];
        vector <int>JobSeq = JobSeqinFam[d];
        int seqlen = JobSeq.size();
        if (seqlen == 1) {
            continue;
        }
        vector <int>JobPermutation = JobSeq;
        //��������ȡ��
        for (int i = 0; i < seqlen; i++) {
            int insertjob = JobPermutation[i];
            vector<int>::iterator it = find(JobSeq.begin(), JobSeq.end(), insertjob);
            JobSeq.erase(it);
            int positions = JobSeq.size() + 1;
            vector<int>fits;
            fits.clear();
            for (int pos = 0; pos < positions; pos++) {
                vector <int >TempJobSeq = JobSeq;
                TempJobSeq.insert(TempJobSeq.begin() + pos, insertjob);
                JobSeqinFam[d] = TempJobSeq;
                vector<vector<int> > JFDTime(m_Jobs, vector<int>(m_Machines, 0));
                int Camx = GetJFDTime_Forward_InFactory(FamSeq, JobSeqinFam, JFDTime);
                fits.push_back(Camx);
            }
            int bestpos = min_element(fits.begin(), fits.end()) - fits.begin();
            JobSeq.insert(JobSeq.begin() + bestpos, insertjob);
        }
        JobSeqinFam[d] = JobSeq;
    }
    Span = GetSpan(FamSeq, JobSeqinFam);
}

void NOperator::JobSwap(vector <int>FamSeq, vector <vector <int>>& JobSeqinFam, int& Span) {
    for (int d = 0; d < FamSeq.size(); d++) {
        vector <int >JobSeq = JobSeqinFam[d];
        int len = JobSeq.size();
        if (len == 1) {
            continue;
        }
        int x = -1, y = -1;
        vector <vector <int>>TempJobSeqinFam = JobSeqinFam;
        for (int i = 0; i < len - 1; i++) {
            for (int j = i + 1; j < len; j++) {
                vector <int >Tempjobseq = JobSeq;

                int job = Tempjobseq[i];
                Tempjobseq[i] = Tempjobseq[j];
                Tempjobseq[j] = job;

                TempJobSeqinFam[d] = Tempjobseq;
                int TempSpan = GetSpan(FamSeq, TempJobSeqinFam);
                if (TempSpan < Span) {
                    Span = TempSpan;
                    x = i;
                    y = j;
                }

            }
        }
        if (x != -1 && y != -1) {
            int job = JobSeq[x];
            JobSeq[x] = JobSeq[y];
            JobSeq[y] = job;
            JobSeqinFam[d] = JobSeq;
        }
    }
}


void NOperator::SwapJob(vector<vector<int>>& JobSeqinFam)
{
    int r = rand() % this->m_Families;
    if (JobSeqinFam[r].size() > 2)
    {
        int pt1 = rand() % JobSeqinFam[r].size();
        int pt2;
        do {
            pt2 = rand() % JobSeqinFam[r].size();
        } while (pt1 == pt2);

        int Job = JobSeqinFam[r][pt1];
        JobSeqinFam[r][pt1] = JobSeqinFam[r][pt2];
        JobSeqinFam[r][pt2] = Job;
    }
    else if (JobSeqinFam[r].size() == 2)
    {
        int Job = JobSeqinFam[r][0];
        JobSeqinFam[r][0] = JobSeqinFam[r][1];
        JobSeqinFam[r][1] = Job;
    }
}
void NOperator::FamInsert(vector <int>& FamSeq, vector <vector <int>> JobSeqinFam, int& Span)
{
    vector<vector<int>> JBDTime(this->m_Jobs), JFDTime;
    for (int j = 0; j < JBDTime.size(); j++)
        JBDTime[j].resize(this->m_Machines);
    JFDTime = JBDTime;
    int  bestPos;
    for (int CurFam = 0; CurFam < FamSeq.size(); CurFam++)
    {
        // Extract a Fam from FamSeq;

        int InsFam = FamSeq[CurFam];
        vector<int>::iterator it = find(FamSeq.begin(), FamSeq.end(), InsFam);
        int OrgPos = it - FamSeq.begin();
        FamSeq.erase(it);

        //Insert the Fam to the best position
        this->GetJFDTime_Forward_InFactory(FamSeq, JobSeqinFam,  JFDTime);
        this->GetJBDTime_Backward_InFactory(FamSeq, JobSeqinFam, JBDTime);
        int TempSpan = this->FindBestPosToInsertFam_InFactory(FamSeq, JobSeqinFam, JFDTime,JBDTime, InsFam, bestPos);
        if (TempSpan < Span)
        {
            Span = TempSpan;
            FamSeq.insert(FamSeq.begin() + bestPos, InsFam);
        }
        else
        {
            FamSeq.insert(FamSeq.begin() + OrgPos, InsFam);
        }

    }
}


//xin
void NOperator::IG_DR(vector <vector<int>>FamSeqinFac, vector <vector <int>>& JobSeqinFam, vector<int> FacSpan)
{
    int cirFac = max_element(FacSpan.begin(), FacSpan.end()) - FacSpan.begin();
    //double r = double(rand()) / RAND_MAX;
    int Pos = -1;
    int FamExt = -1;
    do {
        Pos = rand() % FamSeqinFac[cirFac].size();
        FamExt = FamSeqinFac[cirFac][Pos];
    } while (JobSeqinFam[FamExt].size() < 2);

    vector<int> JobsExtracted;
    JobsExtracted.clear();

    for (int j = 0; j < JobSeqinFam[FamExt].size() / 2; j++)
    {
        int JobPos = rand() % JobSeqinFam[FamExt].size();
        int Job = JobSeqinFam[FamExt][JobPos];	//��ѡ����
        JobSeqinFam[FamExt].erase(JobSeqinFam[FamExt].begin() + JobPos);	//ɾ������
        JobsExtracted.push_back(Job);	//����Extracted
    }

    vector<int> JobSeq = JobSeqinFam[FamExt];

    for (int j = 0; j < JobsExtracted.size(); j++)
    {
        vector<int>fits;
        fits.clear();
        int positions = JobSeq.size() + 1;
        for (int pos = 0; pos < positions; pos++) {
            vector <int >TempJobSeq = JobSeq;
            TempJobSeq.insert(TempJobSeq.begin() + pos, JobsExtracted[j]);

            vector<vector<int> > JFDTime(m_Jobs, vector<int>(m_Machines, 0));
            int Camx = GetJFDTime_Forward_InFactory(FamSeqinFac[cirFac], JobSeqinFam, JFDTime);
            fits.push_back(Camx);
        }
        int bestpos = min_element(fits.begin(), fits.end()) - fits.begin();
        JobSeq.insert(JobSeq.begin() + bestpos, JobsExtracted[j]);
    }

    JobSeqinFam[FamExt] = JobSeq;

    //int Span = GetSpan(FamSeqinFac, JobSeqinFam);
}


//xin
void NOperator::InsertFamBetweenFac(vector <vector <int>>& FamSeqinFac, vector <vector <int>> JobSeqinFam)
{
    int Fac = -1;
    do
    {
        Fac = rand() % this->m_Factories;
    } while (FamSeqinFac[Fac].size() < 2);

    int pt = rand() % FamSeqinFac[Fac].size();

    int InsFam = FamSeqinFac[Fac][pt];

    FamSeqinFac[Fac].erase(FamSeqinFac[Fac].begin() + pt);

    vector<vector<int>> JBDTime(this->m_Jobs), JFDTime;
    for (int j = 0; j < JBDTime.size(); j++)
        JBDTime[j].resize(this->m_Machines);
    JFDTime = JBDTime;

    for (int fac = 0; fac < FamSeqinFac.size(); fac++)
    {
        this->GetJFDTime_Forward_InFactory(FamSeqinFac[fac], JobSeqinFam, JFDTime);  //ǰ�����
        this->GetJBDTime_Backward_InFactory(FamSeqinFac[Fac], JobSeqinFam, JBDTime); //�������
    }

    int bestFac = -1, bestPos = -1;  //��õĹ�����λ��
    int Span = this-> FindBestPosToInsertFam(FamSeqinFac, JobSeqinFam,JFDTime,JBDTime,  InsFam, bestFac, bestPos);  //�����й������ҵ���õ�λ�ò�����õ��깤ʱ��

    // ����õĹ��������λ���в��뵱ǰ�� Insert CurFam to bestPos at bestFac
    FamSeqinFac[bestFac].insert(FamSeqinFac[bestFac].begin() + bestPos, InsFam);
}



//�� ����ָ����������й����������λ��
float NOperator::FindBestPosToInsertFamForAllFac_Ind(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& JCTime, const vector<vector<int>>& JSTime,
                                                     int InsFam, int& BestFac, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, int OrgMS, int OrgTEC, vector<Individual>& CMOEAPopulation)
{
    float minInd = INT_MAX;
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        int Pos = -1;
        float Ind = this->FindBestPosToInsertFamForPerFac_Ind(Fac, FacFamSeq, FacFamSeq[Fac], JobSeqInFam, JCTime, JSTime,
                                                              InsFam, Pos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, OrgMS, OrgTEC, CMOEAPopulation);
        if (Ind < minInd)
        {
            minInd = Ind;
            BestFac = Fac;
            BestPos = Pos;
        }
    }
    return minInd;
}


//�� ����ָ������ڵ�ǰ���������λ��
float NOperator::FindBestPosToInsertFamForPerFac_Ind(int Fac, const vector<vector<int>>& FacFamSeq, const vector<int>& NewFamSeqInFac, const vector<vector<int>>& JobSeqInFam,
                                                     const vector<vector<int>>& JCTime, const vector<vector<int>>& JSTime, int InsFam, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, int OrgMS, int OrgTEC, vector<Individual>& CMOEAPopulation)
{
    float minFacInd = INT_MAX;

    for (int Pos = 0; Pos <= NewFamSeqInFac.size(); Pos++)
    {

        float FacInd = GetIndForPerFacAfterInsertFam(Fac, FacFamSeq, NewFamSeqInFac, JobSeqInFam, JCTime, JSTime,
                                                     InsFam, Pos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, OrgMS, OrgTEC, CMOEAPopulation);

        if (FacInd < minFacInd)
        {
            minFacInd = FacInd;
            BestPos = Pos;
        }
    }
    return minFacInd;
}

float NOperator::FindBestPosToInsertFamForAllFac_Ind_NoAC(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam,
                                                          int InsFam, int& BestFac, int& BestPos, int nadirpointMS, float nadirpointTEC,
                                                          int idealpointMS, float idealpointTEC, int OrgMS, int OrgTEC, vector<Individual>& CMOEAPopulation)
{
    float minInd = INT_MAX;
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        int Pos = -1;
        float Ind = this->FindBestPosToInsertFamForPerFac_Ind_NoAC(Fac, FacFamSeq, FacFamSeq[Fac], JobSeqInFam,InsFam, Pos,
                                                                   nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, OrgMS, OrgTEC, CMOEAPopulation);
        if (Ind < minInd)
        {
            minInd = Ind;
            BestFac = Fac;
            BestPos = Pos;
        }
    }
    return minInd;
}

//�� ����ָ������ڵ�ǰ���������λ��
float NOperator::FindBestPosToInsertFamForPerFac_Ind_NoAC(int Fac, const vector<vector<int>>& FacFamSeq, const vector<int>& NewFamSeqInFac, const vector<vector<int>>& JobSeqInFam,
                                                          int InsFam, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, int OrgMS,
                                                          int OrgTEC, vector<Individual>& CMOEAPopulation)
{
    float minFacInd = INT_MAX;

    for (int Pos = 0; Pos <= NewFamSeqInFac.size(); Pos++)
    {

        float FacInd = GetIndForPerFacAfterInsertFam_NoAC(Fac, FacFamSeq, NewFamSeqInFac, JobSeqInFam,
                                                          InsFam, Pos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, OrgMS, OrgTEC, CMOEAPopulation);

        if (FacInd < minFacInd)
        {
            minFacInd = FacInd;
            BestPos = Pos;
        }
    }
    return minFacInd;
}

//�� ����ָ����������й����������λ��
float NOperator::FindBestPosToInsertFamForAllFac_Ind_DR(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& JFDTime, const vector<vector<int>>& JBDTime,
                                                        int InsFam, int& BestFac, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC)
{
    float minInd = INT_MAX;
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        int Pos = -1;
        float Ind = this->FindBestPosToInsertFamForPerFac_Ind_DR(Fac, FacFamSeq, FacFamSeq[Fac], JobSeqInFam, JFDTime, JBDTime,
                                                                 InsFam, Pos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);
        if (Ind < minInd)
        {
            minInd = Ind;
            BestFac = Fac;
            BestPos = Pos;
        }
    }
    return minInd;
}

//�� ����ָ������ڵ�ǰ���������λ��DR
float NOperator::FindBestPosToInsertFamForPerFac_Ind_DR(int Fac, const vector<vector<int>>& FacFamSeq, const vector<int>& NewFamSeqInFac, const vector<vector<int>>& JobSeqInFam,
                                                        const vector<vector<int>>& JFDTime, const vector<vector<int>>& JBDTime, int InsFam, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC)
{
    float minFacInd = INT_MAX;

    for (int Pos = 0; Pos <= NewFamSeqInFac.size(); Pos++)
    {

        float FacInd = GetIndForPerFacAfterInsertFam_DR(Fac, FacFamSeq, NewFamSeqInFac, JobSeqInFam, JFDTime, JBDTime,
                                                        InsFam, Pos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC);

        if (FacInd < minFacInd)
        {
            minFacInd = FacInd;
            BestPos = Pos;
        }
    }
    return minFacInd;
}


//�� ĳ��λ�ò���������ָ��Ind
float NOperator::GetIndForPerFacAfterInsertFam_DR(int Fac, const vector<vector<int>>& FacFamSeq, const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam,
                                                  const vector<vector<int>>& JFDTime, const vector<vector<int>>& JBDTime, int InsFam, int Pos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC)
{
    if (Pos == 0)
    {
        for (int m = 0; m < this->m_Machines; m++)
        {
            MachReadyTime[m] = this->m_SetupTime[m][InsFam][InsFam]; //MachReadyTime ��Ϊ׼��ʱ��
        }
    }
    else
    {
        int PreFam = FamSeqInFac[Pos - 1]; //ǰһ����
        int LastJobinPreFam = JobSeqInFam[PreFam][JobSeqInFam[PreFam].size() - 1]; //ǰһ��������һ������
        for (int m = 0; m < this->m_Machines; m++)
        {
            //ǰһ�������һ���������깤ʱ��+׼��ʱ��
            MachReadyTime[m] = JFDTime[LastJobinPreFam][m] + this->m_SetupTime[m][PreFam][InsFam];
        }
    }

    // ��������﹤������ Scheduling jobs in InsFam
    for (int j = 0; j < JobSeqInFam[InsFam].size(); j++)
    {
        int CurJob = JobSeqInFam[InsFam][j]; //��ǰ����
        MachReadyTime[0] = max(MachReadyTime[0] + this->m_TureJobOpertime[CurJob][0],MachReadyTime[1]);
        for (int m = 1; m < this->m_Machines; m++)
        {
            if (m == this->m_Machines - 1)
            {
                //���һ̨����
                MachReadyTime[m] = MachReadyTime[m - 1] + this->m_TureJobOpertime[CurJob][m];
            } else
            {
                MachReadyTime[m] = max(MachReadyTime[m - 1] + this->m_TureJobOpertime[CurJob][m],MachReadyTime[m + 1]);
            }
        }
    }

    int Span = 0;
    if (Pos < FamSeqInFac.size())
    {
        int NextFam = FamSeqInFac[Pos]; //��һ����

        int FirstJobinNextFam = JobSeqInFam[NextFam][0]; //��һ����ĵ�һ������

        for (int m = 0; m < this->m_Machines; m++)
        {
            if (MachReadyTime[m] + this->m_SetupTime[m][InsFam][NextFam] + JBDTime[FirstJobinNextFam][m] >= Span)
            {
                Span = MachReadyTime[m] + this->m_SetupTime[m][InsFam][NextFam] + JBDTime[FirstJobinNextFam][m];
            }
        }
    }
    else
    {
        Span = MachReadyTime[this->m_Machines - 1];
    }



    //����δ����֮ǰ��������span
    vector <int> FacSpan(FacFamSeq.size(), 0);
    //cout << endl;
    for (int fac = 0; fac < FacFamSeq.size(); fac++)
    {
        int LastFam = *FacFamSeq[fac].rbegin();
        int LastJobInLstFam = *JobSeqInFam[LastFam].rbegin(); //������һ������
        FacSpan[fac] = JFDTime[LastJobInLstFam][this->m_Machines - 1];
        //cout << "δ����ʱ����" << fac << "��Span��" << FacSpan[fac] << endl;
    }

    FacSpan[Fac] = Span;

    //cout << "�ڹ���" << Fac << "��" << Pos << "λ�ò���" << endl;
    //cout << "����" << Fac << "�������Span��" << FacSpan[Fac] << endl;
    int Makespan = *max_element(FacSpan.begin(), FacSpan.end());
    //cout << "����" << Fac << "�������MakeSpan��" << Makespan << endl;

    //�ж��Ƿ��滻��͵�������
    if (Makespan > nadirpointMS)
        nadirpointMS = Makespan;
    if (Makespan < idealpointMS)
        idealpointMS = Makespan;

    vector<float> FacEC(FacFamSeq.size(), 0);
    float TEC;
    //�����ܺ�
    vector<vector<int>> TempFacFamSeq;
    TempFacFamSeq.clear();
    TempFacFamSeq.resize(FacFamSeq.size());
    TempFacFamSeq = FacFamSeq;
    TempFacFamSeq[Fac].insert(TempFacFamSeq[Fac].begin() + Pos, InsFam);
    for (int fac = 0; fac < TempFacFamSeq.size(); fac++)
    {
        //�����ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
            //int CurFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
            int Fam = TempFacFamSeq[fac][0];
            for (int g = 0; g < TempFacFamSeq[fac].size(); g++)
            {
                Fam = TempFacFamSeq[fac][g];

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

    TEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);
    //cout << "���ܺģ�" << TEC << endl;

    //�ж��Ƿ��滻��͵�������
    if (TEC > nadirpointTEC)
        nadirpointTEC = TEC;
    if (TEC < idealpointTEC)
        idealpointTEC = TEC;


    //��һ��
    float normalMS = -1;
    float normalTEC = -1;
    //cout << endl << "normalize" << endl;

    normalMS = (static_cast<float>(Makespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    normalTEC = ((TEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
    //cout << "�����һ�����MS��" << normalMS << "\tTEC��" << normalTEC << endl;

    //����ָ��Ind
    float convergence_ind = 0;
    convergence_ind += (normalMS - 1.0) * (normalMS - 1.0);
    convergence_ind += (normalTEC - 1.0) * (normalTEC - 1.0);
    convergence_ind = sqrt(convergence_ind);
    convergence_ind = 1 / (convergence_ind + 1);
    //cout << "���������ָ�꣺" << convergence_ind << endl;
    //cout << endl;

    return convergence_ind;
}

//�� ĳ��λ�ò��빤�������ָ��Ind
float NOperator::GetIndForPerFacAfterInsertJob_DR(int fac, const vector<vector<int>>& FacFamSeq, const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam,
                                                  const vector<vector<int>>& JFDTime, const vector<vector<int>>& JBDTime, int InsFam, int FamPos, int InsJob, int JobPos,
                                                  int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC)
{

    // ����λ��֮ǰ����
    for (int f = 0; f < FamPos; ++f)
    {
        int Fam = FamSeqInFac[f];
        int LastJobInFam = JobSeqInFam[Fam][JobSeqInFam[Fam].size() - 1]; //������һ������

    }

    if (FamPos == 0)
    {
        if (JobPos == 0)
        {
            for (int m = 0; m < this->m_Machines; m++)
            {
                MachReadyTime[m] = this->m_SetupTime[m][InsFam][InsFam]; //MachReadyTime ��Ϊ׼��ʱ��
            }
        }
        else
        {
            int PreJob = JobSeqInFam[InsFam][JobPos - 1];
            for (int m = 0; m < this->m_Machines; m++)
            {
                MachReadyTime[m] = JFDTime[PreJob][m];
            }
        }
    }
    else
    {
        if (JobPos == 0)
        {
            int PreFam = FamSeqInFac[FamPos - 1]; //ǰһ����
            int LastJobinPreFam = JobSeqInFam[PreFam][JobSeqInFam[PreFam].size() - 1]; //ǰһ��������һ������
            for (int m = 0; m < this->m_Machines; m++)
            {
                //ǰһ�������һ���������깤ʱ��+׼��ʱ��
                MachReadyTime[m] = JFDTime[LastJobinPreFam][m] + this->m_SetupTime[m][PreFam][InsFam];
            }
        }
        else
        {
            int PreJob = JobSeqInFam[InsFam][JobPos - 1];
            for (int m = 0; m < this->m_Machines; m++)
            {
                MachReadyTime[m] = JFDTime[PreJob][m];
            }
        }
    }

    //����Ĺ���
    MachReadyTime[0] = max(MachReadyTime[0] + this->m_TureJobOpertime[InsJob][0],MachReadyTime[1]);
    for (int m = 1; m < this->m_Machines; m++)
    {
        if (m == this->m_Machines - 1)
        {
            //���һ̨����
            MachReadyTime[m] = MachReadyTime[m - 1] + this->m_TureJobOpertime[InsJob][m];
        } else
        {
            MachReadyTime[m] = max(MachReadyTime[m - 1] + this->m_TureJobOpertime[InsJob][m],MachReadyTime[m + 1]);
        }
    }


    int Span = 0;
    if (FamPos < FamSeqInFac.size() - 1)
    {
        if (JobPos == JobSeqInFam[InsFam].size())
        {
            int NextFam = FamSeqInFac[FamPos + 1]; //��һ����

            int FirstJobInNextFam = JobSeqInFam[NextFam][0]; //��һ����ĵ�һ������

            for (int m = 0; m < this->m_Machines; m++)
            {
                if (MachReadyTime[m] + this->m_SetupTime[m][InsFam][NextFam] + JBDTime[FirstJobInNextFam][m] >= Span)
                {
                    Span = MachReadyTime[m] + this->m_SetupTime[m][InsFam][NextFam] + JBDTime[FirstJobInNextFam][m];

                }
            }
        }
        else
        {
            int NextJob = JobSeqInFam[InsFam][JobPos]; //��һ������
            for (int m = 0; m < this->m_Machines; m++)
            {
                if (MachReadyTime[m] + JBDTime[NextJob][m] >= Span)
                {
                    Span = MachReadyTime[m] + JBDTime[NextJob][m];

                }
            }
        }
    }
    else
    {
        if (JobPos == JobSeqInFam[InsFam].size())
        {
            Span = MachReadyTime[this->m_Machines - 1];

        }
        else
        {
            int NextJob = JobSeqInFam[InsFam][JobPos]; //��һ����
            for (int m = 0; m < this->m_Machines; m++)
            {
                if (MachReadyTime[m] + JBDTime[NextJob][m] >= Span)
                {
                    Span = MachReadyTime[m] + JBDTime[NextJob][m];
                }
            }
        }
    }

    //����δ����֮ǰ��������span
    vector <int> FacSpan(FacFamSeq.size(), 0);
    //cout << endl;
    for (int f = 0; f < FacFamSeq.size(); f++)
    {
        int LastFam = *FacFamSeq[f].rbegin();
        int LastJobInLstFam = *JobSeqInFam[LastFam].rbegin(); //������һ������
        FacSpan[f] = JFDTime[LastJobInLstFam][this->m_Machines - 1];
        //cout << "δ����ʱ����" << f << "��Span��" << FacSpan[f] << endl;
    }

    FacSpan[fac] = Span;

    //cout << "���빤��" << fac << "���Span��" << FacSpan[fac] << endl;
    int Makespan = *max_element(FacSpan.begin(), FacSpan.end());
    //cout << "�����MakeSpan��" << Makespan << endl;

    //�ж��Ƿ��滻��͵�������
    if (Makespan > nadirpointMS)
        nadirpointMS = Makespan;
    if (Makespan < idealpointMS)
        idealpointMS = Makespan;

    vector<float> FacEC(FacFamSeq.size(), 0);
    float TEC;
    //�����ܺ�
    vector<vector<int>> TempFacFamSeq;
    TempFacFamSeq.clear();
    TempFacFamSeq.resize(FacFamSeq.size());
    TempFacFamSeq = FacFamSeq;

    vector<vector<int>> TempJobSeqInFam;
    TempJobSeqInFam.clear();
    TempJobSeqInFam.resize(JobSeqInFam.size());
    TempJobSeqInFam = JobSeqInFam;
    TempJobSeqInFam[InsFam].insert(TempJobSeqInFam[InsFam].begin() + JobPos, InsJob);

    for (int fac = 0; fac < TempFacFamSeq.size(); fac++)
    {
        //�����ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
            //int CurFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
            int Fam = TempFacFamSeq[fac][0];
            for (int g = 0; g < TempFacFamSeq[fac].size(); g++)
            {
                Fam = TempFacFamSeq[fac][g];

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
                for (int j = 0; j < TempJobSeqInFam[Fam].size(); j++)
                {
                    int Job = TempJobSeqInFam[Fam][j];
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

    TEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);
    //cout << "���ܺģ�" << TEC << endl;

    //�ж��Ƿ��滻��͵�������
    if (TEC > nadirpointTEC)
        nadirpointTEC = TEC;
    if (TEC < idealpointTEC)
        idealpointTEC = TEC;


    //��һ��
    float normalMS = -1;
    float normalTEC = -1;
    //cout << endl << "normalize" << endl;

    normalMS = (static_cast<float>(Makespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    normalTEC = ((TEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
    //cout << "�����һ�����MS��" << normalMS << "\tTEC��" << normalTEC << endl;

    //����ָ��Ind
    float convergence_ind = 0;
    convergence_ind += (normalMS - 1.0) * (normalMS - 1.0);
    convergence_ind += (normalTEC - 1.0) * (normalTEC - 1.0);
    convergence_ind = sqrt(convergence_ind);
    convergence_ind = 1 / (convergence_ind + 1);
    //cout << "���������ָ�꣺" << convergence_ind << endl;
    //cout << endl;

    return convergence_ind;
}

//�� ĳ��λ�ò���������ָ��Ind
float NOperator::GetIndForPerFacAfterInsertFam(int Fac, const vector<vector<int>>& FacFamSeq, const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam,
                                               const vector<vector<int>>& JCTime, const vector<vector<int>>& JSTime, int InsFam, int Pos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, int OrgMS, int OrgTEC, vector<Individual>& CMOEAPopulation)
{
    vector<int> ForwardMachReadyTime(this->m_Machines, 0);

    if (Pos == 0)
    {
        for (int m = 0; m < this->m_Machines; m++)
        {
            ForwardMachReadyTime[m] = this->m_SetupTime[m][InsFam][InsFam]; //MachReadyTime ��Ϊ׼��ʱ��
        }
    }
    else
    {
        int PreFam = FamSeqInFac[Pos - 1]; //ǰһ����
        int LastJobinPreFam = JobSeqInFam[PreFam][JobSeqInFam[PreFam].size() - 1]; //ǰһ��������һ������
        for (int m = 0; m < this->m_Machines; m++)
        {
            //ǰһ�������һ���������깤ʱ��+׼��ʱ��
            ForwardMachReadyTime[m] = JCTime[LastJobinPreFam][m] + this->m_SetupTime[m][PreFam][InsFam];
        }
    }

    // ��������﹤������ Scheduling jobs in InsFam
    for (int j = 0; j < JobSeqInFam[InsFam].size(); j++)
    {
        int CurJob = JobSeqInFam[InsFam][j]; //��ǰ����
        ForwardMachReadyTime[0] += this->m_TureJobOpertime[CurJob][0]; //�ڵ�һ��������
        for (int m = 1; m < this->m_Machines; m++)
        {
            ForwardMachReadyTime[m] = max(ForwardMachReadyTime[m - 1], ForwardMachReadyTime[m]) + this->m_TureJobOpertime[CurJob][m];
        }
    }

    int Span = 0;
    if (Pos < FamSeqInFac.size())
    {
        int NextFam = FamSeqInFac[Pos]; //��һ����

        int FirstJobinNextFam = JobSeqInFam[NextFam][0]; //��һ����ĵ�һ������

        for (int m = 0; m < this->m_Machines; m++)
        {
            if (ForwardMachReadyTime[m] + this->m_SetupTime[m][InsFam][NextFam] + JSTime[FirstJobinNextFam][m] >= Span)
            {
                Span = ForwardMachReadyTime[m] + this->m_SetupTime[m][InsFam][NextFam] + JSTime[FirstJobinNextFam][m];
            }
        }
    }
    else
    {
        Span = ForwardMachReadyTime[this->m_Machines - 1];
    }



    //����δ����֮ǰ��������span
    vector <int> FacSpan(FacFamSeq.size(), 0);
    //cout << endl;
    for (int fac = 0; fac < FacFamSeq.size(); fac++)
    {
        int LastFam = *FacFamSeq[fac].rbegin();
        int LastJobInLstFam = *JobSeqInFam[LastFam].rbegin(); //������һ������
        FacSpan[fac] = JCTime[LastJobInLstFam][this->m_Machines - 1];
        //cout << "δ����ʱ����" << fac << "��Span��" << FacSpan[fac] << endl;
    }

    FacSpan[Fac] = Span;

    //cout << "�ڹ���" << Fac << "��" << Pos << "λ�ò���" << endl;
    //cout << "����" << Fac << "�������Span��" << FacSpan[Fac] << endl;
    int Makespan = *max_element(FacSpan.begin(), FacSpan.end());
    //cout << "����" << Fac << "�������MakeSpan��" << Makespan << endl;

    //�ж��Ƿ��滻��͵�������
    if (Makespan > nadirpointMS)
        nadirpointMS = Makespan;
    if (Makespan < idealpointMS)
        idealpointMS = Makespan;

    vector<float> FacEC(FacFamSeq.size(), 0);
    float TEC;
    //�����ܺ�
    vector<vector<int>> TempFacFamSeq;
    TempFacFamSeq.clear();
    TempFacFamSeq.resize(FacFamSeq.size());
    TempFacFamSeq = FacFamSeq;
    TempFacFamSeq[Fac].insert(TempFacFamSeq[Fac].begin() + Pos, InsFam);
    for (int fac = 0; fac < TempFacFamSeq.size(); fac++)
    {
        //�����ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
            //int CurFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
            int Fam = TempFacFamSeq[fac][0];
            for (int g = 0; g < TempFacFamSeq[fac].size(); g++)
            {
                Fam = TempFacFamSeq[fac][g];

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

    TEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);
    //cout << "���ܺģ�" << TEC << endl;

    //�ж��Ƿ��滻��͵�������
    if (TEC > nadirpointTEC)
        nadirpointTEC = TEC;
    if (TEC < idealpointTEC)
        idealpointTEC = TEC;

    //�ж��Ƿ�Ȳ���ǰ�Ľ������Ľ������ο���
    bool flag = true;
    if (((Makespan < OrgMS) && (TEC < OrgTEC)) || ((Makespan < OrgMS) && (TEC == OrgTEC)) || ((Makespan == OrgMS) && (TEC < OrgTEC)) || (Makespan < idealpointMS) || (TEC < idealpointTEC))
    {
        for (int i = 0; i < CMOEAPopulation.size(); i++)
        {
            if ((Makespan == CMOEAPopulation[i].MS) && (TEC == CMOEAPopulation[i].TEC))
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
            tempIndi.MS = Makespan;
            tempIndi.TEC = TEC;
            tempIndi.m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���
            CMOEAPopulation.push_back(tempIndi);
        }
    }

    //��һ��
    float normalMS = -1;
    float normalTEC = -1;
    //cout << endl << "normalize" << endl;

    normalMS = (static_cast<float>(Makespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    normalTEC = ((TEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
    //cout << "�����һ�����MS��" << normalMS << "\tTEC��" << normalTEC << endl;

    //����ָ��Ind
    float convergence_ind = 0;
    convergence_ind += (normalMS - 1.0) * (normalMS - 1.0);
    convergence_ind += (normalTEC - 1.0) * (normalTEC - 1.0);
    convergence_ind = sqrt(convergence_ind);
    convergence_ind = 1 / (convergence_ind + 1);
    //cout << "���������ָ�꣺" << convergence_ind << endl;
    //cout << endl;

    return convergence_ind;
}


//�� ĳ��λ�ò���������ָ��Ind
float NOperator::GetIndForPerFacAfterInsertFam_NoAC(int Fac, const vector<vector<int>>& FacFamSeq, const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam,
                                                    int InsFam, int Pos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,
                                                    int OrgMS, int OrgTEC, vector<Individual>& CMOEAPopulation)

{
    //����δ����֮ǰ��������span
    vector <int> FacSpan(FacFamSeq.size(), 0);
    vector<vector<int>> TempFacFamSeq;
    TempFacFamSeq.clear();
    TempFacFamSeq.resize(FacFamSeq.size());
    TempFacFamSeq = FacFamSeq;
    TempFacFamSeq[Fac].insert(TempFacFamSeq[Fac].begin() + Pos, InsFam);

    int Makespan = GetSpan_Forward(TempFacFamSeq, JobSeqInFam,FacSpan);

    //�ж��Ƿ��滻��͵�������
    if (Makespan > nadirpointMS)
        nadirpointMS = Makespan;
    if (Makespan < idealpointMS)
        idealpointMS = Makespan;

    vector<float> FacEC(FacFamSeq.size(), 0);
    float TEC;
    for (int fac = 0; fac < TempFacFamSeq.size(); fac++)
    {
        //�����ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
            //int CurFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
            if (TempFacFamSeq[fac].empty()) {
                continue; // ֱ��������ǰ����
            }
            int Fam = TempFacFamSeq[fac][0];
            for (int g = 0; g < TempFacFamSeq[fac].size(); g++)
            {
                Fam = TempFacFamSeq[fac][g];

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

    TEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);
    //cout << "���ܺģ�" << TEC << endl;

    //�ж��Ƿ��滻��͵�������
    if (TEC > nadirpointTEC)
        nadirpointTEC = TEC;
    if (TEC < idealpointTEC)
        idealpointTEC = TEC;

    //�ж��Ƿ�Ȳ���ǰ�Ľ������Ľ������ο���
    bool flag = true;
    if (((Makespan < OrgMS) && (TEC < OrgTEC)) || ((Makespan < OrgMS) && (TEC == OrgTEC)) || ((Makespan == OrgMS) && (TEC < OrgTEC)) || (Makespan < idealpointMS) || (TEC < idealpointTEC))
    {
        for (int i = 0; i < CMOEAPopulation.size(); i++)
        {
            if ((Makespan == CMOEAPopulation[i].MS) && (TEC == CMOEAPopulation[i].TEC))
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
            tempIndi.MS = Makespan;
            tempIndi.TEC = TEC;
            tempIndi.m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���
            CMOEAPopulation.push_back(tempIndi);
        }
    }

    //��һ��
    float normalMS = -1;
    float normalTEC = -1;
    //cout << endl << "normalize" << endl;

    normalMS = (static_cast<float>(Makespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    normalTEC = ((TEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
    //cout << "�����һ�����MS��" << normalMS << "\tTEC��" << normalTEC << endl;

    //����ָ��Ind
    float convergence_ind = 0;
    convergence_ind += (normalMS - 1.0) * (normalMS - 1.0);
    convergence_ind += (normalTEC - 1.0) * (normalTEC - 1.0);
    convergence_ind = sqrt(convergence_ind);
    convergence_ind = 1 / (convergence_ind + 1);
    //cout << "���������ָ�꣺" << convergence_ind << endl;
    //cout << endl;

    return convergence_ind;
}

//�� ĳ��λ�ò��빤�������ָ��Ind
float NOperator::GetIndForPerFacAfterInsertJob(int fac, const vector<vector<int>>& FacFamSeq, const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam, const vector<vector<int>>& JFDTime, const vector<vector<int>>& JBDTime, int InsFam, int FamPos, int InsJob, int JobPos,
                                               int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, int OrgMS, float OrgTEC, vector<Individual>& CMOEAPopulation)
{
    // ����λ��֮ǰ����
    for (int f = 0; f < FamPos; ++f)
    {
        int Fam = FamSeqInFac[f];
        int LastJobInFam = JobSeqInFam[Fam][JobSeqInFam[Fam].size() - 1]; //������һ������

    }

    if (FamPos == 0)
    {
        if (JobPos == 0)
        {
            for (int m = 0; m < this->m_Machines; m++)
            {
                MachReadyTime[m] = this->m_SetupTime[m][InsFam][InsFam]; //MachReadyTime ��Ϊ׼��ʱ��
            }
        }
        else
        {
            int PreJob = JobSeqInFam[InsFam][JobPos - 1];
            for (int m = 0; m < this->m_Machines; m++)
            {
                MachReadyTime[m] = JFDTime[PreJob][m];
            }
        }
    }
    else
    {
        if (JobPos == 0)
        {
            int PreFam = FamSeqInFac[FamPos - 1]; //ǰһ����
            int LastJobinPreFam = JobSeqInFam[PreFam][JobSeqInFam[PreFam].size() - 1]; //ǰһ��������һ������
            for (int m = 0; m < this->m_Machines; m++)
            {
                //ǰһ�������һ���������깤ʱ��+׼��ʱ��
                MachReadyTime[m] = JFDTime[LastJobinPreFam][m] + this->m_SetupTime[m][PreFam][InsFam];
            }
        }
        else
        {
            int PreJob = JobSeqInFam[InsFam][JobPos - 1];
            for (int m = 0; m < this->m_Machines; m++)
            {
                MachReadyTime[m] = JFDTime[PreJob][m];
            }
        }
    }

    //����Ĺ���
    MachReadyTime[0] = max(MachReadyTime[0] + this->m_TureJobOpertime[InsJob][0],MachReadyTime[1]);
    for (int m = 1; m < this->m_Machines; m++)
    {
        if (m == this->m_Machines - 1)
        {
            //���һ̨����
            MachReadyTime[m] = MachReadyTime[m - 1] + this->m_TureJobOpertime[InsJob][m];
        } else
        {
            MachReadyTime[m] = max(MachReadyTime[m - 1] + this->m_TureJobOpertime[InsJob][m], MachReadyTime[m + 1]);
        }
    }


    int Span = 0;
    if (FamPos < FamSeqInFac.size() - 1)
    {
        if (JobPos == JobSeqInFam[InsFam].size())
        {
            int NextFam = FamSeqInFac[FamPos + 1]; //��һ����

            int FirstJobInNextFam = JobSeqInFam[NextFam][0]; //��һ����ĵ�һ������

            for (int m = 0; m < this->m_Machines; m++)
            {
                if (MachReadyTime[m] + this->m_SetupTime[m][InsFam][NextFam] + JBDTime[FirstJobInNextFam][m] >= Span)
                {
                    Span = MachReadyTime[m] + this->m_SetupTime[m][InsFam][NextFam] + JBDTime[FirstJobInNextFam][m];

                }
            }
        }
        else
        {
            int NextJob = JobSeqInFam[InsFam][JobPos]; //��һ������
            for (int m = 0; m < this->m_Machines; m++)
            {
                if (MachReadyTime[m] + JBDTime[NextJob][m] >= Span)
                {
                    Span = MachReadyTime[m] + JBDTime[NextJob][m];

                }
            }
        }
    }
    else
    {
        if (JobPos == JobSeqInFam[InsFam].size())
        {
            Span = MachReadyTime[this->m_Machines - 1];

        }
        else
        {
            int NextJob = JobSeqInFam[InsFam][JobPos]; //��һ����
            for (int m = 0; m < this->m_Machines; m++)
            {
                if (MachReadyTime[m] + JBDTime[NextJob][m] >= Span)
                {
                    Span = MachReadyTime[m] + JBDTime[NextJob][m];
                }
            }
        }
    }

    //����δ����֮ǰ��������span
    vector <int> FacSpan(FacFamSeq.size(), 0);
    //cout << endl;
    for (int f = 0; f < FacFamSeq.size(); f++)
    {
        int LastFam = *FacFamSeq[f].rbegin();
        int LastJobInLstFam = *JobSeqInFam[LastFam].rbegin(); //������һ������
        FacSpan[f] = JFDTime[LastJobInLstFam][this->m_Machines - 1];
        //cout << "δ����ʱ����" << f << "��Span��" << FacSpan[f] << endl;
    }

    FacSpan[fac] = Span;

    //cout << "���빤��" << fac << "���Span��" << FacSpan[fac] << endl;
    int Makespan = *max_element(FacSpan.begin(), FacSpan.end());
    //cout << "�����MakeSpan��" << Makespan << endl;

    //�ж��Ƿ��滻��͵�������
    if (Makespan > nadirpointMS)
        nadirpointMS = Makespan;
    if (Makespan < idealpointMS)
        idealpointMS = Makespan;

    vector<float> FacEC(FacFamSeq.size(), 0);
    float TEC;
    //�����ܺ�
    vector<vector<int>> TempFacFamSeq;
    TempFacFamSeq.clear();
    TempFacFamSeq.resize(FacFamSeq.size());
    TempFacFamSeq = FacFamSeq;

    vector<vector<int>> TempJobSeqInFam;
    TempJobSeqInFam.clear();
    TempJobSeqInFam.resize(JobSeqInFam.size());
    TempJobSeqInFam = JobSeqInFam;
    TempJobSeqInFam[InsFam].insert(TempJobSeqInFam[InsFam].begin() + JobPos, InsJob);

    for (int fac = 0; fac < TempFacFamSeq.size(); fac++)
    {
        //�����ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
            //int CurFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
            int Fam = TempFacFamSeq[fac][0];
            for (int g = 0; g < TempFacFamSeq[fac].size(); g++)
            {
                Fam = TempFacFamSeq[fac][g];

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
                for (int j = 0; j < TempJobSeqInFam[Fam].size(); j++)
                {
                    int Job = TempJobSeqInFam[Fam][j];
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

    TEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);
    //cout << "���ܺģ�" << TEC << endl;

    //�ж��Ƿ��滻��͵�������
    if (TEC > nadirpointTEC)
        nadirpointTEC = TEC;
    if (TEC < idealpointTEC)
        idealpointTEC = TEC;

    //�ж��Ƿ�Ȳ���ǰ�Ľ������Ľ������ο���
    bool flag = true;
    if (((Makespan < OrgMS) && (TEC < OrgTEC)) || ((Makespan < OrgMS) && (TEC == OrgTEC)) || ((Makespan == OrgMS) && (TEC < OrgTEC)) || (Makespan < idealpointMS) || (TEC < idealpointTEC))
    {
        for (int i = 0; i < CMOEAPopulation.size(); i++)
        {
            if ((Makespan == CMOEAPopulation[i].MS) && (TEC == CMOEAPopulation[i].TEC))
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
            tempIndi.m_JobSeqInFamArray = TempJobSeqInFam;
            tempIndi.MS = Makespan;
            tempIndi.TEC = TEC;
            tempIndi.m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���
            CMOEAPopulation.push_back(tempIndi);
        }
    }

    //��һ��
    float normalMS = -1;
    float normalTEC = -1;
    //cout << endl << "normalize" << endl;

    normalMS = (static_cast<float>(Makespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    normalTEC = ((TEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
    //cout << "�����һ�����MS��" << normalMS << "\tTEC��" << normalTEC << endl;

    //����ָ��Ind
    float convergence_ind = 0;
    convergence_ind += (normalMS - 1.0) * (normalMS - 1.0);
    convergence_ind += (normalTEC - 1.0) * (normalTEC - 1.0);
    convergence_ind = sqrt(convergence_ind);
    convergence_ind = 1 / (convergence_ind + 1);
    //cout << "���������ָ�꣺" << convergence_ind << endl;
    //cout << endl;

    return convergence_ind;
}


//�� ����ָ��Ĺ�������
void NOperator::Basedind_JobsInFam(int fac, const vector<vector<int>>& FacFamSeq, const vector<int>& FamSeqInFac, vector<vector<int>>& JobSeqInFam, vector<vector<int>>& JCTime, vector<vector<int>>& JSTime,
                                   int Fam, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual>& CMOEAPopulation, int OrgMS, float OrgTEC, int& ObjectMS, float& ObjectTEC)
{
    vector<int> JobSeqsInFamForExtracted = JobSeqInFam[Fam]; //��ǰ��Ĺ�������
    random_shuffle(JobSeqsInFamForExtracted.begin(), JobSeqsInFamForExtracted.end()); //����
    int FamPos = find(begin(FamSeqInFac), end(FamSeqInFac), Fam) - begin(FamSeqInFac);


    int nCnt = 0, CurPos = 0, BestPos = -1;
    int Makespan = -1;
    //���� ���й����� ��
    while (nCnt < JobSeqInFam[Fam].size())
    {

        CurPos = CurPos % JobSeqsInFamForExtracted.size(); //���ѡһ��λ��
        int CurJob = JobSeqsInFamForExtracted[CurPos]; //ȡ��λ�õĹ���
        vector<int>::iterator it = find(JobSeqInFam[Fam].begin(), JobSeqInFam[Fam].end(), CurJob);
        int OrgJobPos = it - JobSeqInFam[Fam].begin();
        JobSeqInFam[Fam].erase(it);
        RefreshJFDJBDTime_InFactory(FamSeqInFac, JobSeqInFam, JCTime, JSTime, FamPos, FamPos, OrgJobPos, OrgJobPos - 1);

        float minFacInd = INT_MAX;
        //vector<int> TempFamSeqInFac = NewFamSeqInFac;
        for (int Pos = 0; Pos <= JobSeqInFam[Fam].size(); Pos++)
        {
            float FacInd = GetIndForPerFacAfterInsertJob(fac, FacFamSeq, FamSeqInFac, JobSeqInFam, JCTime, JSTime, Fam, FamPos, CurJob, Pos,
                                                         nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, OrgMS, OrgTEC, CMOEAPopulation);

            if (FacInd < minFacInd)
            {
                minFacInd = FacInd;
                BestPos = Pos;
            }

        }
        JobSeqInFam[Fam].insert(JobSeqInFam[Fam].begin() + BestPos, CurJob);
        RefreshJFDJBDTime_InFactory(FamSeqInFac, JobSeqInFam, JCTime, JSTime, FamPos, FamPos, BestPos, BestPos);

        //cout << "�ڹ���" << fac << "��" << Fam << "���" << BestPos << "λ�ò���, ָ��Ϊ��" << minFacInd << endl;

        vector<int> AfterInsertSpanFac(FacFamSeq.size(), 0);
        vector<float> AfterInsertECFac(FacFamSeq.size(), 0);
        int AfterInsertMS = -1;
        float AfterInsertTEC = -1;
        GetMSandTECForPerandToalFacByJFD(FacFamSeq, JobSeqInFam, JCTime, AfterInsertSpanFac, AfterInsertECFac, AfterInsertMS, AfterInsertTEC);
        //cout << "����֮���MS��" << AfterInsertMS << "\tTEC��" << AfterInsertTEC << endl;

        ObjectMS = AfterInsertMS;
        ObjectTEC = AfterInsertTEC;

        CurPos++;
        nCnt++;
        //��麯��
        //CheckSol(FacFamSeq, JobSeqInFam, AfterInsertMS);
        //CheckSolTEC(FacFamSeq, JobSeqInFam, AfterInsertTEC);
    }

}


//Heu
void
NOperator::SortJobsInFam(int SortMethod, vector<vector<int>> &JobSeqInFam) //0:LPT,   1:SPT,  2:JobWeightTotalPTime,
{

    JobSeqInFam = this->m_JobsInEachFamily; //initialize memeory
    if (SortMethod == 0) //LPT
    {
        for (auto &JobSeq: JobSeqInFam) {
            sort(begin(JobSeq), end(JobSeq), [this](int j1, int j2) {
                return this->m_JobTotalPTime[j1] > this->m_JobTotalPTime[j2];
            });
        }
    }

    if (SortMethod == 1) //SPT
    {
        for (auto &JobSeq: JobSeqInFam) {
            sort(begin(JobSeq), end(JobSeq), [this](int j1, int j2) {
                return this->m_JobTotalPTime[j1] < this->m_JobTotalPTime[j2];
            });
        }
    }

    if (SortMethod == 2) //in non-increasing order JobWeightTotalPTime
    {
        for (auto &JobSeq: JobSeqInFam) {
            sort(begin(JobSeq), end(JobSeq), [this](int j1, int j2) {
                return this->m_JobWeightTotalPTime[j1] > this->m_JobWeightTotalPTime[j2];
            });
        }
    }

    if (SortMethod == 3) //in non-increasing order (job ptime on the first machine)
    {
        for (auto &JobSeq: JobSeqInFam) {
            sort(begin(JobSeq), end(JobSeq), [this](int j1, int j2) {
                return this->m_JobOperPTime[j1][0] > this->m_JobOperPTime[j2][0];
            });
        }
    }

    if (SortMethod == 4) //in non-increasing order (job ptime on the last machine)
    {
        for (auto &JobSeq: JobSeqInFam) {
            sort(begin(JobSeq), end(JobSeq), [this](int j1, int j2) {
                return this->m_JobOperPTime[j1][this->m_Machines - 1] > this->m_JobOperPTime[j2][this->m_Machines - 1];
            });
        }
    }

    if (SortMethod == 5) //in non-increasing order (job ptime on the first machine - job ptime on the last machine)
    {
        for (auto &JobSeq: JobSeqInFam) {
            sort(begin(JobSeq), end(JobSeq), [this](int j1, int j2) {
                return this->m_JobOperPTime[j1][0] - this->m_JobOperPTime[j1][this->m_Machines - 1] >
                       this->m_JobOperPTime[j2][0] - this->m_JobOperPTime[j2][this->m_Machines - 1];
            });
        }
    }

    if (SortMethod == 6) //��ӱӱʦ���������������ɲ�Ʒ�Ĺ�������������
    {
        for (int Fam = 0; Fam < JobSeqInFam.size(); Fam++) {
            Base::Pair<int> *ch = new Base::Pair<int>[JobSeqInFam[Fam].size()];
            vector<int> temp(m_Families, 0);

            for (int j = 0; j < JobSeqInFam[Fam].size(); j++) {
                for (int i = 1; i <= m_Machines; i++)
                    temp[j] = (2 * i - this->m_Machines - 1) * this->m_JobOperPTime[i][JobSeqInFam[Fam][j]];

                ch[j].dim = JobSeqInFam[Fam][j];
                ch[j].value = temp[j];
            }
            sort(ch, ch + JobSeqInFam[Fam].size(), Base::PairGreater<int>());
            for (int j = 0; j < JobSeqInFam[Fam].size(); j++)
                JobSeqInFam[Fam][j] = ch[j].dim;
            delete[] ch;
        }
    }
}



void NOperator::CombineSortedSequences(int SortMethod1, int SortMethod2, vector<int> &FamPermu) {
    vector<int> SortedSeq1, SortedSeq2;

    // ʹ��SortFam��SortedSeq1��SortedSeq2��������
    SortFam(SortMethod1, SortedSeq1);
    SortFam(SortMethod2, SortedSeq2);

    // �ϲ�SortedSeq1��SortedSeq2��ȥ��
    FamPermu = MergeSequences(SortedSeq1, SortedSeq2);
}

vector<int> NOperator::MergeSequences(const vector<int> &SortedSeq1, const vector<int> &SortedSeq2) {
    vector<int> CombinedSeq = SortedSeq1;  // Initialize with the FirstHalf

    // Iterate over the SecondHalf (SortedSeq2) and add unique elements
    for (int i = 0; i < SortedSeq2.size(); ++i) {
        // Check if the current element in SortedSeq2 is already in CombinedSeq
        // If it is not found, add it to CombinedSeq
        if (find(CombinedSeq.begin(), CombinedSeq.end(), SortedSeq2[i]) == CombinedSeq.end()) {
            CombinedSeq.push_back(SortedSeq2[i]);
        }
    }

    // Log output for debugging to show the combined sequence
//    cout << "CombinedSeq: ";
//    for (int elem : CombinedSeq) {
//        cout << elem << " ";
//    }
//    cout << endl;

    return CombinedSeq;
}

/**
 * �������ĳһ��������
 * @param SortMethod
 * @param FamPermu
 */
void NOperator::SortFam(int SortMethod, vector<int> &FamPermu) {

    FamPermu.clear();
    FamPermu.resize(this->m_Families);
    for (int Fam = 0; Fam < this->m_Families; Fam++)
        FamPermu[Fam] = Fam;
    Base::Pair<double> *ch = new Base::Pair<double>[FamPermu.size()];

    if (SortMethod == 0) //LPT
    {
        for (int Fam = 0; Fam < FamPermu.size(); Fam++) {
            ch[Fam].dim = FamPermu[Fam];
            ch[Fam].value = this->m_FamTotalPTime[FamPermu[Fam]];
        }
        sort(ch, ch + this->m_Families, Base::PairGreater<double>());
    }
    if (SortMethod == 1) //SPT
    {
        for (int Fam = 0; Fam < FamPermu.size(); Fam++) {
            ch[Fam].dim = FamPermu[Fam];
            ch[Fam].value = this->m_FamTotalPTime[FamPermu[Fam]];
        }
        sort(ch, ch + this->m_Families, Base::PairLess<double>());
    }

    if (SortMethod == 2) // ���ܼӹ�ʱ��+ƽ�����׼��ʱ��+ƫ��
    {
        for (int Fam = 0; Fam < FamPermu.size(); Fam++) {
            ch[Fam].dim = FamPermu[Fam];
            ch[Fam].value = this->m_FamTotalPTime[FamPermu[Fam]] + this->m_FamAvgSetupTime[FamPermu[Fam]] +
                            this->m_FamTotalSkewness[FamPermu[Fam]];
            //       cout<<"---"<<m_FamTotalPTime[FamPermu[Fam]] << "---"<< m_FamAvgSetupTime[FamPermu[Fam]]<<"---"<<m_FamTotalSkewness[FamPermu[Fam]]<<endl;
        }
        sort(ch, ch + this->m_Families, Base::PairGreater<double>());
    }

    if (SortMethod == 3)// ��Ȩ��+ƽ�����׼��ʱ��+ƫ��
    {
        for (int Fam = 0; Fam < FamPermu.size(); Fam++) {
            ch[Fam].dim = FamPermu[Fam];
            ch[Fam].value = this->m_FamWeightTotalPTime[FamPermu[Fam]] + this->m_FamAvgSetupTime[FamPermu[Fam]] +
                            this->m_FamTotalSkewness[FamPermu[Fam]];
            //     cout<<"---"<<m_FamWeightTotalPTime[FamPermu[Fam]] << "---"<< m_FamAvgSetupTime[FamPermu[Fam]]<<"---"<<m_FamTotalSkewness[FamPermu[Fam]]<<endl;
        }
        sort(ch, ch + this->m_Families,Base:: PairGreater<double>());
    }

    if (SortMethod == 4) {
        for (int Fam = 0; Fam < FamPermu.size(); Fam++) {
            ch[Fam].dim = FamPermu[Fam];
            ch[Fam].value =
                    this->m_FamTotalPTimeOnLastMachine[FamPermu[Fam]] + this->m_FamTotalPTimeOnFirstMachine[FamPermu[Fam]] +
                    this->m_FamTotalSkewness[FamPermu[Fam]];
            //  cout<<"---"<<m_FamTotalPTimeOnLastMachine[FamPermu[Fam]] << "---"<< m_FamAvgSetupTime[FamPermu[Fam]]<<"---"<<m_FamTotalSkewness[FamPermu[Fam]]<<endl;
        }
        sort(ch, ch + this->m_Families, Base::PairGreater<double>());
    }
    if (SortMethod == 5) {
        for (int Fam = 0; Fam < FamPermu.size(); Fam++) {
            ch[Fam].dim = FamPermu[Fam];
            ch[Fam].value =
                    this->m_FamTotalPTimeOnLastMachine[FamPermu[Fam]] - this->m_FamTotalPTimeOnFirstMachine[FamPermu[Fam]] +
                    this->m_FamTotalSkewness[FamPermu[Fam]];
            //     cout<<"---"<<m_FamTotalPTimeOnLastMachine[FamPermu[Fam]]<< "---"<< m_FamAvgSetupTime[FamPermu[Fam]]<<"---"<<m_FamTotalSkewness[FamPermu[Fam]]<<endl;
        }
        sort(ch, ch + this->m_Families, Base::PairGreater<double>());
    }
    for (int Fam = 0; Fam < this->m_Families; Fam++)
        FamPermu[Fam] = ch[Fam].dim;
    delete[]ch;
}


//�� ����ָ����������й����������λ��
float NOperator::FindBestPosToInsertFamForAllFac_Ind_DR_NoAC(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam,

                                                             int InsFam, int& BestFac, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,vector<Individual> &CCEAPopulation)
{
    float minInd = INT_MAX;
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        int Pos = -1;
        float Ind = this->FindBestPosToInsertFamForPerFac_Ind_DR_NoAC(Fac, FacFamSeq, FacFamSeq[Fac], JobSeqInFam,
                                                                      InsFam, Pos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC,CCEAPopulation);
        if (Ind < minInd)
        {
            minInd = Ind;
            BestFac = Fac;
            BestPos = Pos;
        }
    }
    return minInd;
}

//�� ����ָ������ڵ�ǰ���������λ��DR
float NOperator::FindBestPosToInsertFamForPerFac_Ind_DR_NoAC(int Fac, const vector<vector<int>>& FacFamSeq, const vector<int>& NewFamSeqInFac, const vector<vector<int>>& JobSeqInFam,
                                                             int InsFam, int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC,vector<Individual> &CCEAPopulation)
{
    float minFacInd = INT_MAX;

    for (int Pos = 0; Pos <= NewFamSeqInFac.size(); Pos++)
    {

        float FacInd = GetIndForPerFacAfterInsertFam_DR_NoAC(Fac, FacFamSeq, NewFamSeqInFac, JobSeqInFam,
                                                             InsFam, Pos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC,CCEAPopulation);

        if (FacInd < minFacInd)
        {
            minFacInd = FacInd;
            BestPos = Pos;
        }
    }
    return minFacInd;
}


//�� ĳ��λ�ò���������ָ��Ind
float NOperator::GetIndForPerFacAfterInsertFam_DR_NoAC(int Fac, const vector<vector<int>>& FacFamSeq, const vector<int>& FamSeqInFac, const vector<vector<int>>& JobSeqInFam,
                                                       int InsFam, int Pos, int nadirpointMS, float nadirpointTEC, int idealpointMS, float idealpointTEC, vector<Individual> &CCEAPopulation)
{
    //����δ����֮ǰ��������span
    vector <int> FacSpan(FacFamSeq.size(), 0);
    vector<vector<int>> TempFacFamSeq;
    TempFacFamSeq.clear();
    TempFacFamSeq.resize(FacFamSeq.size());
    TempFacFamSeq = FacFamSeq;
    TempFacFamSeq[Fac].insert(TempFacFamSeq[Fac].begin() + Pos, InsFam);

    int Makespan = GetSpan_Forward(TempFacFamSeq, JobSeqInFam,FacSpan);


    //�ж��Ƿ��滻��͵�������
    if (Makespan > nadirpointMS)
        nadirpointMS = Makespan;
    if (Makespan < idealpointMS)
        idealpointMS = Makespan;

    vector<float> FacEC(FacFamSeq.size(), 0);
    float TEC;
    for (int fac = 0; fac < TempFacFamSeq.size(); fac++)
    {
        //�����ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
            int Fam = -1;  // �ȳ�ʼ��
            if (!TempFacFamSeq[fac].empty()) {
                Fam = TempFacFamSeq[fac][0];  // ֻ���ڷǿ�����²Ÿ�ֵ
            } else {
                FacEC[fac] = 0;  // û����Ĺ�����ֱ����Ϊ 0
                continue;
            }
            //  int Fam = TempFacFamSeq[fac][0];
            for (int g = 0; g < TempFacFamSeq[fac].size(); g++)
            {
                Fam = TempFacFamSeq[fac][g];

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

    TEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);
    //cout << "���ܺģ�" << TEC << endl;

    //�ж��Ƿ��滻��͵�������
    if (TEC > nadirpointTEC)
        nadirpointTEC = TEC;
    if (TEC < idealpointTEC)
        idealpointTEC = TEC;

    //��һ��
    float normalMS = -1;
    float normalTEC = -1;
    //cout << endl << "normalize" << endl;

    normalMS = (static_cast<float>(Makespan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    normalTEC = ((TEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
    //cout << "�����һ�����MS��" << normalMS << "\tTEC��" << normalTEC << endl;

    //����ָ��Ind
    float convergence_ind = 0;
    convergence_ind += (normalMS - 1.0) * (normalMS - 1.0);
    convergence_ind += (normalTEC - 1.0) * (normalTEC - 1.0);
    convergence_ind = sqrt(convergence_ind);
    convergence_ind = 1 / (convergence_ind + 1);
    //cout << "���������ָ�꣺" << convergence_ind << endl;
    //cout << endl;
    // ���� calc_distribution_ind ����������ֲ���ָ��
    Individual::calc_distribution_ind(CCEAPopulation);
    float distribution_ind = 0;
    for (int i = 0; i < CCEAPopulation.size(); i++) {
        for (int j = i + 1; j < CCEAPopulation.size(); j++) {
            distribution_ind += CCEAPopulation[i].distribution_ind[j];
        }
    }
    distribution_ind = distribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);
    // �������ռ�Ȩָ��
    float combined_ind = 0.5f * convergence_ind + 0.5f * distribution_ind;

    return combined_ind;
}

float NOperator::FindBestPosToInsertJobForPerFac_Ind_NoAC(int fac ,const vector<vector<int>> &FacFamSeq,const vector<int>& FamSeqInFac, const vector<vector<int>> JobSeqInFam,
                                                          int Fam, int FamPos, int InsJob,
                                                          int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                                          float idealpointTEC,int OrgMS, float OrgTEC,
                                                          vector<Individual> &CCEAPopulation)
{
    float minInd = INT_MAX;
    for (int Pos = 0; Pos <= JobSeqInFam[Fam].size(); Pos++)
    {
        auto re = this->GetIndForPerFacAfterInsertJob_NoAC(fac, FacFamSeq, FamSeqInFac, JobSeqInFam,  Fam, FamPos, InsJob, Pos,
                                                           nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, OrgMS, OrgTEC, CCEAPopulation);

        if (re < minInd)
        {
            minInd = re;
            BestPos = Pos;
        }
    }
    return minInd;
}

//�� ĳ��λ�ò��빤�������ָ��Ind
float NOperator::GetIndForPerFacAfterInsertJob_NoAC(int fac, const vector<vector<int>> &FacFamSeq,
                                                    const vector<int> &FamSeqInFac,
                                                    const vector<vector<int>> &JobSeqInFam,
                                                    int InsFam, int FamPos,
                                                    int InsJob, int JobPos,
                                                    int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                                    float idealpointTEC, int OrgMS, float OrgTEC,
                                                    vector<Individual> &CCEAPopulation) {

    vector<vector<int>> TempFacFamSeq;
    TempFacFamSeq.clear();
    TempFacFamSeq.resize(FacFamSeq.size());
    TempFacFamSeq = FacFamSeq;

    vector<vector<int>> TempJobSeqInFam;
    TempJobSeqInFam.clear();
    TempJobSeqInFam.resize(JobSeqInFam.size());
    TempJobSeqInFam = JobSeqInFam;
    TempJobSeqInFam[InsFam].insert(TempJobSeqInFam[InsFam].begin() + JobPos, InsJob);

    vector <int> FacSpan(FacFamSeq.size(), 0);

    int MakeSpan = GetSpan_Forward(TempFacFamSeq,  TempJobSeqInFam,FacSpan);
    //cout << "�����MakeSpan��" << MakeSpan << endl;

    //�ж��Ƿ��滻��͵�������
    if (MakeSpan > nadirpointMS)
        nadirpointMS = MakeSpan;
    if (MakeSpan < idealpointMS)
        idealpointMS = MakeSpan;


    vector<float> FacEC(FacFamSeq.size(), 0);
    for (int fac = 0; fac < TempFacFamSeq.size(); fac++)
    {
        //�����ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
            //int CurFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
           // int Fam = TempFacFamSeq[fac][0];
            // �жϵ�ǰ�������������Ƿ�Ϊ��
            if ( TempFacFamSeq[fac].empty()) {
                continue; // ֱ��������ǰ����
            }

            int Fam = TempFacFamSeq[fac][0];// ���ڿ��԰�ȫ����
           for (int g = 0; g < TempFacFamSeq[fac].size(); g++)
            {
                Fam = TempFacFamSeq[fac][g];

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
                for (int j = 0; j < TempJobSeqInFam[Fam].size(); j++)
                {
                    int Job = TempJobSeqInFam[Fam][j];
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
    float TEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);
    //cout << "���ܺģ�" << TEC << endl;

    //�ж��Ƿ��滻��͵�������
    if (TEC > nadirpointTEC)
        nadirpointTEC = TEC;
    if (TEC < idealpointTEC)
        idealpointTEC = TEC;
    //�ж��Ƿ�Ȳ���ǰ�Ľ������Ľ������ο���
    bool flag = true;
    if (((MakeSpan < OrgMS) && (TEC < OrgTEC)) || ((MakeSpan < OrgMS) && (TEC == OrgTEC)) ||
        ((MakeSpan == OrgMS) && (TEC < OrgTEC)) || (MakeSpan < idealpointMS) || (TEC < idealpointTEC))
    {
        for (int i = 0; i < CCEAPopulation.size(); i++)
        {
            if ((MakeSpan == CCEAPopulation[i].MS) && (TEC == CCEAPopulation[i].TEC))
            {
                flag = false;
                break;
            }
        }
        if (flag)
        {

            Individual tempIndi;
            tempIndi.m_FacFamSeqArray = TempFacFamSeq;
            tempIndi.m_JobSeqInFamArray = TempJobSeqInFam;
            tempIndi.MS = MakeSpan;
            tempIndi.TEC = TEC;
            tempIndi.m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���
            CCEAPopulation.push_back(tempIndi);
        }
    }

    //��һ��
    float normalMS = -1;
    float normalTEC = -1;
    //cout << endl << "normalize" << endl;

    normalMS = (static_cast<float>(MakeSpan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    normalTEC = ((TEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
    // ���� normalMS
    if (nadirpointMS == idealpointMS) {
        normalMS = 0.0f; // ������Ĭ��ֵ
    } else {
        normalMS = static_cast<float>(MakeSpan - idealpointMS) / (nadirpointMS - idealpointMS);
    }
    // ���� normalTEC
    if (nadirpointTEC == idealpointTEC) {
        normalTEC = 0.0f; // ������Ĭ��ֵ
    } else {
        normalTEC = (TEC - idealpointTEC) / (nadirpointTEC - idealpointTEC);
    }



    //����ָ��Ind
    float convergence_ind = 0;
    convergence_ind += (normalMS - 1.0) * (normalMS - 1.0);
    convergence_ind += (normalTEC - 1.0) * (normalTEC - 1.0);
    convergence_ind = sqrt(convergence_ind);
    convergence_ind = 1 / (convergence_ind + 1);

    Individual::calc_distribution_ind(CCEAPopulation);
    float Distribution_ind = 0;
    for (int i = 0; i < CCEAPopulation.size(); i++) {
        for (int j = i + 1; j < CCEAPopulation.size(); j++) {
            Distribution_ind += CCEAPopulation[i].distribution_ind[j];
        }
    }
    Distribution_ind = Distribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

    //  cout << "Distribution_ind: " << Distribution_ind << endl;`
    // �������ռ�Ȩָ��
    float combined_ind = 0.5f * convergence_ind + 0.5f * Distribution_ind;
    return combined_ind;
}


//�� ĳ��λ�ò��빤�������ָ��Ind
float NOperator::GetIndForPerFacAfterInsertJob_DR_NoAC(int fac, const vector<vector<int>> &FacFamSeq,
                                                       const vector<int> &FamSeqInFac,
                                                       const vector<vector<int>> &JobSeqInFam,
                                                       int InsFam, int FamPos,
                                                       int InsJob, int JobPos,
                                                       int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                                       float idealpointTEC ,vector<Individual> &CCEAPopulation) {


    vector<vector<int>> TempFacFamSeq;
    TempFacFamSeq.clear();
    TempFacFamSeq.resize(FacFamSeq.size());
    TempFacFamSeq = FacFamSeq;

    vector<vector<int>> TempJobSeqInFam;
    TempJobSeqInFam.clear();
    TempJobSeqInFam.resize(JobSeqInFam.size());
    TempJobSeqInFam = JobSeqInFam;
    TempJobSeqInFam[InsFam].insert(TempJobSeqInFam[InsFam].begin() + JobPos, InsJob);

    vector <int> FacSpan(FacFamSeq.size(), 0);

    int MakeSpan = GetSpan_Forward(TempFacFamSeq,  TempJobSeqInFam,FacSpan);
    //cout << "�����MakeSpan��" << MakeSpan << endl;

    //�ж��Ƿ��滻��͵�������
    if (MakeSpan > nadirpointMS)
        nadirpointMS = MakeSpan;
    if (MakeSpan < idealpointMS)
        idealpointMS = MakeSpan;

    vector<float> FacEC(FacFamSeq.size(), 0);
    float TEC;
    for (int fac = 0; fac < TempFacFamSeq.size(); fac++) {
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++) {
            int preFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;

            int Fam = -1;  // �ȳ�ʼ��
            if (!TempFacFamSeq[fac].empty()) {
                Fam = TempFacFamSeq[fac][0];  // ֻ���ڷǿ�����²Ÿ�ֵ
            } else {
                FacEC[fac] = 0;  // û����Ĺ�����ֱ����Ϊ 0
                continue;
            }

            for (int g = 0; g < TempFacFamSeq[fac].size(); g++) {
                Fam = TempFacFamSeq[fac][g];

                if (preFam == -1) {
                    TSetupTimeEC += m_SetupTime[m][Fam][Fam] * UnitSetupEC;
                    Tsetuptime += m_SetupTime[m][Fam][Fam];
                } else {
                    TSetupTimeEC += m_SetupTime[m][preFam][Fam] * UnitSetupEC;
                    Tsetuptime += m_SetupTime[m][preFam][Fam];
                }
                preFam = Fam;

                for (int j = 0; j < TempJobSeqInFam[Fam].size(); j++) {
                    int Job = TempJobSeqInFam[Fam][j];
                    TpTimeEC += m_TureJobOpertime[Job][m] * UnitPEC[Job][m];
                    Tptime += m_TureJobOpertime[Job][m];
                }
            }

            Tidletime += (FacSpan[fac] - Tptime - Tsetuptime);
            TIdleTimeEC += Tidletime * UnitIdleEC;
        }

        FacEC[fac] = TSetupTimeEC + TIdleTimeEC + TpTimeEC;
    }

    TEC = accumulate(FacEC.begin(), FacEC.end(), 0);
    //cout << "���ܺģ�" << TEC << endl;
    //�ж��Ƿ��滻��͵�������
    if (TEC > nadirpointTEC)
        nadirpointTEC = TEC;
    if (TEC < idealpointTEC)
        idealpointTEC = TEC;


    //��һ��
    float normalMS = -1;
    float normalTEC = -1;
    //cout << endl << "normalize" << endl;

    normalMS = (static_cast<float>(MakeSpan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    normalTEC = ((TEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));
    //cout << "�����һ�����MS��" << normalMS << "\tTEC��" << normalTEC << endl;

    //����ָ��Ind
    float convergence_ind = 0;
    convergence_ind += (normalMS - 1.0) * (normalMS - 1.0);
    convergence_ind += (normalTEC - 1.0) * (normalTEC - 1.0);
    convergence_ind = sqrt(convergence_ind);
    convergence_ind = 1 / (convergence_ind + 1);
    //cout << "���������ָ�꣺" << convergence_ind << endl;

    Individual::calc_distribution_ind(CCEAPopulation);
    float distribution_ind = 0;
    for (int i = 0; i < CCEAPopulation.size(); i++) {
        for (int j = i + 1; j < CCEAPopulation.size(); j++) {
            distribution_ind += CCEAPopulation[i].distribution_ind[j];
        }
    }
    distribution_ind = distribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

    // �������ռ�Ȩָ��
    float combined_ind = 0.5f * convergence_ind + 0.5f * distribution_ind;
    return combined_ind;
}