#include <cfloat>
#include "NOperator_New.h"

NOperator_New::NOperator_New() {

}

NOperator_New::~NOperator_New() {

}

//Pareto_relation
void NOperator_New::Pareto_relation(vector<Individual> &Population) {
    int pflag, pflag1;
    for (int i = 0; i < Population.size(); i++) {
        Population[i].flag = 0;
        for (int j = 0; j < Population.size(); j++) {
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
                Population[j].pareto_rel[i] = 1;    //iռ��j
            else
                Population[j].pareto_rel[i] = 0;
        }

    }
}

//Heu
void
NOperator_New::SortJobsInFam(int SortMethod, vector<vector<int>> &JobSeqInFam) //0:LPT,   1:SPT,  2:JobWeightTotalPTime,
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
            Pair<int> *ch = new Pair<int>[JobSeqInFam[Fam].size()];
            vector<int> temp(m_Families, 0);

            for (int j = 0; j < JobSeqInFam[Fam].size(); j++) {
                for (int i = 1; i <= m_Machines; i++)
                    temp[j] = (2 * i - this->m_Machines - 1) * this->m_JobOperPTime[i][JobSeqInFam[Fam][j]];

                ch[j].dim = JobSeqInFam[Fam][j];
                ch[j].value = temp[j];
            }
            sort(ch, ch + JobSeqInFam[Fam].size(), PairGreater<int>());
            for (int j = 0; j < JobSeqInFam[Fam].size(); j++)
                JobSeqInFam[Fam][j] = ch[j].dim;
            delete[] ch;
        }
    }
}



void NOperator_New::CombineSortedSequences(int SortMethod1, int SortMethod2, vector<int> &FamPermu) {
    vector<int> SortedSeq1, SortedSeq2;

    // ʹ��SortFam��SortedSeq1��SortedSeq2��������
    SortFam(SortMethod1, SortedSeq1);
    SortFam(SortMethod2, SortedSeq2);

    // �ϲ�SortedSeq1��SortedSeq2��ȥ��
    FamPermu = MergeSequences(SortedSeq1, SortedSeq2);
}

vector<int> NOperator_New::MergeSequences(const vector<int> &SortedSeq1, const vector<int> &SortedSeq2) {
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
void NOperator_New::SortFam(int SortMethod, vector<int> &FamPermu) {

    FamPermu.clear();
    FamPermu.resize(this->m_Families);
    for (int Fam = 0; Fam < this->m_Families; Fam++)
        FamPermu[Fam] = Fam;
    Pair<double> *ch = new Pair<double>[FamPermu.size()];

    if (SortMethod == 0) //LPT
    {
        for (int Fam = 0; Fam < FamPermu.size(); Fam++) {
            ch[Fam].dim = FamPermu[Fam];
            ch[Fam].value = this->m_FamTotalPTime[FamPermu[Fam]];
        }
        sort(ch, ch + this->m_Families, PairGreater<double>());
    }
    if (SortMethod == 1) //SPT
    {
        for (int Fam = 0; Fam < FamPermu.size(); Fam++) {
            ch[Fam].dim = FamPermu[Fam];
            ch[Fam].value = this->m_FamTotalPTime[FamPermu[Fam]];
        }
        sort(ch, ch + this->m_Families, PairLess<double>());
    }

    if (SortMethod == 2) // ���ܼӹ�ʱ��+ƽ�����׼��ʱ��+ƫ��
    {
        for (int Fam = 0; Fam < FamPermu.size(); Fam++) {
            ch[Fam].dim = FamPermu[Fam];
            ch[Fam].value = this->m_FamTotalPTime[FamPermu[Fam]] + this->m_FamAvgSetupTime[FamPermu[Fam]] +
                            this->m_FamTotalSkewness[FamPermu[Fam]];
                     //       cout<<"---"<<m_FamTotalPTime[FamPermu[Fam]] << "---"<< m_FamAvgSetupTime[FamPermu[Fam]]<<"---"<<m_FamTotalSkewness[FamPermu[Fam]]<<endl;
        }
        sort(ch, ch + this->m_Families, PairGreater<double>());
    }

    if (SortMethod == 3)// ��Ȩ��+ƽ�����׼��ʱ��+ƫ��
    {
        for (int Fam = 0; Fam < FamPermu.size(); Fam++) {
            ch[Fam].dim = FamPermu[Fam];
            ch[Fam].value = this->m_FamWeightTotalPTime[FamPermu[Fam]] + this->m_FamAvgSetupTime[FamPermu[Fam]] +
                            this->m_FamTotalSkewness[FamPermu[Fam]];
      //     cout<<"---"<<m_FamWeightTotalPTime[FamPermu[Fam]] << "---"<< m_FamAvgSetupTime[FamPermu[Fam]]<<"---"<<m_FamTotalSkewness[FamPermu[Fam]]<<endl;
        }
        sort(ch, ch + this->m_Families, PairGreater<double>());
    }

    if (SortMethod == 4) {
        for (int Fam = 0; Fam < FamPermu.size(); Fam++) {
            ch[Fam].dim = FamPermu[Fam];
            ch[Fam].value =
                    this->m_FamTotalPTimeOnLastMachine[FamPermu[Fam]] + this->m_FamTotalPTimeOnFirstMachine[FamPermu[Fam]] +
                    this->m_FamTotalSkewness[FamPermu[Fam]];
          //  cout<<"---"<<m_FamTotalPTimeOnLastMachine[FamPermu[Fam]] << "---"<< m_FamAvgSetupTime[FamPermu[Fam]]<<"---"<<m_FamTotalSkewness[FamPermu[Fam]]<<endl;
        }
        sort(ch, ch + this->m_Families, PairGreater<double>());
    }
    if (SortMethod == 5) {
        for (int Fam = 0; Fam < FamPermu.size(); Fam++) {
            ch[Fam].dim = FamPermu[Fam];
            ch[Fam].value =
                    this->m_FamTotalPTimeOnLastMachine[FamPermu[Fam]] - this->m_FamTotalPTimeOnFirstMachine[FamPermu[Fam]] +
                    this->m_FamTotalSkewness[FamPermu[Fam]];
       //     cout<<"---"<<m_FamTotalPTimeOnLastMachine[FamPermu[Fam]]<< "---"<< m_FamAvgSetupTime[FamPermu[Fam]]<<"---"<<m_FamTotalSkewness[FamPermu[Fam]]<<endl;
        }
        sort(ch, ch + this->m_Families, PairGreater<double>());
    }
    for (int Fam = 0; Fam < this->m_Families; Fam++)
        FamPermu[Fam] = ch[Fam].dim;
    delete[]ch;
}

//NOperator

void NOperator_New::CheckSol(vector<int> FamSeq, vector<vector<int>> JobSeqinFam, int Span) {
    if (!FamSeq.size()) {
        cout << "Error in CheckSol: FamSeq.size()=0" << endl;
        getchar();
        exit(0);
    }

    //Check Families
    vector<bool> bFamArray(this->m_Families, false);

    for (int Fam = 0; Fam < FamSeq.size(); Fam++)
        bFamArray[FamSeq[Fam]] = true;
    for (int i = 0; i < this->m_Families; i++) {
        if (!bFamArray[i]) {
            cout << "Family " << i << " is lost!" << endl;
            getchar();
            exit(0);
        }
    }

    //Check Jobs
    vector<bool> bJobArray(this->m_Jobs, false);
    for (int Fam = 0; Fam < JobSeqinFam.size(); Fam++) {
        for (int j = 0; j < JobSeqinFam[Fam].size(); j++) {
            int Job = JobSeqinFam[Fam][j];
            bJobArray[Job] = true;
        }
    }

    for (int i = 0; i < this->m_Jobs; i++) {
        if (!bJobArray[i]) {
            cout << "Job " << i << " is lost!" << endl;
            getchar();
            exit(0);
        }
    }

    //Check Jobs in each Family
    for (int Fam = 0; Fam < this->m_JobSeqInFam.size(); Fam++) {
        for (int j = 0; j < this->m_JobSeqInFam[Fam].size(); j++) {
            int Job = this->m_JobSeqInFam[Fam][j];
            if (JobSeqinFam[Fam].end() == find(JobSeqinFam[Fam].begin(), JobSeqinFam[Fam].end(), Job)) {
                cout << "Job is not found in Family:\t" << Fam << endl;
                getchar();
                exit(0);
            }
        }
    }

    //----- Check makespan---------
    int TSpan = this->GetSpan(FamSeq, JobSeqinFam);
    if (TSpan != Span) {
        cout << "Span is Erro!" << Span << "\t" << TSpan << endl;
        getchar();
        exit(0);
    }
    cout << "��ȷ��" << endl;
}

//��
void NOperator_New::CheckSol(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam, int Span) {
    if (!FacFamSeq.size()) {
        cout << "Error in CheckSol: FamSeq.size()=0" << endl;
        getchar();
        exit(0);
    }
    vector<bool> bFamArray(this->m_Families, false); //Check Families
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
        for (int Fam = 0; Fam < FacFamSeq[Fac].size(); Fam++)
            bFamArray[FacFamSeq[Fac][Fam]] = true;
    for (int i = 0; i < this->m_Families; i++) {
        if (!bFamArray[i]) {
            cout << "Family " << i << " is lost!" << endl;
            getchar();
            exit(0);
        }
    }
    vector<bool> bJobArray(this->m_Jobs, false); //Check Jobs
    for (int Fam = 0; Fam < JobSeqInFam.size(); Fam++) {
        for (int j = 0; j < JobSeqInFam[Fam].size(); j++) {
            int Job = JobSeqInFam[Fam][j];
            bJobArray[Job] = true;
        }
    }
    for (int i = 0; i < this->m_Jobs; i++) {
        if (!bJobArray[i]) {
            cout << "Job " << i << " is lost!" << endl;
            getchar();
            exit(0);
        }
    }
    for (int Fam = 0; Fam < this->m_JobSeqInFam.size(); Fam++) //Check Jobs in each Family
    {
        for (int j = 0; j < this->m_JobSeqInFam[Fam].size(); j++) {
            int Job = this->m_JobSeqInFam[Fam][j];
            if (JobSeqInFam[Fam].end() == find(JobSeqInFam[Fam].begin(), JobSeqInFam[Fam].end(), Job)) {
                cout << "Job is not found in Family:\t" << Fam << endl;
                getchar();
                exit(0);
            }
        }
    }

    //----- Check makespan---------
    int TSpan = this->GetSpan(FacFamSeq, JobSeqInFam);  //�õ��깤ʱ��
    if (TSpan != Span) {
        cout << "Span is Erro!" << Span << "\t" << TSpan << endl;
        getchar();
        exit(0);
    }

}

//xin
int NOperator_New::GetSpan(vector <vector <int>> FacFamSeq, vector <vector <int>> JobSeqInFam)
{
    vector<vector<int>> JFDTime(this->m_Jobs);
    for (int j = 0; j < JFDTime.size(); j++)
        JFDTime[j].resize(this->m_Machines, 0);
    vector<int> FacSpan(FacFamSeq.size(), 0);
    int mSpan = this->GetJFDTime_Forward(FacFamSeq, JobSeqInFam, JFDTime, FacSpan);   //����GetJCTime_Forward
    return mSpan;
}

void NOperator_New::CheckSolTEC(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam, float TEC) {
    if (!FacFamSeq.size()) {
        cout << "Error in CheckSol: FamSeq.size()=0" << endl;
        getchar();
        exit(0);
    }
    vector<bool> bFamArray(this->m_Families, false); //Check Families
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
        for (int Fam = 0; Fam < FacFamSeq[Fac].size(); Fam++)
            bFamArray[FacFamSeq[Fac][Fam]] = true;
    for (int i = 0; i < this->m_Families; i++) {
        if (!bFamArray[i]) {
            cout << "Family " << i << " is lost!" << endl;
            getchar();
            exit(0);
        }
    }
    vector<bool> bJobArray(this->m_Jobs, false); //Check Jobs
    for (int Fam = 0; Fam < JobSeqInFam.size(); Fam++) {
        for (int j = 0; j < JobSeqInFam[Fam].size(); j++) {
            int Job = JobSeqInFam[Fam][j];
            bJobArray[Job] = true;
        }
    }
    for (int i = 0; i < this->m_Jobs; i++) {
        if (!bJobArray[i]) {
            cout << "Job " << i << " is lost!" << endl;
            getchar();
            exit(0);
        }
    }
    for (int Fam = 0; Fam < this->m_JobSeqInFam.size(); Fam++) //Check Jobs in each Family
    {
        for (int j = 0; j < this->m_JobSeqInFam[Fam].size(); j++) {
            int Job = this->m_JobSeqInFam[Fam][j];
            if (JobSeqInFam[Fam].end() == find(JobSeqInFam[Fam].begin(), JobSeqInFam[Fam].end(), Job)) {
                cout << "Job is not found in Family:\t" << Fam << endl;
                getchar();
                exit(0);
            }
        }
    }

    //----- Check TEC---------
    float tempTEC = this->GetTEC(FacFamSeq, JobSeqInFam);  //�õ��깤ʱ��

    if (tempTEC != TEC) {
        cout.precision(10);
        cout << "TEC is Erro!" << "TEC:" << TEC << "\t" << "tempTEC:" << tempTEC << endl;
        getchar();
        exit(0);
    }
}


int NOperator_New::GetSpan(vector<int> FamSeq, vector<vector<int>> JobSeqinFam) {
    vector<vector<int>> JFDTime(this->m_Jobs);
    for (int j = 0; j < JFDTime.size(); j++)
        JFDTime[j].resize(this->m_Machines, 0);
    int mSpan = GetJFDTime_Forward_InFactory(FamSeq, JobSeqinFam, JFDTime);
    return mSpan;
}


/**
 * ���ݹ����������뿪ʱ����㵥��������MAKESPAN
 * @param FamSeqInFac
 * @param JobSeqInFam
 * @param JFDTime
 * @return
 */
int NOperator_New::GetSpanPerFacByJFD(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                      const vector<vector<int>> &JFDTime) {
    int Span = 0;
    for (int f = 0; f < FamSeqInFac.size(); ++f) {
        int Fam = FamSeqInFac[f];
        int LastJobInFam = JobSeqInFam[Fam][JobSeqInFam[Fam].size() - 1];//������һ������
        Span = JFDTime[LastJobInFam][this->m_Machines - 1];
    }
    return Span;
}

int NOperator_New::GetSpanForAllFacByJFD(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam,
                                           const vector<vector<int>>& JFDTime)
{
    int Makespan = 0;
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        int TempMakespan = GetSpanPerFacByJFD(FacFamSeq[Fac], JobSeqInFam, JFDTime);
        if (Makespan < TempMakespan)
        {
            Makespan = TempMakespan;
        }
    }
    return Makespan;
}


//�� �õ����ܺ�
float NOperator_New::GetTEC(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam) {
    double TEC = 0;
    vector<vector<int>> JFDTime(this->m_Jobs);
    for (int j = 0; j < JFDTime.size(); j++)
        JFDTime[j].resize(this->m_Machines, 0);
    vector<float> FacEC(FacFamSeq.size(), 0);
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++) {
        FacEC[Fac] = GetPerFacEC(FacFamSeq[Fac], JobSeqInFam, JFDTime);//�õ�ÿ��������EC
        TEC += FacEC[Fac];
    }
    return TEC;
}

float
NOperator_New::GetPerFacEC(vector<int> FamSeqInFac, vector<vector<int>> JobSeqInFam, vector<vector<int>> &JFDTime) {

    int Span = GetJFDTime_Forward_InFactory(FamSeqInFac, JobSeqInFam, JFDTime);//�õ�ÿ���������깤ʱ��

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
        Tidletime += (Span - Tptime - Tsetuptime);
        TIdleTimeEC += Tidletime * UnitIdleEC;
    }
    float FacEC = TSetupTimeEC + TIdleTimeEC + TpTimeEC;
    return FacEC;
}

int NOperator_New::GetDelayTime_Forward(vector<vector<int>> FacFamSeq, vector<vector<int>> JobSeqInFam,
                                        vector<vector<int>> &JFDTime, vector<int> &FacSpan,
                                        vector<vector<int>> &DelayTime)//Forward pass calculation
{
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
        FacSpan[Fac] = GetDelayTime_Forward_InFactory(FacFamSeq[Fac], JobSeqInFam, JFDTime, DelayTime);//�õ�ÿ���������깤ʱ��
    return *max_element(FacSpan.begin(), FacSpan.end()); //�õ��ؼ�����
}


int NOperator_New::GetDelayTime_Forward_InFactory(vector<int> FamSeqInFac, vector<vector<int>> JobSeqinFam,
                                                  vector<vector<int>> &JFDTime, vector<vector<int>> &DelayTime) {
    for (int Fam = 0; Fam < FamSeqInFac.size(); Fam++) {
        int CurFam = FamSeqInFac[Fam];  // ��ǰ��
        if (Fam == 0)  // ��һ����
        {
            for (int m = 0; m < this->m_Machines; m++)
                MachReadyTime[m] = this->m_SetupTime[m][CurFam][CurFam];
        } else  // �ӵڶ��鵽�����
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

            for (int m = 1; m < this->m_Machines; m++) {
                // �������һ̨������ֱ�Ӽ����깤ʱ��
                if (m == this->m_Machines - 1) {
                    JFDTime[CurJob][m] = JFDTime[CurJob][m - 1] + this->m_TureJobOpertime[CurJob][m];
                } else {
                    // �����������Ͽ������������
                    JFDTime[CurJob][m] = max(JFDTime[CurJob][m - 1] + m_TureJobOpertime[CurJob][m],
                                             MachReadyTime[m + 1]);
                }

                // ��¼ DelayTime[PreJob][m] �ڴ˴�
                if (j > 0)  // ������ǵ�һ������
                {
                    int PreJob = JobSeqinFam[CurFam][j - 1];  // ǰһ������
                    // ���ǰһ���������깤ʱ�����ڻ�����׼��ʱ�䣬���¼�ӳ�ʱ��
                    if (JFDTime[CurJob][m - 1] > MachReadyTime[m]) {
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
 * ���ݹ�����ǰ���뿪ʱ��������ܺ�TEC
 * @param FamSeqInFac
 * @param JobSeqInFam
 * @param JFDTime
 * @return
 */
float NOperator_New::GetTECForAllFacByJFD(const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                          const vector<vector<int>> &JFDTime) {
    float TEC = 0;
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++) {
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
float NOperator_New::GetTECForPerFacByJFD(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                          const vector<vector<int>> &JFDTime) {

    if (FamSeqInFac.empty()) {
        return 0; // ���߷���һ���ʵ���ֵ
    }
    int MakeSpan = GetSpanPerFacByJFD(FamSeqInFac,JobSeqInFam,JFDTime);

    //�����ܺ�
    float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

    for (int m = 0; m < m_Machines; m++) {
        int preFam = -1;
        int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
        int Fam = FamSeqInFac[0];
        for (int g = 0; g < FamSeqInFac.size(); g++) {
            Fam = FamSeqInFac[g];

            if (preFam == -1) {
                TSetupTimeEC += m_SetupTime[m][Fam][Fam] * UnitSetupEC;
                Tsetuptime += m_SetupTime[m][Fam][Fam];
            } else {
                TSetupTimeEC += m_SetupTime[m][preFam][Fam] * UnitSetupEC;
                Tsetuptime += m_SetupTime[m][preFam][Fam];
            }
            preFam = Fam;
            for (int j = 0; j < JobSeqInFam[Fam].size(); j++) {
                int Job = JobSeqInFam[Fam][j];
                TpTimeEC += m_TureJobOpertime[Job][m]* UnitPEC[Job][m];
                Tptime += m_TureJobOpertime[Job][m];
            }

        }
        Tidletime += (MakeSpan - Tptime - Tsetuptime);
        TIdleTimeEC += Tidletime * UnitIdleEC;
    }
    float FacEC = TSetupTimeEC + TIdleTimeEC + TpTimeEC;
    return FacEC;
}


void NOperator_New::GetMSandTECForPerandToalFacByJFD(const vector<vector<int>>& FacFamSeq, const vector<vector<int>>& JobSeqInFam,
                                                     const vector<vector<int>>& JFDTime, vector<int>& FacSpan, vector<float>& FacEC, int& Makespan, float& TotalEC)
{

    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        FacSpan[Fac]= GetSpanPerFacByJFD(FacFamSeq[Fac],JobSeqInFam,JFDTime);
        FacEC[Fac]= GetTECForPerFacByJFD(FacFamSeq[Fac],JobSeqInFam,JFDTime);
    }
    Makespan = *max_element(FacSpan.begin(), FacSpan.end()); //�õ��ؼ�����
    TotalEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);
   // cout<<"Get:"<<Makespan<<"+"<<TotalEC<< endl;

}

/**
 * ����ÿ�������������뿪ʱ��
 * @param FamSeq
 * @param JobSeqinFam
 * @param JFDTime
 * @param Span
 */
int NOperator_New::GetJFDTime_Forward(const vector<vector<int>> &FamSeq, const vector<vector<int>> &JobSeqInFam,
                                      vector<vector<int>> &JFDTime, vector<int> &FacSpan ) {
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
NOperator_New::GetJFDTime_Forward_InFactory(const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                            vector<vector<int>> &JFDTime) {

    for (int Fam = 0; Fam < FamSeqInFac.size(); Fam++) {
        int CurFam = FamSeqInFac[Fam];//��ǰ��
        if (Fam == 0)//the first group of jobs ��һ����
        {
            for (int m = 0; m < this->m_Machines; m++)
                MachReadyTime[m] = this->m_SetupTime[m][CurFam][CurFam];
        } else //�ڶ�������� from the second group of jobs to the end;
        {
            int PreFam = FamSeqInFac[Fam - 1];
            for (int m = 0; m < this->m_Machines; m++)
                MachReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
        }

        for (int j = 0; j < JobSeqInFam[CurFam].size(); j++)//��ǰ��Ĺ������� Scheduling Jobs in CurFam
        {
            int CurJob = JobSeqInFam[CurFam][j];//��ǰ����
            JFDTime[CurJob][0] = max(MachReadyTime[0] + this->m_TureJobOpertime[CurJob][0], MachReadyTime[1]);

            for (int m = 1; m < this->m_Machines; m++) {
                if (m == this->m_Machines - 1) {
                    JFDTime[CurJob][m] = JFDTime[CurJob][m - 1] + this->m_TureJobOpertime[CurJob][m];
                } else {
                    JFDTime[CurJob][m] = max(JFDTime[CurJob][m - 1] + m_TureJobOpertime[CurJob][m],
                                             MachReadyTime[m + 1]);
                   // cout << JFDTime[CurJob][m] << endl;
                }
            }
            MachReadyTime = JFDTime[CurJob];//����
        }
    }
    return MachReadyTime[this->m_Machines - 1];
}


/**
 * �¿���
 * ����һ������
 * ��ͨ����
 * @param job
 * @param machine_ready ����׼��ʱ��
 */
void NOperator_New::scheduling_one_job(int job, vector<int> &machine_ready) {

    machine_ready[0] += m_TureJobOpertime[job][0];

    for (int m = 1; m < m_Machines; ++m) {
        if (machine_ready[m - 1] >= machine_ready[m]) {
            machine_ready[m] = machine_ready[m - 1] + m_TureJobOpertime[job][m];
        } else {
            machine_ready[m - 1] = machine_ready[m];
            machine_ready[m] = machine_ready[m] + m_TureJobOpertime[job][m];
        }
    }

}

/**
 * �¿���
 * ����һ������
 * ��¼��κ��뿪ʱ��
 * @param job
 * @param machine_ready
 * @param one_job_hierarchy
 */
void NOperator_New::scheduling_one_job(int job, vector<int> &machine_ready, vector<int> &one_job_hierarchy) {

        int job_start_time = machine_ready[0];
        machine_ready[0] += m_TureJobOpertime[job][0];

        for (int m = 1; m < m_Machines; ++m)
        {
            if (machine_ready[m - 1] >= machine_ready[m])
            {
                one_job_hierarchy[m - 1] = machine_ready[m - 1] - job_start_time;
                job_start_time = machine_ready[m - 1];
                machine_ready[m] = machine_ready[m - 1] + m_TureJobOpertime[job][m];
            }
            else
            {
                one_job_hierarchy[m - 1] = machine_ready[m] - job_start_time;
                job_start_time = machine_ready[m];
                machine_ready[m - 1] = machine_ready[m];
                machine_ready[m] = machine_ready[m] + m_TureJobOpertime[job][m];
            }
        }
}

bool
NOperator_New::scheduling_one_job_with_Update_Return(int job, vector<int> &machine_ready,
                                                     vector<int> &one_job_hierarchy) {

    bool is_job_hierarchy_same = true;

    int job_start_time = machine_ready[0];

    machine_ready[0] += m_TureJobOpertime[job][0];

    for (int m = 1; m < m_Machines; ++m) {
        if (machine_ready[m - 1] >= machine_ready[m]) {
            if (machine_ready[m - 1] - job_start_time != one_job_hierarchy[m - 1]) {
                is_job_hierarchy_same = false;
                one_job_hierarchy[m - 1] = machine_ready[m - 1] - job_start_time;
            }
            job_start_time = machine_ready[m - 1];
            machine_ready[m] = machine_ready[m - 1] + m_TureJobOpertime[job][m];
        } else {
            if (machine_ready[m] - job_start_time != one_job_hierarchy[m - 1]) {
                is_job_hierarchy_same = false;
                one_job_hierarchy[m - 1] = machine_ready[m] - job_start_time;
            }
            job_start_time = machine_ready[m];
            machine_ready[m - 1] = machine_ready[m];
            machine_ready[m] = machine_ready[m] + m_TureJobOpertime[job][m];
        }
    }
    return is_job_hierarchy_same;
}

/**
 * �¿���
 * ����һ������
 * �жϲ���Ƿ��ԭ����ͬ
 * @param job
 * @param machine_ready
 * @param one_job_hierarchy
 * @param scenario
 * @return
 */
bool
NOperator_New::scheduling_one_job(int job, vector<int> &machine_ready, const vector<int> &one_job_hierarchy) {

    bool is_job_hierarchy_same = true;

    int job_start_time = machine_ready[0];

    machine_ready[0] += m_TureJobOpertime[job][0];

    for (int m = 1; m < m_Machines; ++m) {
        if (machine_ready[m - 1] >= machine_ready[m]) {
            if (machine_ready[m - 1] - job_start_time != one_job_hierarchy[m - 1]) {
                is_job_hierarchy_same = false;
            }
            job_start_time = machine_ready[m - 1];
            machine_ready[m] = machine_ready[m - 1] + m_TureJobOpertime[job][m];
        } else {
            if (machine_ready[m] - job_start_time != one_job_hierarchy[m - 1]) {
                is_job_hierarchy_same = false;
            }
            job_start_time = machine_ready[m];
            machine_ready[m - 1] = machine_ready[m];
            machine_ready[m] = machine_ready[m] + m_TureJobOpertime[job][m];
        }
    }
    return is_job_hierarchy_same;
}

/**
 * �¿���
 * �ڲ�������߹���֮��
 * ���±仯�Ĺ����Ĳ�κ��뿪ʱ�䣨�����ڲ�����򹤼��������
 * @param FamSeqInFac
 * @param JobSeqInFam
 * @param JFDTime
 * @param hierarchy
 * @param StartFamPos
 * @param StartJobPos
 * @param scenario
 */
void NOperator_New::RefreshJFDTimeHierarchy_InFactory_Insert(const vector<int> &FamSeqInFac,
                                                             const vector<vector<int>> &JobSeqInFam,
                                                             vector<vector<int>> &JFDTime,
                                                             vector<vector<int>> &hierarchy,
                                                             int StartFamPos, int StartJobPos, int EndJobPos) {
    if (FamSeqInFac.empty() || StartFamPos >= FamSeqInFac.size()) {
        return;
    }

    int StartFam = FamSeqInFac[StartFamPos];

    if (StartJobPos == 0) {
        //�����ڵ�һ��λ��
        if (StartFamPos == 0) {
            //��һ����ĵ�һ��λ��
            for (int m = 0; m < this->m_Machines; m++) {
                //ReadyTime��Ϊ׼��ʱ��
                MachReadyTime[m] = this->m_SetupTime[m][StartFam][StartFam];
            }
        } else {
            //������ĵ�һ��λ��
            int PreFam = FamSeqInFac[StartFamPos - 1];//ǰһ����
            int LastJobInPreFam = *JobSeqInFam[PreFam].rbegin();//ǰһ��������һ������
            for (int m = 0; m < this->m_Machines; ++m) {
                //ǰһ�������һ���������뿪ʱ��+׼��ʱ��
                MachReadyTime[m] = JFDTime[LastJobInPreFam][m] + this->m_SetupTime[m][PreFam][StartFam];
            }
        }
    } else {
        //�������ڵ�һ��λ��
        int PreJob = JobSeqInFam[StartFam][StartJobPos - 1]; //ǰһ������
        for (int m = 0; m < this->m_Machines; m++) {
            MachReadyTime[m] = JFDTime[PreJob][m];
        }
    }

    //������������빤�����֣���һ���ֲ������ÿ�������Ϊ������֮ǰ�����ͬ��������������β���
    //vector<vector<int>> temp_hierarchy(this->m_Jobs, vector<int>(this->m_Machines - 1)); //��ʱ�Ĳ�� //wangyuting 0330

    for (int j = StartJobPos; j <= EndJobPos; j++) {
        int CurJob = JobSeqInFam[StartFam][j];
        scheduling_one_job(CurJob, MachReadyTime, hierarchy[CurJob]);
        //hierarchy[CurJob] = temp_hierarchy[CurJob]; //wangyuting 0330
        JFDTime[CurJob] = MachReadyTime;
    }

    //����ʼ�仯�Ĺ���
    bool is_hierarchy_same = false; //����Ƿ���ͬ
    int offset = 0; //ƫ����
    for (int j = EndJobPos + 1; j < JobSeqInFam[StartFam].size(); j++) {
        int CurJob = JobSeqInFam[StartFam][j];

        //��ǰ��β�ͬ������ͨ�������п��ܻ���ͬ
        is_hierarchy_same = scheduling_one_job_with_Update_Return(CurJob, MachReadyTime, hierarchy[CurJob]);
        //scheduling_one_job(CurJob, MachReadyTime, temp_hierarchy[CurJob], scenario);
        if (is_hierarchy_same) {
            offset = MachReadyTime[this->m_Machines - 1] - JFDTime[CurJob][this->m_Machines - 1];
            for (int j1 = j; j1 < JobSeqInFam[StartFam].size(); ++j1) {
                CurJob = JobSeqInFam[StartFam][j1];
                for (int m = 0; m < this->m_Machines; ++m) {
                    JFDTime[CurJob][m] += offset;
                }
            }
            break;
        } else {
            //����б仯
            //hierarchy[CurJob] = temp_hierarchy[CurJob];
            JFDTime[CurJob] = MachReadyTime;
        }
    }

    //������һ����֮�����Ͳ��ü���ÿ��������
    if (is_hierarchy_same) {
        for (int Fam = StartFamPos + 1; Fam < FamSeqInFac.size(); ++Fam) {
            int CurFam = FamSeqInFac[Fam];
            for (int job: JobSeqInFam[CurFam]) {
                //��ǰ����
                for (int m = 0; m < this->m_Machines; m++) {
                    JFDTime[job][m] += offset;
                }
            }
        }
    } else {
        //Ŀǰ��λ���һ�������������ܻ�һ��
        int PreFam = StartFam;
        for (int Fam = StartFamPos + 1; Fam < FamSeqInFac.size(); ++Fam) {
            int CurFam = FamSeqInFac[Fam];
            for (int m = 0; m < this->m_Machines; m++) {
                MachReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
            }
            //�������ڵĹ���
            for (int j = 0; j < JobSeqInFam[CurFam].size(); ++j) {
                int CurJob = JobSeqInFam[CurFam][j];

                //��β�ͬ����ͨ�����ȹ�������ͬ�Ŀ���
                is_hierarchy_same = scheduling_one_job_with_Update_Return(CurJob, MachReadyTime, hierarchy[CurJob]);
                //scheduling_one_job(CurJob, MachReadyTime, temp_hierarchy[CurJob], scenario);
                if (is_hierarchy_same) {
                    offset = MachReadyTime[this->m_Machines - 1] - JFDTime[CurJob][this->m_Machines - 1];
                    for (int j1 = j; j1 < JobSeqInFam[CurFam].size(); ++j1) {
                        CurJob = JobSeqInFam[CurFam][j1];
                        for (int m = 0; m < this->m_Machines; ++m) {
                            JFDTime[CurJob][m] += offset;
                        }
                    }
                    for (int Fam1 = Fam + 1; Fam1 < FamSeqInFac.size(); ++Fam1) {
                        CurFam = FamSeqInFac[Fam1];
                        for (int job: JobSeqInFam[CurFam]) {
                            for (int m = 0; m < this->m_Machines; ++m) {
                                JFDTime[job][m] += offset;
                            }
                        }
                    }
                    return;
                } else {
                    //����б仯
                    //hierarchy[CurJob] = temp_hierarchy[CurJob];
                    JFDTime[CurJob] = MachReadyTime;
                }
            }
            PreFam = CurFam;
        }
    }
}

/**
 * �¿���
 * ��Ĩȥ����߹���֮��
 * ���±仯�Ĺ����Ĳ�κ��뿪ʱ�䣨������Ĩȥ��򹤼��������
 * @param FamSeqInFac
 * @param JobSeqInFam
 * @param JFDTime
 * @param hierarchy
 * @param StartFamPos
 * @param StartJobPos
 * @param scenario
 */
void
NOperator_New::RefreshJFDTimeHierarchy_InFactory_Erase(const vector<int> &FamSeqInFac,
                                                       const vector<vector<int>> &JobSeqInFam,
                                                       vector<vector<int>> &JFDTime, vector<vector<int>> &hierarchy,
                                                       int StartFamPos, int StartJobPos) {
    if (FamSeqInFac.empty() || StartFamPos >= FamSeqInFac.size()) {
        return;
    }

    //vector<int> MachReadyTime(this->m_Machines, 0);//����׼��ʱ���ʼ��Ϊ0

    int StartFam = FamSeqInFac[StartFamPos];
    //int StartJob = JobSeqInFam[StartFam][StartJobPos];

    if (StartJobPos == 0) {
        //�����ڵ�һ��λ��
        if (StartFamPos == 0) {
            //��һ����ĵ�һ��λ��
            for (int m = 0; m < this->m_Machines; m++) {
                //ReadyTime��Ϊ׼��ʱ��
                MachReadyTime[m] = this->m_SetupTime[m][StartFam][StartFam];
            }
        } else {
            //������ĵ�һ��λ��
            int PreFam = FamSeqInFac[StartFamPos - 1];//ǰһ����
            int LastJobInPreFam = *JobSeqInFam[PreFam].rbegin();//ǰһ��������һ������
            for (int m = 0; m < this->m_Machines; ++m) {
                //ǰһ�������һ���������뿪ʱ��+׼��ʱ��
                MachReadyTime[m] = JFDTime[LastJobInPreFam][m] + this->m_SetupTime[m][PreFam][StartFam];
            }
        }
    } else {
        //�������ڵ�һ��λ��
        int PreJob = JobSeqInFam[StartFam][StartJobPos - 1]; //ǰһ������
        for (int m = 0; m < this->m_Machines; m++) {
            MachReadyTime[m] = JFDTime[PreJob][m];
        }
    }

    //����ʼ�仯�Ĺ���
    bool is_hierarchy_same = false; //����Ƿ���ͬ
    int offset = 0; //ƫ����
    //vector<vector<int>> temp_hierarchy(this->m_Jobs, vector<int>(this->m_Machines - 1)); //��ʱ�Ĳ�Σ��жϺ�֮ǰ�Ĳ���Ƿ�һ��
    for (int j = StartJobPos; j < JobSeqInFam[StartFam].size(); ++j) {
        int CurJob = JobSeqInFam[StartFam][j];

        //��ǰ��β�ͬ������ͨ�������п��ܻ���ͬ
        is_hierarchy_same = scheduling_one_job_with_Update_Return(CurJob, MachReadyTime, hierarchy[CurJob]);
        //scheduling_one_job(CurJob, MachReadyTime, temp_hierarchy[CurJob], scenario);
        if (is_hierarchy_same) {
            offset = MachReadyTime[this->m_Machines - 1] - JFDTime[CurJob][this->m_Machines - 1];
            for (int j1 = j; j1 < JobSeqInFam[StartFam].size(); ++j1) {
                CurJob = JobSeqInFam[StartFam][j1];
                for (int m = 0; m < this->m_Machines; ++m) {
                    JFDTime[CurJob][m] += offset;
                }
            }
            break;
        } else {
            //����б仯
            //hierarchy[CurJob] = temp_hierarchy[CurJob];
            JFDTime[CurJob] = MachReadyTime;
        }
    }

    //������һ����֮�����Ͳ��ü���ÿ��������
    if (is_hierarchy_same) {
        for (int Fam = StartFamPos + 1; Fam < FamSeqInFac.size(); ++Fam) {
            int CurFam = FamSeqInFac[Fam];
            for (int j = 0; j < JobSeqInFam[CurFam].size(); ++j) {
                int CurJob = JobSeqInFam[CurFam][j]; //��ǰ����
                for (int m = 0; m < this->m_Machines; m++) {
                    JFDTime[CurJob][m] += offset;
                }
            }
        }
    } else {
        //Ŀǰ��λ���һ�������������ܻ�һ��
        int PreFam = StartFam;
        for (int Fam = StartFamPos + 1; Fam < FamSeqInFac.size(); ++Fam) {
            int CurFam = FamSeqInFac[Fam];
            for (int m = 0; m < this->m_Machines; m++) {
                MachReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
            }
            //�������ڵĹ���
            for (int j = 0; j < JobSeqInFam[CurFam].size(); ++j) {
                int CurJob = JobSeqInFam[CurFam][j];
                //��β�ͬ����ͨ�����ȹ�������ͬ�Ŀ���
                is_hierarchy_same = scheduling_one_job_with_Update_Return(CurJob, MachReadyTime, hierarchy[CurJob]);
                //scheduling_one_job(CurJob, MachReadyTime, temp_hierarchy[CurJob], scenario);
                if (is_hierarchy_same) {
                    offset = MachReadyTime[this->m_Machines - 1] - JFDTime[CurJob][this->m_Machines - 1];
                    for (int j1 = j; j1 < JobSeqInFam[CurFam].size(); ++j1) {
                        CurJob = JobSeqInFam[CurFam][j1];
                        for (int m = 0; m < this->m_Machines; ++m) {
                            JFDTime[CurJob][m] += offset;
                        }
                    }
                    for (int Fam1 = Fam + 1; Fam1 < FamSeqInFac.size(); ++Fam1) {
                        CurFam = FamSeqInFac[Fam1];
                        for (int job: JobSeqInFam[CurFam]) {
                            for (int m = 0; m < this->m_Machines; ++m) {
                                JFDTime[job][m] += offset;
                            }
                        }
                    }
                    return;
                } else {
                    //����б仯
                    //hierarchy[CurJob] = temp_hierarchy[CurJob];
                    JFDTime[CurJob] = MachReadyTime;
                }
            }
            PreFam = CurFam;
        }
    }
}


/**
 * �¿���
 * �����й������ҵ�������makespan��С��λ��
 * @param FacFamSeq
 * @param JobSeqInFam
 * @param JFDTime
 * @param InsFam
 * @param BestFac
 * @param BestPos
 * @return
 */
double NOperator_New::FindBestPosToInsertFamForAllFacs_MC_New(const vector<vector<int>> &FacFamSeq,
                                                              const vector<vector<int>> &JobSeqInFam,
                                                              const vector<vector<int>> &JFDTime,
                                                              const vector<vector<int>> &Hierarchy,
                                                              int InsFam, int &BestFac, int &BestPos) {
    int minMC = INT_MAX;
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++) {
        int Pos;
        int MC = FindBestPosToInsertFam_InFactory_MC_New(Fac, FacFamSeq, FacFamSeq[Fac], JobSeqInFam, JFDTime, Hierarchy, InsFam,Pos);
        if (MC < minMC) {
            minMC = MC;
            BestFac = Fac;
            BestPos = Pos;
        }
    }
    return minMC;
}


double NOperator_New::FindBestPosToInsertFam_InFactory_MC_New(int Fac , const vector<vector<int>> &FacFamSeq,const vector<int> &NewFamSeqInFac,
                                                              const vector<vector<int>> &JobSeqInFam,
                                                              const vector<vector<int>> &JFDTime,
                                                              const vector<vector<int>> &Hierarchy, int InsFam,
                                                              int &BestPos) {
    int minMC = INT_MAX;
    int MC = 0;
    for (int Pos = 0; Pos <= NewFamSeqInFac.size(); ++Pos) {
        auto re = GetMCForPerFacAfterInsertFam_New(Fac,FacFamSeq, NewFamSeqInFac, JobSeqInFam, JFDTime, Hierarchy, InsFam, Pos);
        MC = re;
        if (MC < minMC) {
            minMC = MC;
            BestPos = Pos;
        }
    }

    return minMC;
}

double
NOperator_New::GetMCForPerFacAfterInsertFam_New(int Fac , const vector<vector<int>> &FacFamSeq,const vector<int> &FamSeqInFac, const vector<vector<int>> &JobSeqInFam,
                                                const vector<vector<int>> &JFDTime,
                                                const vector<vector<int>> &Hierarchy,
                                                int InsFam, int Pos) {

    //����δ����֮ǰ��������span
    vector <int> FacSpan(FacFamSeq.size(), 0);
    int Span = GetSpanForPerFacAfterInsertFam_New(FamSeqInFac,JobSeqInFam,JFDTime,Hierarchy,InsFam,Pos);
    for (int f = 0; f < FacFamSeq.size(); f++) {
        FacSpan[f]= GetSpanPerFacByJFD(FacFamSeq[f],JobSeqInFam,JFDTime);
    }
    //�������󹤳�Fac��span
    FacSpan[Fac]= Span;
    int Makespan = *max_element(FacSpan.begin(), FacSpan.end());

    vector <int> FacEC(FacFamSeq.size(), 0);
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
    }

    float TEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);

    double MC = 0.5f * Makespan + 0.5f * TEC ;
    return MC;
}

float NOperator_New::FindBestPosInsertFamAllFacs_TEC_New(const vector<vector<int>> &FacFamSeq,
                                                         const vector<vector<int>> &JobSeqInFam,
                                                         const vector<vector<int>> &JFDTime,
                                                         const vector<vector<int>> &Hierarchy, int InsFam, int &BestFac,
                                                         int &BestPos) {
    float minTEC = INT_MAX;

    //����δ����֮ǰ��������EC
    vector<float> FacEC(FacFamSeq.size(), 0.0);
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++) {
        FacEC[Fac] = GetTECForPerFacByJFD(FacFamSeq[Fac], JobSeqInFam, JFDTime);

    }

    vector<float> TempFacEC(FacFamSeq.size(), 0.0);
    TempFacEC = FacEC;

    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++) {
        int Pos = -1;
        TempFacEC[Fac] = this->FindBestPosToInsertFam_InFactory_TEC_New(FacFamSeq[Fac], FacFamSeq,JobSeqInFam, JFDTime, Hierarchy,
                                                                        InsFam, Pos);

        float TEC = accumulate(TempFacEC.begin(), TempFacEC.end(), 0.0);

        if (TEC < minTEC) {
            minTEC = TEC;
            BestFac = Fac;
            BestPos = Pos;
        }
        TempFacEC[Fac] = FacEC[Fac];
    }
    return minTEC;
}

float NOperator_New::FindBestPosToInsertFam_InFactory_TEC_New(const vector<int> &NewFamSeqInFac,const vector<vector<int>> &FacFamSeq,
                                                              const vector<vector<int>> &JobSeqInFam,
                                                              const vector<vector<int>> &JFDTime,
                                                              const vector<vector<int>> &Hierarchy, int InsFam,
                                                              int &BestPos) {
    float minFacEC = INT_MAX;

    for (int Pos = 0; Pos <= NewFamSeqInFac.size(); Pos++) {
        float FacEC = GetTECForPerFacAfterInsertFam_New(NewFamSeqInFac, FacFamSeq,JobSeqInFam, JFDTime, Hierarchy, InsFam, Pos);

        if (FacEC < minFacEC) {
            minFacEC = FacEC;
            BestPos = Pos;
        }

    }
    return minFacEC;
}


float
NOperator_New::GetTECForPerFacAfterInsertJob_New(const vector<int> &FamSeqInFac,  const vector<vector<int>> &FacFamSeq, vector<vector<int>> &JobSeqInFam,
                                                 const vector<vector<int>> &JFDTime,
                                                 const vector<vector<int>> &hierarchy,
                                                 int InsFam, int FamPos,int InsJob,int JobPos) {
    int Span = 0;
    Span= GetSpanForPerFacAfterInsertJob_New(FamSeqInFac,JobSeqInFam,JFDTime,hierarchy,InsFam,FamPos,InsJob,JobPos);
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

    float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

    for (int f = 0; f < TempFacFamSeq.size(); ++f) {
        //�����ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++) {
            int preFam = -1;
            int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
            int Fam = FamSeqInFac[0];
            for (int g = 0; g < FamSeqInFac.size(); g++) {
                Fam = FamSeqInFac[g];

                if (preFam == -1) {
                    TSetupTimeEC += m_SetupTime[m][Fam][Fam] * UnitSetupEC;
                    Tsetuptime += m_SetupTime[m][Fam][Fam];
                } else {
                    TSetupTimeEC += m_SetupTime[m][preFam][Fam] * UnitSetupEC;
                    Tsetuptime += m_SetupTime[m][preFam][Fam];
                }
                preFam = Fam;
                for (int j = 0; j < TempJobSeqInFam[Fam].size(); j++) {
                    int Job =  TempJobSeqInFam[Fam][j];
                    TpTimeEC += m_TureJobOpertime[Job][m]* UnitPEC[Job][m];
                    Tptime += m_TureJobOpertime[Job][m];
                }

            }
            Tidletime += (Span - Tptime - Tsetuptime);
            TIdleTimeEC += Tidletime * UnitIdleEC;
        }
    }
    float FacEC = TSetupTimeEC + TIdleTimeEC + TpTimeEC;
    return FacEC;
}

float
NOperator_New::GetTECForPerFacAfterInsertFam_New(const vector<int> &FamSeqInFac,  const vector<vector<int>> &FacFamSeq, const vector<vector<int>> &JobSeqInFam,
                                                 const vector<vector<int>> &JFDTime,
                                                 const vector<vector<int>> &Hierarchy,
                                                 int InsFam, int FamPos) {
    int Span = 0;
    Span= GetSpanForPerFacAfterInsertFam_New(FamSeqInFac,JobSeqInFam,JFDTime,Hierarchy,InsFam, FamPos);
    //�����ܺ�
    float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

    vector<int> TempFamSeqInFac;
    TempFamSeqInFac.clear();
    TempFamSeqInFac.resize(FacFamSeq.size());
    TempFamSeqInFac = FamSeqInFac;
    TempFamSeqInFac.insert(TempFamSeqInFac.begin() + FamPos, InsFam);

    for (int m = 0; m < m_Machines; m++) {
        int preFam = -1;
        int Tptime = 0, Tsetuptime = 0, Tidletime = 0;
        int Fam = TempFamSeqInFac[0];
        for (int g = 0; g < TempFamSeqInFac.size(); g++) {
            Fam = TempFamSeqInFac[g];

            if (preFam == -1) {
                TSetupTimeEC += m_SetupTime[m][Fam][Fam] * UnitSetupEC;
                Tsetuptime += m_SetupTime[m][Fam][Fam];
            } else {
                TSetupTimeEC += m_SetupTime[m][preFam][Fam] * UnitSetupEC;
                Tsetuptime += m_SetupTime[m][preFam][Fam];
            }
            preFam = Fam;
            for (int j = 0; j < JobSeqInFam[Fam].size(); j++) {
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


void
NOperator_New::GetJFDHierarchy_Forward_New(const vector<vector<int>> &FacFamSeq,
                                           const vector<vector<int>> &JobSeqInFam,
                                           vector<vector<int>> &JFDTime, vector<vector<int>> &Hierarchy)
{
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
        GetJFDHierarchy_Forward_InFactory_New(FacFamSeq[Fac], JobSeqInFam, JFDTime,Hierarchy);//�õ�ÿ���������깤ʱ��
}


/**
 * �¿���
 * ���㲿�����еĲ�κ͹����뿪ʱ��
 * @param FamSeqInFac
 * @param JobSeqInFam
 * @param JFDTime
 * @param hierarchy
 */
void
NOperator_New::GetJFDHierarchy_Forward_InFactory_New(const vector<int> &FamSeqInFac,
                                                     const vector<vector<int>> &JobSeqInFam,
                                                     vector<vector<int>> &JFDTime, vector<vector<int>> &hierarchy) {

    int PreFam; //ǰһ����
    for (int Fam = 0; Fam < FamSeqInFac.size(); ++Fam) {
        //��ĵ���
        int CurFam = FamSeqInFac[Fam]; //��ǰ��
        int FirstJobInFam = JobSeqInFam[Fam][0];
        if (Fam == 0) {
            for (int m = 0; m < this->m_Machines; m++)
                MachReadyTime[m] = this->m_SetupTime[m][CurFam][CurFam];
        } else {
            int LastJobInPreFam = *JobSeqInFam[PreFam].rbegin();
            for (int m = 0; m < this->m_Machines; m++)
                MachReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
        }
        //���ڹ����ĵ���
        for (int j = 0; j < JobSeqInFam[CurFam].size(); ++j) {
            int CurJob = JobSeqInFam[CurFam][j];//��ǰ����
            //����һ����������¼��κ��뿪ʱ��
            scheduling_one_job(CurJob, MachReadyTime, hierarchy[CurJob]);
            JFDTime[CurJob] = MachReadyTime;
        }
        PreFam = CurFam;
    }
}


/**
 * �¿���
 * ���Խ�����뵽�����ڵ�ĳ��λ��
 * @param FamSeqInFac
 * @param JobSeqInFam
 * @param JFDTime
 * @param InsFam
 * @param Pos
 * @param Hierarchy
 * @param scenario
 * @return
 */
int NOperator_New::GetSpanForPerFacAfterInsertFam_New(const vector<int> &FamSeqInFac,
                                                      const vector<vector<int>> &JobSeqInFam,
                                                      const vector<vector<int>> &JFDTime,
                                                      const vector<vector<int>> &Hierarchy,
                                                      int InsFam, int Pos) {

    int Span = 0;
    //����λ��֮ǰ����
    for (int f = 0; f < Pos; f++) {
        int Fam = FamSeqInFac[f];
        int LastJobInFam = *JobSeqInFam[Fam].rbegin();
        Span = JFDTime[LastJobInFam][this->m_Machines - 1];
    }

    //����Ҫ�������
    int PreFam;
    if (Pos == 0) {
        //����λ���ǵ�һ��
        for (int m = 0; m < this->m_Machines; ++m) {
            MachReadyTime[m] = this->m_SetupTime[m][InsFam][InsFam]; //MachReadyTime��Ϊ׼��ʱ��
        }
    } else {
        //���뵽����λ��
        PreFam = FamSeqInFac[Pos - 1];//ǰһ����
        int LastJobInPreFam = *JobSeqInFam[PreFam].rbegin();
        for (int m = 0; m < this->m_Machines; m++) {
            //ǰһ�������һ���������뿪ʱ��+׼��ʱ��
            MachReadyTime[m] = JFDTime[LastJobInPreFam][m] + this->m_SetupTime[m][PreFam][InsFam];
        }
    }

    //��������Ĺ�������
    for (int CurJob: JobSeqInFam[InsFam]) {
        scheduling_one_job(CurJob, MachReadyTime);
    }
    //�������
    Span = MachReadyTime[this->m_Machines - 1];

    //�������λ��֮�����
    bool is_hierarchy_same = false;
    int offset;
    PreFam = InsFam;
    for (int Fam = Pos; Fam < FamSeqInFac.size(); ++Fam) {
        int CurFam = FamSeqInFac[Fam];
        for (int m = 0; m < this->m_Machines; m++) {
            MachReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
        }

        for (int j = 0; j < JobSeqInFam[CurFam].size(); ++j) {
            int CurJob = JobSeqInFam[CurFam][j];
            //����һ���������жϲ���Ƿ��ԭ����ͬ
            is_hierarchy_same = scheduling_one_job(CurJob, MachReadyTime, Hierarchy[CurJob]);
            //���������������еĲ�κ��ڲ��������еĲ��һ��������ƫ����
            if (is_hierarchy_same) {
                offset = MachReadyTime[this->m_Machines - 1] - JFDTime[CurJob][this->m_Machines - 1];
                break;
            }
        }
        if (is_hierarchy_same) {
            //����ƫ�������ٴ���������
            for (int f1 = Fam; f1 < FamSeqInFac.size(); ++f1) {
                CurFam = FamSeqInFac[f1];
                int LastJobInFam = *JobSeqInFam[CurFam].rbegin();
                Span = (JFDTime[LastJobInFam][this->m_Machines - 1] + offset);
            }
            break;
        } else {
            Span = MachReadyTime[this->m_Machines - 1];
            PreFam = CurFam;
        }
    }

    return Span;
}


/**
 * �¿���
 * @param NewFamSeqInFac
 * @param JobSeqInFam
 * @param JFDTime
 * @param Hierarchy
 * @param InsFam
 * @param BestPos
 * @return
 */
int NOperator_New::FindBestPosToInsertFam_InFactory_MakeSpan_New(const vector<int> &NewFamSeqInFac,
                                                                 const vector<vector<int>> &JobSeqInFam,
                                                                 const vector<vector<int>> &JFDTime,
                                                                 const vector<vector<int>> &Hierarchy,
                                                                 int InsFam, int &BestPos) {
    int minSpan = INT_MAX;
    int Span = 0;
    for (int Pos = 0; Pos <= NewFamSeqInFac.size(); ++Pos) {
        auto re = GetSpanForPerFacAfterInsertFam_New(NewFamSeqInFac, JobSeqInFam, JFDTime, Hierarchy, InsFam, Pos);
        Span = re;
        if (Span < minSpan) {
            minSpan = Span;
            BestPos = Pos;
        }
    }

    return minSpan;
}

/**
 * �¿���
 * �����й������ҵ�������makespan��С��λ��
 * @param FacFamSeq
 * @param JobSeqInFam
 * @param JFDTime
 * @param InsFam
 * @param BestFac
 * @param BestPos
 * @return
 */
int NOperator_New::FindBestPosToInsertFamForAllFacs_MakeSpan_New(const vector<vector<int>> &FacFamSeq,
                                                                 const vector<vector<int>> &JobSeqInFam,
                                                                 const vector<vector<int>> &JFDTime,
                                                                 const vector<vector<int>> &Hierarchy,
                                                                 int InsFam, int &BestFac, int &BestPos) {
    int minSpan = INT_MAX;
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++) {
        int Pos;
        int Span = FindBestPosToInsertFam_InFactory_MakeSpan_New(FacFamSeq[Fac], JobSeqInFam, JFDTime, Hierarchy,InsFam,Pos);
        if (Span < minSpan) {
            minSpan = Span;
            BestFac = Fac;
            BestPos = Pos;
        }
    }
    return minSpan;
}


/**
 * �¿���
 * ���Խ��������뵽���ڵ�ĳ��λ��
 * @param FamSeqInFac
 * @param JobSeqInFam
 * @param JFDTime
 * @param InsFam
 * @param FamPos
 * @param InsJob
 * @param hierarchy
 * @return
 */
int NOperator_New::GetSpanForPerFacAfterInsertJob_New(const vector<int> &FamSeqInFac,const vector<vector<int>> &JobSeqInFam,
                                                      const vector<vector<int>> &JFDTime,  const vector<vector<int>> &hierarchy,
                                                      int InsFam, int FamPos,int InsJob,int JobPos) {

    int Span = 0;
    //����λ��֮ǰ����
    for (int f = 0; f < FamPos; ++f) {
        int Fam = FamSeqInFac[f];
        int LastJobInFam = *JobSeqInFam[Fam].rbegin();
        Span = JFDTime[LastJobInFam][this->m_Machines - 1];
    }

    //���빤��λ�õ���
    if (JobPos == 0) {
        if (FamPos == 0) {
            //��һ����ĵ�һ��λ��
            for (int m = 0; m < this->m_Machines; m++) {
                //MachReadyTime��Ϊ׼��ʱ��
                MachReadyTime[m] = this->m_SetupTime[m][InsFam][InsFam];
            }
        } else {
            //������ĵ�һ��λ��
            int PreFam = FamSeqInFac[FamPos - 1];//ǰһ����
            int LastJobInPreFam = *JobSeqInFam[PreFam].rbegin();//ǰһ��������һ������
            for (int m = 0; m < this->m_Machines; ++m) {
                //ǰһ�������һ���������뿪ʱ��+׼��ʱ��
                MachReadyTime[m] = JFDTime[LastJobInPreFam][m] + this->m_SetupTime[m][PreFam][InsFam];
            }
        }
    } else {
        //���ǵ�һ��λ��
        int PreJob = JobSeqInFam[InsFam][JobPos - 1];
        for (int m = 0; m < this->m_Machines; m++) {
            MachReadyTime[m] = JFDTime[PreJob][m];
        }
    }

    //���Ȳ���Ĺ���
    scheduling_one_job(InsJob, MachReadyTime);

    bool is_hierarchy_same = false;
    int offset;

    //���빤�������ڲ��빤��λ��֮��Ĺ���
    if (JobPos < JobSeqInFam[InsFam].size()) {
        for (int j = JobPos; j < JobSeqInFam[InsFam].size(); j++) {
            int CurJob = JobSeqInFam[InsFam][j];
            //����һ���������жϲ���Ƿ��ԭ����ͬ
            is_hierarchy_same = scheduling_one_job(CurJob, MachReadyTime, hierarchy[CurJob]);
            //���������������еĲ�κ��ڲ��������еĲ��һ��������ƫ����
            if (is_hierarchy_same) {
                offset = MachReadyTime[this->m_Machines - 1] - JFDTime[CurJob][this->m_Machines - 1];
                break;
            }
        }

        if (is_hierarchy_same) {
            //����ƫ�������ٴ���������ڲ��빤������Ĺ���
            int LastJobInFam = *JobSeqInFam[InsFam].rbegin();
            //JFDTimeOfFamOnLastMachine[InsFam] = JFDTime[LastJobInFam][this->m_Machines - 1] + offset;
            Span = (JFDTime[LastJobInFam][this->m_Machines - 1] + offset);
        } else {
            Span = MachReadyTime[this->m_Machines - 1];
            //JFDTimeOfFamOnLastMachine[InsFam] = MachReadyTime[this->m_Machines - 1];
        }
    } else {
        //����Ĺ����ڲ���������һ��λ��
        Span = MachReadyTime[this->m_Machines - 1];
        //JFDTimeOfFamOnLastMachine[InsFam] = MachReadyTime[this->m_Machines - 1];
    }

    //������һ����������֮�����Ͳ��ü���ÿ��������
    if (is_hierarchy_same) {
        for (int Fam = FamPos + 1; Fam < FamSeqInFac.size(); ++Fam) {
            int CurFam = FamSeqInFac[Fam];
            int LastJobInFam = *JobSeqInFam[CurFam].rbegin();
            //JFDTimeOfFamOnLastMachine[CurFam] = JFDTime[LastJobInFam][this->m_Machines - 1] + offset;
            Span = (JFDTime[LastJobInFam][this->m_Machines - 1] + offset);
        }
    } else {
        //Ŀǰ��β�һ�������������ܻ�һ��
        int PreFam = InsFam;
        for (int Fam = FamPos + 1; Fam < FamSeqInFac.size(); ++Fam) {
            int CurFam = FamSeqInFac[Fam];
            for (int m = 0; m < this->m_Machines; m++) {
                MachReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
            }

            //�������ڵĹ���
            for (int j = 0; j < JobSeqInFam[CurFam].size(); ++j) {
                int CurJob = JobSeqInFam[CurFam][j];
                //����һ���������жϲ���Ƿ��ԭ����ͬ
                is_hierarchy_same = scheduling_one_job(CurJob, MachReadyTime, hierarchy[CurJob]);
                //���������������еĲ�κ��ڲ��������еĲ��һ��������ƫ����
                if (is_hierarchy_same) {
                    offset = MachReadyTime[this->m_Machines - 1] - JFDTime[CurJob][this->m_Machines - 1];
                    break;
                }
            }
            if (is_hierarchy_same) {
                //����ƫ��������������
                for (int f1 = Fam; f1 < FamSeqInFac.size(); ++f1) {
                    CurFam = FamSeqInFac[f1];
                    int LastJobInFam = *JobSeqInFam[CurFam].rbegin();
                    //JFDTimeOfFamOnLastMachine[CurFam] = JFDTime[LastJobInFam][this->m_Machines - 1] + offset;
                    Span = (JFDTime[LastJobInFam][this->m_Machines - 1] + offset);
                }
                break;
            } else {
                Span = MachReadyTime[this->m_Machines - 1];
                //JFDTimeOfFamOnLastMachine[CurFam] = MachReadyTime[this->m_Machines - 1];
                PreFam = CurFam;
            }
        }
    }
    return Span;
}


void NOperator_New::Speed_mutation(vector<Individual> &CCEAPopulation, vector<Individual> &Temp,
                                   vector<Individual> &tureCCEAPopulation) {

    for (int PS = 0; PS < CCEAPopulation.size(); PS++) {
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
        Makespan = GetJFDTime_Forward(Temp[PS].m_FacFamSeqArray, Temp[PS].m_JobSeqInFamArray, m_TempJFDTime,FacSpan);
        TotalEC = GetTECForAllFacByJFD(Temp[PS].m_FacFamSeqArray, Temp[PS].m_JobSeqInFamArray, m_TempJFDTime);
        Temp[PS].MS = Makespan;
        Temp[PS].TEC = TotalEC;
        if (((Temp[PS].MS < CCEAPopulation[PS].MS) && (Temp[PS].TEC < CCEAPopulation[PS].TEC)) ||
            ((Temp[PS].MS < CCEAPopulation[PS].MS) && (Temp[PS].TEC == CCEAPopulation[PS].TEC)) ||
            ((Temp[PS].MS == CCEAPopulation[PS].MS) && (Temp[PS].TEC < CCEAPopulation[PS].TEC))) {
            CCEAPopulation[PS].m_SpeedVector = Temp[PS].m_SpeedVector;
            CCEAPopulation[PS].MS = Temp[PS].MS;
            CCEAPopulation[PS].TEC = Temp[PS].TEC;
            tureCCEAPopulation.push_back(Temp[PS]);  //����ο���
        }

    }

}

/************************************************************************************************************************/
/****************************  �Թ���InsJob�ڹ�����JobSeqInFam[Fam]��Ѱ�����λ��  **************************************/
float NOperator_New::FindBestPosToInsertJobForPerFac_Ind_New(int fac ,const vector<vector<int>> &FacFamSeq,const vector<int>& FamSeqInFac, const vector<vector<int>> JobSeqInFam,
                                                             const vector<vector<int>>& JFDTime, const vector<vector<int>>& Hierarchy, int Fam, int FamPos, int InsJob,
                                                             int& BestPos, int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                                             float idealpointTEC,int OrgMS, float OrgTEC,
                                                             vector<Individual> &CCEAPopulation)
{
    float minInd = INT_MAX;
    for (int Pos = 0; Pos <= JobSeqInFam[Fam].size(); Pos++)
    {
        auto re = this->GetIndForPerFacAfterInsertJob_New(fac, FacFamSeq, FamSeqInFac, JobSeqInFam, JFDTime, Hierarchy, Fam, FamPos, InsJob, Pos,
                                                          nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, OrgMS, OrgTEC, CCEAPopulation);

        if (re < minInd)
        {
            minInd = re;
            BestPos = Pos;
        }
    }
    return minInd;
}

//����ָ����������й����������λ��
float NOperator_New::FindBestPosToInsertFamForAllFac_Ind_New(const vector<vector<int>> &FacFamSeq,
                                                             const vector<vector<int>> &JobSeqInFam,
                                                             const vector<vector<int>> &JFDTime,
                                                             const vector<vector<int>> &Hierarchy,
                                                             int InsFam, int &BestFac, int &BestPos, int nadirpointMS,
                                                             float nadirpointTEC, int idealpointMS, float idealpointTEC,vector<int>& FacSpan, vector<float>& FacEC,
                                                             int OrgMS, float OrgTEC,
                                                             vector<Individual> &CCEAPopulation)
 {
    float minInd = INT_MAX;
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++) {
        int Pos = -1;
        float Ind = this->FindBestPosToInsertFamForPerFac_Ind_New(Fac, FacFamSeq, FacFamSeq[Fac], JobSeqInFam, JFDTime,
                                                                  Hierarchy, InsFam, Pos, nadirpointMS, nadirpointTEC,
                                                                  idealpointMS, idealpointTEC, FacSpan, FacEC, OrgMS, OrgTEC, CCEAPopulation);
        if (Ind < minInd) {
            minInd = Ind;
            BestFac = Fac;
            BestPos = Pos;
        }
    }
    return minInd;
}


//�� ����ָ������ڵ�ǰ���������λ��
float NOperator_New::FindBestPosToInsertFamForPerFac_Ind_New(int Fac, const vector<vector<int>> &FacFamSeq,
                                                             const vector<int> &NewFamSeqInFac,
                                                             const vector<vector<int>> &JobSeqInFam,
                                                             const vector<vector<int>> &JFDTime,
                                                             const vector<vector<int>> &Hierarchy, int InsFam,
                                                             int &BestPos, int nadirpointMS,
                                                             float nadirpointTEC, int idealpointMS, float idealpointTEC,vector<int>& FacSpan, vector<float>& FacEC,
                                                             int OrgMS, float OrgTEC,
                                                             vector<Individual> &CCEAPopulation) {
    float minFacInd = INT_MAX;

    for (int Pos = 0; Pos <= NewFamSeqInFac.size(); Pos++) {

        float FacInd = GetIndForPerFacAfterInsertFam_New(Fac, FacFamSeq, NewFamSeqInFac, JobSeqInFam, JFDTime,
                                                         Hierarchy, InsFam, Pos,
                                                         nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, FacSpan, FacEC,
                                                         OrgMS, OrgTEC, CCEAPopulation);

        if (FacInd < minFacInd) {
            minFacInd = FacInd;
            BestPos = Pos;
        }
    }
    return minFacInd;
}

//ĳ��λ�ò���������ָ��Ind
float NOperator_New::GetIndForPerFacAfterInsertFam_New(int Fac, const vector<vector<int>> &FacFamSeq,
                                                       const vector<int> &FamSeqInFac,
                                                       const vector<vector<int>> &JobSeqInFam,
                                                       const vector<vector<int>> &JFDTime,
                                                       const vector<vector<int>> &Hierarchy, int InsFam, int FamPos,
                                                       int nadirpointMS,float nadirpointTEC, int idealpointMS, float idealpointTEC,vector<int>& FacSpan, vector<float>& FacEC,
                                                       int OrgMS, float OrgTEC, vector<Individual> &CCEAPopulation) {
    //�������󹤳�Fac��span
    int Span = GetSpanForPerFacAfterInsertFam_New(FamSeqInFac,JobSeqInFam,JFDTime,Hierarchy,InsFam,FamPos);

    //����δ����֮ǰ��������span
    for (int f = 0; f < FacFamSeq.size(); f++) {
        FacSpan[f]= GetSpanPerFacByJFD(FacFamSeq[f],JobSeqInFam,JFDTime);
    }

    FacSpan[Fac]= Span ;
    int MakeSpan = *max_element(FacSpan.begin(), FacSpan.end());

    //�ж��Ƿ��滻��͵�������
    if (MakeSpan > nadirpointMS)
        nadirpointMS = MakeSpan;
    if (MakeSpan < idealpointMS)
        idealpointMS = MakeSpan;

    vector<vector<int>> TempFacFamSeq;
    TempFacFamSeq.clear();
    TempFacFamSeq.resize(FacFamSeq.size());
    TempFacFamSeq = FacFamSeq;
    TempFacFamSeq[Fac].insert(TempFacFamSeq[Fac].begin() + FamPos, InsFam);
    for (int fac = 0; fac < TempFacFamSeq.size(); fac++)
    {
        //�����ܺ�
        float TpTimeEC = 0, TSetupTimeEC = 0, TIdleTimeEC = 0;

        for (int m = 0; m < m_Machines; m++)
        {
            int preFam = -1;
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

    float TEC = accumulate(FacEC.begin(), FacEC.end(), 0.0);
    //�ж��Ƿ��滻��͵�������
    if (TEC > nadirpointTEC)
        nadirpointTEC = TEC;
    if (TEC < idealpointTEC)
        idealpointTEC = TEC;

    //�ж��Ƿ�Ȳ���ǰ�Ľ������Ľ������ο���
    bool flag = true;
    if (((MakeSpan < OrgMS) && (TEC < OrgTEC)) || ((MakeSpan < OrgMS) && (TEC == OrgTEC)) ||
        ((MakeSpan == OrgMS) && (TEC < OrgTEC)) || (MakeSpan < idealpointMS) || (TEC < idealpointTEC)) {
        for (int i = 0; i < CCEAPopulation.size(); i++) {
            if ((MakeSpan == CCEAPopulation[i].MS) && (TEC == CCEAPopulation[i].TEC)) {
                flag = false;
                break;
            }
        }
        if (flag) {
            Individual tempIndi;
            tempIndi.m_FacFamSeqArray =  TempFacFamSeq;
            tempIndi.m_JobSeqInFamArray = JobSeqInFam;
            tempIndi.MS = MakeSpan;
            tempIndi.TEC = TEC;
            tempIndi.m_SpeedVector = m_SpeedMatrix;  //�ٶȾ���
            CCEAPopulation.push_back(tempIndi);
        }
    }

    //��һ��
    float normalMS = -1;
    float normalTEC = -1;

    normalMS = (static_cast<float>(MakeSpan - idealpointMS) / static_cast<float>(nadirpointMS - idealpointMS));
    normalTEC = ((TEC - idealpointTEC) / (nadirpointTEC - idealpointTEC));

    //����ָ��Ind
    float convergence_ind = 0;
    convergence_ind += (normalMS - 1.0) * (normalMS - 1.0);
    convergence_ind += (normalTEC - 1.0) * (normalTEC - 1.0);
    convergence_ind = sqrt(convergence_ind);
    convergence_ind = 1 / (convergence_ind + 1);


    // ���� calc_distribution_ind ����������ֲ���ָ��
    Individual::calc_distribution_ind(CCEAPopulation);
    float Distribution_ind = 0;
    for (int i = 0; i < CCEAPopulation.size(); i++) {
        for (int j = i + 1; j < CCEAPopulation.size(); j++) {
            Distribution_ind += CCEAPopulation[i].distribution_ind[j];
        }
    }
    Distribution_ind = Distribution_ind / (CCEAPopulation.size() * (CCEAPopulation.size() - 1) / 2);

   // cout<<" Distribution_ind: "<< Distribution_ind<<"--convergence_ind:"<<convergence_ind<<endl;
    // �������ռ�Ȩָ��
    float combined_ind = 0.5f * convergence_ind + 0.5f * Distribution_ind;

    return combined_ind;
}
//�� ĳ��λ�ò��빤�������ָ��Ind
float NOperator_New::GetIndForPerFacAfterInsertJob_New(int fac, const vector<vector<int>> &FacFamSeq,
                                                       const vector<int> &FamSeqInFac,
                                                       const vector<vector<int>> &JobSeqInFam,
                                                       const vector<vector<int>> &JFDTime,
                                                       const vector<vector<int>> &Hierarchy, int InsFam, int FamPos,
                                                       int InsJob, int JobPos,
                                                       int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                                       float idealpointTEC, int OrgMS, float OrgTEC,
                                                       vector<Individual> &CCEAPopulation) {

    vector <int> FacSpan(FacFamSeq.size(), 0);
    int Span = GetSpanForPerFacAfterInsertJob_New(FamSeqInFac,JobSeqInFam,JFDTime,Hierarchy,InsFam,FamPos,InsJob,JobPos);

    //����δ����֮ǰ��������span
    for (int f = 0; f < FacFamSeq.size(); f++) {
        FacSpan[f] = GetSpanPerFacByJFD(FacFamSeq[f],JobSeqInFam,JFDTime);
    }
    FacSpan[fac] = Span;
    int MakeSpan = *max_element(FacSpan.begin(), FacSpan.end());
    //cout << "�����MakeSpan��" << MakeSpan << endl;

    //�ж��Ƿ��滻��͵�������
    if (MakeSpan > nadirpointMS)
        nadirpointMS = MakeSpan;
    if (MakeSpan < idealpointMS)
        idealpointMS = MakeSpan;

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


//�� ����ָ����������й����������λ��
float NOperator_New::FindBestPosToInsertFamForAllFac_Ind_DR_New( const vector<vector<int>> &FacFamSeq,
                                                                const vector<vector<int>> &JobSeqInFam,
                                                                const vector<vector<int>> &JFDTime,
                                                                const vector<vector<int>> &Hierarchy,
                                                                int InsFam, int &BestFac, int &BestPos,
                                                                int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                                                float idealpointTEC, vector<Individual> &CCEAPopulation) {
    float minInd = INT_MAX;
    for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
    {
        int Pos = -1;
        float Ind = this->FindBestPosToInsertFamForPerFac_Ind_DR_New(Fac, FacFamSeq, FacFamSeq[Fac], JobSeqInFam, JFDTime, Hierarchy,
                                                                     InsFam, Pos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, CCEAPopulation);
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
float NOperator_New::FindBestPosToInsertFamForPerFac_Ind_DR_New(int Fac,const vector<vector<int>> &FacFamSeq,
                                                                const vector<int> &NewFamSeqInFac,
                                                                const vector<vector<int>> &JobSeqInFam,
                                                                const vector<vector<int>> &JFDTime,
                                                                const vector<vector<int>> &Hierarchy, int InsFam,
                                                                int &BestPos, int nadirpointMS, float nadirpointTEC,
                                                                int idealpointMS, float idealpointTEC, vector<Individual> &CCEAPopulation) {
    float minFacInd = INT_MAX;
    for (int Pos = 0; Pos <= NewFamSeqInFac.size(); Pos++)
    {

        float FacInd = GetIndForPerFacAfterInsertFam_DR_New(Fac, FacFamSeq, NewFamSeqInFac, JobSeqInFam, JFDTime, Hierarchy,
                                                            InsFam, Pos, nadirpointMS, nadirpointTEC, idealpointMS, idealpointTEC, CCEAPopulation);

        if (FacInd < minFacInd)
        {
            minFacInd = FacInd;
            BestPos = Pos;
        }
    }
    return minFacInd;
}


//�� ĳ��λ�ò���������ָ��Ind
float NOperator_New::GetIndForPerFacAfterInsertFam_DR_New(int Fac,  const  vector<vector<int>> &FacFamSeq,
                                                          const vector<int> &FamSeqInFac,
                                                          const vector<vector<int>> &JobSeqInFam,
                                                          const vector<vector<int>> &JFDTime,
                                                          const vector<vector<int>> &Hierarchy, int InsFam, int Pos,
                                                          int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                                          float idealpointTEC ,vector<Individual> &CCEAPopulation) {

    //����δ����֮ǰ��������span
    vector<int> FacSpan(FacFamSeq.size(), 0);
    int Span = GetSpanForPerFacAfterInsertFam_New(FamSeqInFac,JobSeqInFam,JFDTime,Hierarchy,InsFam,Pos);

    for (int f = 0; f < FacFamSeq.size(); f++) {
        FacSpan[f]= GetSpanPerFacByJFD(FacFamSeq[f],JobSeqInFam,JFDTime);
      //  cout << "δ����ʱ����" << f << "��Span��" << FacSpan[f] << endl;
    }

    FacSpan[Fac] = Span;
    int MakeSpan = *max_element(FacSpan.begin(), FacSpan.end());

    //�ж��Ƿ��滻��͵�������
    if (MakeSpan > nadirpointMS)
        nadirpointMS = MakeSpan;
    if (MakeSpan < idealpointMS)
        idealpointMS = MakeSpan;

    vector<float> FacEC(FacFamSeq.size(), 0);
    float TEC;
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

            int Fam = -1;  // �ȳ�ʼ��
            if (!TempFacFamSeq[fac].empty()) {
                Fam = TempFacFamSeq[fac][0];  // ֻ���ڷǿ�����²Ÿ�ֵ
            } else {
                FacEC[fac] = 0;  // û����Ĺ�����ֱ����Ϊ 0
                continue;
            }
//            int Fam = TempFacFamSeq[fac][0];
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

//�� ĳ��λ�ò��빤�������ָ��Ind
float NOperator_New::GetIndForPerFacAfterInsertJob_DR_New(int fac, const vector<vector<int>> &FacFamSeq,
                                                          const vector<int> &FamSeqInFac,
                                                          const vector<vector<int>> &JobSeqInFam,
                                                          const vector<vector<int>> &JFDTime,
                                                          const vector<vector<int>> &Hierarchy, int InsFam, int FamPos,
                                                          int InsJob, int JobPos,
                                                          int nadirpointMS, float nadirpointTEC, int idealpointMS,
                                                          float idealpointTEC ,vector<Individual> &CCEAPopulation) {

    //����δ����֮ǰ��������span
    vector<int> FacSpan(FacFamSeq.size(), 0);
    int Span =   GetSpanForPerFacAfterInsertJob_New(FamSeqInFac,JobSeqInFam,JFDTime,Hierarchy,InsFam,FamPos,InsJob,JobPos);
    //  cout << "���빤��" << fac << "���Span��" << FacSpan[fac] << endl;

    for (int f = 0; f < FacFamSeq.size(); ++f) {
        FacSpan[f] = GetSpanPerFacByJFD(FacFamSeq[f],JobSeqInFam,JFDTime);
    }
    FacSpan[fac]= Span;
    int Makespan = *max_element(FacSpan.begin(), FacSpan.end());
    //cout << "�����MakeSpan��" << Makespan << endl;

    //�ж��Ƿ��滻��͵�������
    if (Makespan > nadirpointMS)
        nadirpointMS = Makespan;
    if (Makespan < idealpointMS)
        idealpointMS = Makespan;
    //�����ܺ�
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

void NOperator_New::CopyJFDHierarchy(const vector<vector<int>> &srcJFDTime, const vector<vector<int>> &srcHierarchy,
                                     vector<vector<int>> &dstJFDTime, vector<vector<int>> &dstHierarchy) {
    for (int j = 0; j < this->m_Jobs; ++j) {
        copy(begin(srcJFDTime[j]), end(srcJFDTime[j]), begin(dstJFDTime[j]));
        copy(begin(srcHierarchy[j]), end(srcHierarchy[j]), begin(dstHierarchy[j]));
    }
}