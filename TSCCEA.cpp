#include "TSCCEA.h"
using namespace std;

TSCCEA::TSCCEA()
{
}

TSCCEA::~TSCCEA()
{
}

void TSCCEA::SetParameters(int RefSize, int PS1, int PS2, int AgeLimit, int DesLen, long TimeLimit) //PS1>RefSize; PS2>RefSize;
{
    m_RefSize = RefSize;
    m_PS1 = PS1;
    m_PS2 = PS2;
    m_AgeLimit = AgeLimit;
    m_DesLen = DesLen;
    m_TimeLimit = TimeLimit;

    TSCCEAPopulation.clear();
    TSCCEAPopulation.resize(m_RefSize);

    m_nadirpointMS = -1;
    m_nadirpointTEC = -1;

    m_idealpointMS = -1;
    m_idealpointTEC = -1;

    //档案集 Reference set
    m_RefSpanArray.clear();//档案集完工时间
    m_RefSpanArray.resize(m_RefSize);

    m_RefTECArray.clear();//档案集总能耗
    m_RefTECArray.resize(m_RefSize);

    m_RefSpeedVector.clear();
    m_RefSpeedVector.resize(m_RefSize);


    m_RefFacFamSeqArray.clear();//档案集工厂组
    m_RefFacFamSeqArray.resize(m_RefSize);

    m_RefJobSeqinFamArray.clear();//档案集组工件
    m_RefJobSeqinFamArray.resize(m_RefSize);

    //m_ArchiveFacSpanArray.clear();//工厂完工时间
    //m_ArchiveFacSpanArray.resize(m_RefSize);

    m_bFlag1.clear();//
    m_bFlag1.resize(m_RefSize);

    m_bFlag2.clear();//
    m_bFlag2.resize(m_RefSize);

    //组
    m_TECArray1.clear(); //能耗
    m_TECArray1.resize(m_PS1);
    m_SpanArray1.clear();
    m_SpanArray1.resize(m_PS1);
    m_Map1.clear();
    m_Map1.resize(m_PS1);
    m_Age1.clear();
    m_Age1.resize(m_PS1);
    m_FacFamSeqArray.clear();//组序列
    m_FacFamSeqArray.resize(m_PS1);

    //工件
    m_TECArray2.clear(); //能耗
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

void TSCCEA::InitPop()
{
    //初始化档案集 Initialize Reference set and Set flags
    for (int PS = 0; PS < m_RefSize; PS++)
    {
        //变速 及单位加工能耗
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                m_SpeedMatrix[j][i] = rand() % 3;
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }


        vector<int> FamPrmu(m_Families), FacSpan(m_Factories);//组，工厂完工时间
        vector<float> FacEC(m_Factories);
        vector<vector<int>> JobSeqinFam, FacFamSeq, newFacFamSeq(m_Factories); //工件在组的序列，组在工厂的序列，新组在工厂的序列

        if (PS == 0) //第一个参考解用LPT，其它随机生成 Generate job sequence in each family and family sequence using LPT
        {

            GetJobTotalPTime();  //工件总加工时间
            GetFamTotalPTime();  //组总加工时间
            SortJobsinFam(0, 0, JobSeqinFam);  //根据工件总加工时间对组内工件排序
            SortFam(0, 0, FamPrmu);  //根据组加工时间对组排序
        }
        else // 随机生成工件序列和组序列 Generate job sequence in each family and family sequence randomly
        {
            JobSeqinFam = m_JobsInEachFamily;
            for (int fam = 0; fam < JobSeqinFam.size(); fam++)
                random_shuffle(JobSeqinFam[fam].begin(), JobSeqinFam[fam].end());//打乱组内工件顺序
            for (int fam = 0; fam < FamPrmu.size(); fam++)
                FamPrmu[fam] = fam;
            random_shuffle(FamPrmu.begin(), FamPrmu.end());//打乱组顺序
        }

        // 将上面生成的组分配到工厂产生一个解决方案 Generte a solution using assignment rule--include
        GetSol_Include(FamPrmu,JobSeqinFam, FacFamSeq, FacSpan);

        // 对每个工厂里的组进行NEH初始化
        if (PS == 0)
        {
            newFacFamSeq = FacFamSeq;
            FacFamSeq.clear();
            FacFamSeq.resize(m_Factories);
            for (int Fac = 0; Fac < m_Factories; Fac++)
                NEHFam(newFacFamSeq[Fac],JobSeqinFam, FacFamSeq[Fac], FacSpan[Fac]);
        }

        if (PS == 1 || PS == 2)
        {
            int MS0 = 0; float TEC0 = 0;
            //cout << endl << "SAS前种群中第" << PS << "个个体" << endl;
            GetMSandTECForPerandToalFac(FacFamSeq, JobSeqinFam, FacSpan, FacEC, MS0, TEC0);
           // cout << "Makespan：" << MS0 << "\t" << "TEC：" << TEC0 << endl;

            //检查
           // this->CheckSol(FacFamSeq, JobSeqinFam, MS0);
          //  this->CheckSolTEC(FacFamSeq, JobSeqinFam, TEC0);

            vector<vector<int>> JFDTime(this->m_Jobs);
            for (int j = 0; j < JFDTime.size(); j++)
                JFDTime[j].resize(this->m_Machines, 0);

            //得到JFD
            GetJFDTime_Forward(FacFamSeq, JobSeqinFam, JFDTime, FacSpan);

            //SAS--降低makespan，提高能耗
            for (int fac = 0; fac < m_Factories; fac++)
            {
                for (int fam = 0; fam < FacFamSeq[fac].size(); fam++)
                {
                    int CurFam = FacFamSeq[fac][fam];
                    for (int j = 0; j < JobSeqinFam[CurFam].size(); j++)
                    {
                        if ((fam != 0) && (j == 0))
                        {
                            int OrgFam = FacFamSeq[fac][fam - 1];
                            int last_job = JobSeqinFam[OrgFam][JobSeqinFam[OrgFam].size() - 1];
                            int CurJob = JobSeqinFam[CurFam][j];
                            for (int m = 0; m < m_Machines - 1; m++)
                            {
                                if (m == 0)
                                {
                                    if (JFDTime[last_job][m] + m_SetupTime[m][OrgFam][CurFam] + m_TureJobOpertime[CurJob][m] > JFDTime[last_job][m + 1] + m_SetupTime[m + 1][OrgFam][CurFam])
                                    {
                                        if (m_SpeedMatrix[CurJob][m] != 2)
                                        {
                                            int Speedlevel = m_SpeedMatrix[CurJob][m] + 1;
                                            int tempJobOpertime = static_cast<int>(100 * (m_JobOperPTime[CurJob][m] / m_Speed[Speedlevel]));
                                            if ((Speedlevel != 2) && (JFDTime[last_job][m] + m_SetupTime[m][OrgFam][CurFam] + tempJobOpertime > JFDTime[last_job][m + 1] + m_SetupTime[m + 1][OrgFam][CurFam]))
                                            {
                                                Speedlevel++;
                                                tempJobOpertime = static_cast<int>(100 * (m_JobOperPTime[CurJob][m] / m_Speed[Speedlevel]));
                                            }
                                            m_SpeedMatrix[CurJob][m] = Speedlevel;
                                            JFDTime[CurJob][m] = JFDTime[last_job][m] + tempJobOpertime;
                                        }
                                    }
                                }
                                else
                                {
                                    if (JFDTime[CurJob][m - 1] + m_TureJobOpertime[CurJob][m] > JFDTime[last_job][m + 1] + m_SetupTime[m + 1][OrgFam][CurFam])
                                    {
                                        if (m_SpeedMatrix[CurJob][m] != 2)
                                        {
                                            int Speedlevel = m_SpeedMatrix[CurJob][m] + 1;
                                            int tempJobOpertime = static_cast<int>(100 * (m_JobOperPTime[CurJob][m] / m_Speed[Speedlevel]));
                                            if ((Speedlevel != 2) && (JFDTime[CurJob][m - 1] + tempJobOpertime > JFDTime[last_job][m + 1] + m_SetupTime[m + 1][OrgFam][CurFam]))
                                            {
                                                Speedlevel++;
                                                tempJobOpertime = static_cast<int>(100 * (m_JobOperPTime[CurJob][m] / m_Speed[Speedlevel]));
                                            }
                                            m_SpeedMatrix[CurJob][m] = Speedlevel;
                                            //JFDTime[CurJob][m] = max(JFDTime[last_job][m], JFDTime[CurJob][m - 1]) + tempJobOpertime;
                                            JFDTime[CurJob][m] = max(JFDTime[last_job][m - 1] + tempJobOpertime,JFDTime[CurJob][m - 1]);
                                        }
                                    }
                                }
                            }
                        }
                        if (j > 0)
                        {
                            int last_job =JobSeqinFam[CurFam][j - 1];
                            int CurJob = JobSeqinFam[CurFam][j];
                            for (int m = 0; m < m_Machines - 1; m++)
                            {
                                if (m == 0)
                                {
                                    if (JFDTime[last_job][m] + m_TureJobOpertime[CurJob][m] > JFDTime[last_job][m + 1])
                                    {
                                        if (m_SpeedMatrix[CurJob][m] != 2)
                                        {
                                            int Speedlevel = m_SpeedMatrix[CurJob][m] + 1;
                                            int tempJobOpertime = static_cast<int>(100 * (m_JobOperPTime[CurJob][m] / m_Speed[Speedlevel]));
                                            if ((Speedlevel != 2) && (JFDTime[last_job][m] + tempJobOpertime > JFDTime[last_job][m + 1]))
                                            {
                                                Speedlevel++;
                                                tempJobOpertime = static_cast<int>(100 * (m_JobOperPTime[CurJob][m] / m_Speed[Speedlevel]));
                                            }
                                            m_SpeedMatrix[CurJob][m] = Speedlevel;
                                            JFDTime[CurJob][m] = JFDTime[last_job][m] + tempJobOpertime;
                                        }
                                    }
                                }
                                else
                                {
                                    if (JFDTime[CurJob][m - 1] + m_TureJobOpertime[CurJob][m] > JFDTime[last_job][m + 1])
                                    {
                                        if (m_SpeedMatrix[CurJob][m] != 2)
                                        {
                                            int Speedlevel = m_SpeedMatrix[CurJob][m] + 1;
                                            int tempJobOpertime = static_cast<int>(100 * (m_JobOperPTime[CurJob][m] / m_Speed[Speedlevel]));
                                            if ((Speedlevel != 2) && (JFDTime[CurJob][m - 1] + tempJobOpertime > JFDTime[last_job][m + 1]))
                                            {
                                                Speedlevel++;
                                                tempJobOpertime = static_cast<int>(100 * (m_JobOperPTime[CurJob][m] / m_Speed[Speedlevel]));
                                            }
                                            m_SpeedMatrix[CurJob][m] = Speedlevel;
                                           // JFDTime[CurJob][m] = max(JFDTime[last_job][m], JFDTime[CurJob][m - 1]) + tempJobOpertime;
                                            JFDTime[CurJob][m] = max(JFDTime[last_job][m - 1] + tempJobOpertime,JFDTime[CurJob][m - 1]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }

            //变速 及单位加工能耗
            for (int j = 0; j < m_Jobs; j++)
            {
                for (int i = 0; i < m_Machines; i++)
                {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }

        }

        int MS = 0; float TEC = 0;
        //cout << endl << "种群中第" << PS << "个个体" << endl;
        GetMSandTECForPerandToalFac(FacFamSeq, JobSeqinFam, FacSpan, FacEC, MS, TEC);
        //cout << "Makespan：" << MS << "\t" << "TEC：" << TEC << endl;


        // 初始化档案集（赋值） Initialize a reference
        m_RefFacFamSeqArray[PS] = FacFamSeq;   //组序列
        m_RefJobSeqinFamArray[PS] = JobSeqinFam;  //工件序列
        m_RefSpanArray[PS] = MS;  // 最大完工时间
        m_RefTECArray[PS] = TEC;  //总能耗
        m_RefSpeedVector[PS] = m_SpeedMatrix;
        m_bFlag1[PS] = false;
        m_bFlag2[PS] = false;

        //更新参考集
        TSCCEAPopulation[PS].m_FacFamSeqArray = FacFamSeq;
        TSCCEAPopulation[PS].m_JobSeqInFamArray = JobSeqinFam;
        TSCCEAPopulation[PS].MS = MS;
        TSCCEAPopulation[PS].TEC = TEC;
        TSCCEAPopulation[PS].m_SpeedVector = m_SpeedMatrix;  //速度矩阵

        //检查
//        this->CheckSol(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], m_ArchiveSpanArray[PS]);
//        this->CheckSolTEC(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], m_ArchiveTECArray[PS]);
//        CheckSol(TSCCEAPopulation[PS].m_FacFamSeqArray, TSCCEAPopulation[PS].m_JobSeqInFamArray, TSCCEAPopulation[PS].MS);

    }


    //组初始化 Initilize Familiy-sequence population, i.e., PS1

    //前五个种群 --> 用JPA生成的组序列，然后随机配对jobpop
    for (int PS = 0; PS < m_PS1 / 2; PS++)
    {
        //变速 及单位加工能耗
        m_SpeedMatrix = m_RefSpeedVector[PS % 5];
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }


        vector<vector<int>> newFacFamSeq(m_Factories); //每个工厂中组的序列

        //JPA
        JPA_TS(newFacFamSeq);

        m_FacFamSeqArray[PS] = newFacFamSeq;
        m_Map1[PS] = rand() % m_RefSize;  //随机选一个jobpop与当前fampop组成一对

        vector<int> FacSpan(m_Factories, 0);
        vector<float> FacEC(m_Factories, 0);

        GetMSandTECForPerandToalFac(m_FacFamSeqArray[PS], m_RefJobSeqinFamArray[m_Map1[PS]], FacSpan, FacEC, m_SpanArray1[PS], m_TECArray1[PS]);
        //检查
//        this->CheckSol(m_FacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[m_Map1[PS]], m_SpanArray1[PS]);
//        this->CheckSolTEC(m_FacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[m_Map1[PS]], m_TECArray1[PS]);

    }

    //后五个种群 --> 从参考集中随机选组序列，执行组交换(工厂内或者工厂间) 然后随机配对jobpop
    for (int PS = m_PS1 / 2; PS < m_PS1; PS++)
    {
        //变速 及单位加工能耗
        m_SpeedMatrix = m_RefSpeedVector[PS % 5];
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }


        //randomly select a Family Sequence from reference set and Swap
        int r = rand() % m_RefSize;
        m_FacFamSeqArray[PS] = m_RefFacFamSeqArray[r];

        SwapFam(m_FacFamSeqArray[PS]);
        //随机选一个jobpop从档案集中 randomly select a job seqence from reference set
        m_Map1[PS] = rand() % m_RefSize;  //随机选一个jobpop与当前fampop组成一对
        vector<int> FacSpan(m_Factories, 0);
        vector<float> FacEC(m_Factories, 0);

        GetMSandTECForPerandToalFac(m_FacFamSeqArray[PS], m_RefJobSeqinFamArray[m_Map1[PS]], FacSpan, FacEC, m_SpanArray1[PS], m_TECArray1[PS]);


        //检查
//        this->CheckSol(m_FacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[m_Map1[PS]], m_SpanArray1[PS]);
//        this->CheckSolTEC(m_FacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[m_Map1[PS]], m_TECArray1[PS]);
    }


    //Initialize Job-Sequence population, i.e., PS2
    //前五个种群 --> 从 m_JobsInEachFamil选择工件序列(即一开始读到的工件序列)，然后将工件找最好位置i插入， 然后随机配对fampop
    for (int PS = 0; PS < m_PS2 / 2; PS++)
    {
        //变速 及单位加工能耗
        m_SpeedMatrix = m_RefSpeedVector[PS % 5];
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }


        //randomly select a Job Sequence from reference set and Swap
        m_JobSeqinFamArray[PS] = m_JobsInEachFamily;
        m_Map2[PS] = rand() % m_RefSize;
        vector<int> FacSpan(m_Factories, 0);
        vector<float> FacEC(m_Factories, 0);

        for (int fac = 0; fac < m_Factories; fac++)
        {
            JobInsert(m_RefFacFamSeqArray[m_Map2[PS]][fac], m_JobSeqinFamArray[PS], FacSpan[fac]);
        }

        GetMSandTECForPerandToalFac(m_RefFacFamSeqArray[m_Map2[PS]], m_JobSeqinFamArray[PS], FacSpan, FacEC, m_SpanArray2[PS], m_TECArray2[PS]);

        //检查
//        this->CheckSol(m_ArchiveFacFamSeqArray[m_Map2[PS]], m_JobSeqinFamArray[PS], m_SpanArray2[PS]);
//        this->CheckSolTEC(m_ArchiveFacFamSeqArray[m_Map2[PS]], m_JobSeqinFamArray[PS], m_TECArray2[PS]);

    }


    //后五个种群 --> 从 参考及随机选择工件序列，然后将工件交换， 然后随机配对fampop
    for (int PS = m_PS2 / 2; PS < m_PS2; PS++)
    {
        //变速 及单位加工能耗
        m_SpeedMatrix = m_RefSpeedVector[PS % 5];
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }


        //randomly select a Job Sequence from reference set and Swap
        int r = rand() % m_RefSize;
        m_JobSeqinFamArray[PS] = m_RefJobSeqinFamArray[r];
        SwapJob(m_JobSeqinFamArray[PS]);
        m_Map2[PS] = rand() % m_RefSize;

        vector<int> FacSpan(m_Factories);
        vector<float> FacEC(m_Factories, 0);
        GetMSandTECForPerandToalFac(m_RefFacFamSeqArray[m_Map2[PS]], m_JobSeqinFamArray[PS], FacSpan, FacEC, m_SpanArray2[PS], m_TECArray2[PS]);

        //检查
        //this->CheckSol(m_ArchiveFacFamSeqArray[m_Map2[PS]], m_JobSeqinFamArray[PS], m_SpanArray2[PS]);
        //this->CheckSolTEC(m_ArchiveFacFamSeqArray[m_Map2[PS]], m_JobSeqinFamArray[PS], m_TECArray2[PS]);
    }

}


int TSCCEA::Destruction_Construction(int Len, vector<vector<int>>& FacFamSeq, vector<vector<int>> JobSeqInFam, vector<int>& FacSpan)
{
    // 找到关键工厂 Find critical factories
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
        if (ExtractFamSeq.size() < Len / 2) //从关键工厂挑选d/2个组
            Fac = CriFac[rand() % CriFac.size()];
        else  //从其他工厂挑选剩下的d/2个组
            Fac = rand() % m_Factories;
        if (FacFamSeq[Fac].size() > 1)
        {
            int Pos = rand() % FacFamSeq[Fac].size();
            ExtractFamSeq.push_back(FacFamSeq[Fac][Pos]);  //将挑选出的组添加到ExtractFamSeq
            FacFamSeq[Fac].erase(FacFamSeq[Fac].begin() + Pos);  //从原序列中删除pos位置的组
        }
    } while (ExtractFamSeq.size() < Len);

      GetSpan(FacFamSeq, JobSeqInFam, FacSpan);
    // 将挑选出来的d个组插入到最好的位置 Insert the extracted Families into the best Positions
    return NEHFam(ExtractFamSeq, JobSeqInFam, FacFamSeq, FacSpan);
}


void TSCCEA::ReInitPop() //reproduce a solution if age Limit is met
{

    //组重新初始化
    for (int PS = 0; PS < m_PS1; PS++)
    {
        m_SpeedMatrix = m_RefSpeedVector[PS % 5];
        //得到真正的处理时间 及单位加工能耗
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }

        if (m_Age1[PS] < m_AgeLimit) continue;
        m_FacFamSeqArray[PS] = m_RefFacFamSeqArray[m_Map1[PS]];

        vector<int> FacSpan(m_Factories, 0);
        vector<float> FacEC(m_Factories, 0);
        int MS = 0; float TEC = 0;
        //cout << endl << "组进化中第" << PS << "个个体" << endl;
        GetMSandTECForPerandToalFac(m_FacFamSeqArray[PS], m_RefJobSeqinFamArray[m_Map1[PS]], FacSpan, FacEC, MS, TEC);
        //cout << "组重新初始化前的MS：" << MS << "\tTEC：" << TEC << endl;

        m_SpanArray1[PS] = Destruction_Construction(m_DesLen, m_FacFamSeqArray[PS], m_RefJobSeqinFamArray[m_Map1[PS]], FacSpan);//解构——重构
        m_Age1[PS] = 0;

        GetMSandTECForPerandToalFac(m_FacFamSeqArray[PS], m_RefJobSeqinFamArray[m_Map1[PS]], FacSpan, FacEC, MS, TEC);
        //cout << "组重新初始化后的MS：" << MS << "\tTEC：" << TEC << endl;

        //检查
        //CheckSol(m_FacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[m_Map1[PS]], MS);
        //CheckSolTEC(m_FacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[m_Map1[PS]], TEC);

        //加入参考集
        bool flag = true;
        if (((TEC < m_TECArray1[PS]) && (MS < m_SpanArray1[PS])) || ((TEC < m_TECArray1[PS]) && (MS == m_SpanArray1[PS])) || ((TEC == m_TECArray1[PS]) && (MS < m_SpanArray1[PS])))
        {
            for (int i = 0; i < TSCCEAPopulation.size(); i++)
            {
                if ((MS == TSCCEAPopulation[i].MS) && (TEC == TSCCEAPopulation[i].TEC))
                {
                    flag = false;
                    break;
                }
            }
            if (flag)
            {
                //cout << "*********************" << endl;
                Individual tempIndi;
                tempIndi.m_FacFamSeqArray = m_FacFamSeqArray[PS];
                tempIndi.m_JobSeqInFamArray = m_RefJobSeqinFamArray[m_Map1[PS]];
                tempIndi.MS = MS;
                tempIndi.TEC = TEC;
                tempIndi.m_SpeedVector = m_SpeedMatrix;  //速度矩阵
                TSCCEAPopulation.push_back(tempIndi);
            }
        }

    }


    //工件重新初始化
    for (int PS = 0; PS < m_PS2; PS++)
    {
        m_SpeedMatrix = m_RefSpeedVector[PS % 5];
        //得到真正的处理时间 及单位加工能耗
        for (int j = 0; j < m_Jobs; j++)
        {
            for (int i = 0; i < m_Machines; i++)
            {
                m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
            }
        }


        if (m_Age2[PS] < m_AgeLimit) continue;
        m_JobSeqinFamArray[PS] = m_RefJobSeqinFamArray[m_Map2[PS]];

        vector<int> FacSpan(m_Factories, 0);
        vector<float> FacEC(m_Factories, 0);
        int MS = 0; float TEC = 0;
        //cout << endl << "工件进化中第" << PS << "个个体" << endl;
        GetMSandTECForPerandToalFac(m_RefFacFamSeqArray[m_Map2[PS]], m_JobSeqinFamArray[PS], FacSpan, FacEC, MS, TEC);
        //cout << "工件重新初始化前的MS：" << MS << "\tTEC：" << TEC << endl;

        SwapJob(m_JobSeqinFamArray[PS]);
        SwapJob(m_JobSeqinFamArray[PS]);
        m_Age2[PS] = 0;

        GetMSandTECForPerandToalFac(m_RefFacFamSeqArray[m_Map2[PS]], m_JobSeqinFamArray[PS], FacSpan, FacEC, MS, TEC);
        //cout << "工件重新初始化后的MS：" << MS << "\tTEC：" << TEC << endl;

        //检查
        //CheckSol(m_ArchiveFacFamSeqArray[m_Map2[PS]], m_JobSeqinFamArray[PS], MS);
        //CheckSolTEC(m_ArchiveFacFamSeqArray[m_Map2[PS]], m_JobSeqinFamArray[PS], TEC);

        //加入参考集
        int flag = true;
        if (((TEC < m_TECArray2[PS]) && (MS < m_SpanArray2[PS])) || ((TEC < m_TECArray2[PS]) && (MS == m_SpanArray2[PS])) || ((TEC == m_TECArray2[PS]) && (MS < m_SpanArray2[PS])))
        {
            for (int i = 0; i < TSCCEAPopulation.size(); i++)
            {
                if ((MS == TSCCEAPopulation[i].MS) && (TEC == TSCCEAPopulation[i].TEC))
                {
                    flag = false;
                    break;
                }
            }
            if (flag)
            {
                //cout << "*********************" << endl;
                Individual tempIndi;
                tempIndi.m_FacFamSeqArray = m_RefFacFamSeqArray[m_Map2[PS]];
                tempIndi.m_JobSeqInFamArray = m_JobSeqinFamArray[PS];
                tempIndi.MS = MS;
                tempIndi.TEC = TEC;
                tempIndi.m_SpeedVector = m_SpeedMatrix;  //速度矩阵
                TSCCEAPopulation.push_back(tempIndi);
            }
        }
    }
}

void TSCCEA::Evolution()//randomly setup the mapping relationship, compare with the original relationship, considering span and the number of critical pathes
{
    //cout << "HELLO" << endl;
    int count = 0;
    m_InitTime = Base::GetElapsedProcessTime();
    while (true)
    {
        count++;
        //------------Evolution for Job Populaiton-----------
        for (int PS = 0; PS < m_PS2; PS++)
        {
            m_SpeedMatrix = m_RefSpeedVector[PS % 5];
            //得到真正的处理时间 及单位加工能耗
            for (int j = 0; j < m_Jobs; j++)
            {
                for (int i = 0; i < m_Machines; i++)
                {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }

            vector<vector<int>> JobSeqinFam = m_JobSeqinFamArray[PS];
            int Map2 = rand() % m_RefSize;//form a solution by randomly selecting a ReFacFamseq

            vector<int> FacSpan(m_Factories, 0);
            vector<float> FacEC(m_Factories, 0);
            int MS = 0; float TEC = 0;
            //cout << endl << "工件进化中第" << PS << "个个体" << endl;
            GetMSandTECForPerandToalFac(m_RefFacFamSeqArray[Map2], JobSeqinFam, FacSpan, FacEC, MS, TEC);
            //cout << "工件进化前的MS：" << MS << "\tTEC：" << TEC << endl;


            IG_DR(m_RefFacFamSeqArray[Map2], JobSeqinFam, FacSpan);  // IG Destruction(不对) and Construction
            //cout << "1.55555555555" << endl;

            for (int fac = 0; fac < m_Factories; fac++)
            {
                JobInsert(m_RefFacFamSeqArray[Map2][fac], JobSeqinFam, FacSpan[fac]); //FLS
                //JobSwap(m_ArchiveFacFamSeqArray[Map2][fac], JobSeqinFam, FacSpan[fac]);
            }

            GetMSandTECForPerandToalFac(m_RefFacFamSeqArray[Map2], JobSeqinFam, FacSpan, FacEC, MS, TEC);
            //cout << "工件进化后的MS：" << MS << "\tTEC：" << TEC << endl;

            //检查
            //CheckSol(m_ArchiveFacFamSeqArray[Map2], JobSeqinFam, MS);
            //CheckSolTEC(m_ArchiveFacFamSeqArray[Map2], JobSeqinFam, TEC);

            //1: The original Maping Reference FamSeq is not changed; 2. Changed. reset corresponding relations
            if ((!m_bFlag1[m_Map2[PS]] && MS < m_SpanArray2[PS]) || m_bFlag1[m_Map1[PS]])//Update solution
            {
                m_SpanArray2[PS] = MS;
                m_Map2[PS] = Map2;
                m_JobSeqinFamArray[PS] = JobSeqinFam;
                m_Age2[PS] = 0;
            }

            //Update RefSet;
            bool flag = true;
            if (MS < m_RefSpanArray[Map2])
            {
                m_RefSpanArray[Map2] = MS;
                m_RefJobSeqinFamArray[Map2] = JobSeqinFam;
                m_bFlag2[Map2] = true;
                //加入参考集
                if ((TEC < m_RefTECArray[Map2]) || (TEC == m_RefTECArray[Map2]))
                {
                    for (int i = 0; i < TSCCEAPopulation.size(); i++)
                    {
                        if ((MS == TSCCEAPopulation[i].MS) && (TEC == TSCCEAPopulation[i].TEC))
                        {
                            flag = false;
                            break;
                        }
                    }
                    if (flag)
                    {
                        //cout << "*********************" << endl;
                        Individual tempIndi;
                        tempIndi.m_FacFamSeqArray = m_RefFacFamSeqArray[Map2];
                        tempIndi.m_JobSeqInFamArray = JobSeqinFam;
                        tempIndi.MS = MS;
                        tempIndi.TEC = TEC;
                        tempIndi.m_SpeedVector = m_SpeedMatrix;  //速度矩阵
                        TSCCEAPopulation.push_back(tempIndi);
                    }
                }
            }
            m_Age2[PS]++;
        }//end Job Evolution

        //Reset m_bFlag1;
        for (int PS = 0; PS < m_RefSize; PS++)
            m_bFlag1[PS] = false;

        //----------Evolution for Family Populaiton------------
        for (int PS = 0; PS < m_PS1; PS++)
        {
            m_SpeedMatrix = m_RefSpeedVector[PS % 5];
            //得到真正的处理时间 及单位加工能耗
            for (int j = 0; j < m_Jobs; j++)
            {
                for (int i = 0; i < m_Machines; i++)
                {
                    m_TureJobOpertime[j][i] = static_cast<int>(m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]);
                    UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                }
            }

            vector<vector<int>> FacFamSeq = m_FacFamSeqArray[PS];
            //vector<int> FamSeq = m_FamSeqArray[PS];
            int Map1 = rand() % m_RefSize;//form a solution by randomly selecting a RefJobSeq;
            vector<int> FacSpan(m_Factories, 0);
            vector<float> FacEC(m_Factories, 0);
            int MS = 0; float TEC = 0;
            //cout << endl << "组进化中第" << PS << "个个体" << endl;
            GetMSandTECForPerandToalFac(FacFamSeq, m_RefJobSeqinFamArray[Map1], FacSpan, FacEC, MS, TEC);
            //cout << "组进化前的MS：" << MS << "\tTEC：" << TEC << endl;

            int k = 1;
            while (k < 5)
            {
                // 生成随机数，减少 rand() 调用
                int randFactory = rand() % this->m_Factories;
                int randChoice = rand() % 2;  // 随机决定插入或交换

                // 根据 k 的值选择不同的操作
                switch (k)
                {
                    case 1:
                        InsertFamBetweenFac(FacFamSeq, m_RefJobSeqinFamArray[Map1]);
                        break;
                    case 2:
                        SwapFaminFac(FacFamSeq[randFactory]);
                        break;
                    case 3:
                        SwapFamBetweenFacs(FacFamSeq);
                        break;
                    case 4:
                        FamInsert(FacFamSeq[randFactory], m_RefJobSeqinFamArray[Map1], FacSpan[randFactory]);
                        break;
                }

                // 根据 randChoice 做进一步的插入或交换
                if (randChoice)
                    SwapFamBetweenFacs(FacFamSeq);
                else
                    InsertFamBetweenFac(FacFamSeq, m_RefJobSeqinFamArray[Map1]);

                k++;
            }


            GetMSandTECForPerandToalFac(FacFamSeq, m_RefJobSeqinFamArray[Map1], FacSpan, FacEC, MS, TEC);
            //cout << "组进化后的MS：" << MS << "\tTEC：" << TEC << endl;

            //检查
            //CheckSol(FacFamSeq, m_ArchiveJobSeqinFamArray[Map1], MS);
            //CheckSolTEC(FacFamSeq, m_ArchiveJobSeqinFamArray[Map1], TEC);

            //1: The original Maping Reference JobSeq is not changed; 2. Changed. reset corresponding relations
            if ((!m_bFlag2[m_Map1[PS]] && (MS < m_SpanArray1[PS])) || m_bFlag2[m_Map1[PS]])//Update solution
            {
                m_SpanArray1[PS] = MS;
                m_Map1[PS] = Map1;
                m_FacFamSeqArray[PS] = FacFamSeq;
                m_Age1[PS] = 0;

            }

            //Update RefSet;
            bool flag = true;
            if (MS < m_RefSpanArray[Map1])
            {
                m_RefSpanArray[Map1] = MS;
                m_RefFacFamSeqArray[Map1] = FacFamSeq;
                m_bFlag1[Map1] = true;
                //加入参考集
                if ((TEC < m_RefTECArray[Map1]) || (TEC == m_RefTECArray[Map1]))
                {
                    for (int i = 0; i < TSCCEAPopulation.size(); i++)
                    {
                        if ((MS == TSCCEAPopulation[i].MS) && (TEC == TSCCEAPopulation[i].TEC))
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
                        tempIndi.m_JobSeqInFamArray = m_RefJobSeqinFamArray[Map1];
                        tempIndi.MS = MS;
                        tempIndi.TEC = TEC;
                        tempIndi.m_SpeedVector = m_SpeedMatrix;  //速度矩阵
                        TSCCEAPopulation.push_back(tempIndi);
                    }
                }
            }
            m_Age1[PS]++;
        }//end Fam Evolution

        //Reset m_bFlag2;
        for (int PS = 0; PS < m_RefSize; PS++)
            m_bFlag2[PS] = false;

        //重新初始化
        ReInitPop();

        vector<Individual> Archive;
        Archive.clear();
        Archive.resize(m_RefSize);
        for (int PS = 0; PS < Archive.size(); PS++)
        {
            Archive[PS].m_JobSeqInFamArray = m_RefJobSeqinFamArray[PS];  //工件序列
            Archive[PS].m_FacFamSeqArray = m_RefFacFamSeqArray[PS];   //组序列
            Archive[PS].MS = m_RefSpanArray[PS];  // 最大完工时间
            Archive[PS].TEC = m_RefTECArray[PS];  //总能耗
            Archive[PS].m_SpeedVector = m_RefSpeedVector[PS];
        }

        //变速
        vector<Individual> Temp1;
        Temp1.clear();
        for (int i = 0; i < Archive.size(); i++)
        {
            Temp1.push_back(Archive[i]);
        }
        Speed_mutation(Archive, Temp1);

        //cout << "555555555" << endl;

        for (int PS = 0; PS < Archive.size(); PS++)
        {
            //m_ArchiveJobSeqinFamArray[PS] = Archive[PS].m_JobSeqInFamArray;  //工件序列
            //m_ArchiveFacFamSeqArray[PS] = Archive[PS].m_FacFamSeqArray;   //组序列
            m_RefSpanArray[PS] = Archive[PS].MS;  // 最大完工时间
            m_RefTECArray[PS] = Archive[PS].TEC;  //总能耗
            m_RefSpeedVector[PS] = Archive[PS].m_SpeedVector;
        }

        long CurTime = Base::GetElapsedProcessTime();
        long ElapsedTime = CurTime - m_InitTime;
        if (ElapsedTime >= m_TimeLimit) break;
    }

}
void TSCCEA::RunEvolution(int CPUFactor, int Rep, vector<vector<Individual>>& TSCCEAFinalAfterRepParetoSet)
{
    ReadInstanceFileNameList("..\\Benchmark\\");
    int RefSize = 5; //档案集
    int PS1 = 10;
    int PS2 = 10;
    int AgeLimit = 5;
    int DesLen = 4;
    int Instances = 405;
    vector<int> result;

    string FileDirectory = "..\\Result\\";
    ostringstream str;
    str << "TSCCEA_" << CPUFactor << "_experiment" << ".txt"; //不同的算法
    ofstream ofile;
    ofile.open(FileDirectory + str.str());

    for (int ins = 0; ins < Instances; ins++)
    {
        ReadInstance(ins);

        vector<Individual> FinalTSCCEAParetoSet;
        vector<Individual> AfterRepParetoSet;

        vector<Individual> TempAfterRepParetoSet;
        TempAfterRepParetoSet.clear();



        for (int r = 0; r < Rep; r++) {

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
            MachReadyTime.resize(this->m_Machines);

            long StartTime_Ins = Base::GetElapsedProcessTime();

            // 参数设置
            SetParameters(RefSize, PS1, PS2, AgeLimit, DesLen, CPUFactor * m_Families * m_Machines);
            //初始化种群
            InitPop();
            //进化
            Evolution();

            //非支配解
            TSCCEAParetoSet.clear();

            Pareto_relation(TSCCEAPopulation);

            //弱化支配解 de-emphasize dominated solutions
            for (int j = 0; j < TSCCEAPopulation.size(); j++)
            {
                for (int i = 0; i < TSCCEAPopulation.size(); i++)
                {
                    if (TSCCEAPopulation[i].flag == 0)
                    {
                        if (TSCCEAPopulation[i].pareto_rel[j] == 1)
                        {
                            TSCCEAPopulation[i].flag = 999;
                        }
                    }
                }
            }

            for (int i = 0; i < TSCCEAPopulation.size(); i++)
            {
                if (TSCCEAPopulation[i].flag == 0)
                    TSCCEAParetoSet.push_back(TSCCEAPopulation[i]);
            }

            //变速
            vector<Individual> Temp1;
            Temp1.clear();
            for (int i = 0; i < TSCCEAParetoSet.size(); i++)
            {
                Temp1.push_back(TSCCEAParetoSet[i]);
            }
            Speed_mutation(TSCCEAParetoSet, Temp1);

            //去除重复
            FinalTSCCEAParetoSet.clear();
            for (int i = 0; i < TSCCEAParetoSet.size(); i++)
            {
                bool fg = true;
                for (int j = 0; j < FinalTSCCEAParetoSet.size(); j++)
                {
                    if ((FinalTSCCEAParetoSet[j].MS == TSCCEAParetoSet[i].MS) && (FinalTSCCEAParetoSet[j].TEC == TSCCEAParetoSet[i].TEC))
                    {
                        fg = false;
                        break;
                    }
                }
                if (fg)
                {
                    FinalTSCCEAParetoSet.push_back(TSCCEAParetoSet[i]);
                }
            }

            long EndTime_Ins = Base::GetElapsedProcessTime();

            //检查
            for (int PS = 0; PS < FinalTSCCEAParetoSet.size(); PS++)
            {
                m_SpeedMatrix = FinalTSCCEAParetoSet[PS].m_SpeedVector;
                //得到真正的处理时间 及单位加工能耗
                for (int j = 0; j < m_Jobs; j++)
                {
                    for (int i = 0; i < m_Machines; i++)
                    {
                        m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
                        UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                    }
                }

                CheckSol(FinalTSCCEAParetoSet[PS].m_FacFamSeqArray, FinalTSCCEAParetoSet[PS].m_JobSeqInFamArray, FinalTSCCEAParetoSet[PS].MS);
                CheckSolTEC(FinalTSCCEAParetoSet[PS].m_FacFamSeqArray, FinalTSCCEAParetoSet[PS].m_JobSeqInFamArray, FinalTSCCEAParetoSet[PS].TEC);
            }


//            string Str = this->m_InstanceNameList[ins];
//
//            cout << Str << "\t" << this->m_Factories << "\t" << this->m_Machines << "\t" << this->m_Families << "\t"
//                 << this->m_SetupType << "\t" << this->m_Jobs << "\t";

            for (int PS = 0; PS < FinalTSCCEAParetoSet.size(); PS++)
            {
                //cout << FinalTSCCEAParetoSet[PS].MS << "\t" << FinalTSCCEAParetoSet[PS].TEC << "\t";
                TempAfterRepParetoSet.push_back(FinalTSCCEAParetoSet[PS]);
            }
        }

        // 非支配解
        AfterRepParetoSet.clear();

        Pareto_relation(TempAfterRepParetoSet);

        // 弱化支配解 de-emphasize dominated solutions
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

        // 去除重复，确保每个实例只剩下一个解
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

        cout << ins + 1 << "\t" << "Factories:" << this->m_Factories<< "\t" << "Machines :" << this->m_Machines<< "\t" << "Families:" << this->m_Families<< "\t" << "Jobs :" << this->m_Jobs << "\t";

        for (int PS = 0; PS < FinalAfterRepParetoSet.size(); PS++)
        {
            // 输出 MS 和 TEC，所有解都放在同一行
            cout << "MS:" << FinalAfterRepParetoSet[PS].MS
                 << " " << "TEC:" << FinalAfterRepParetoSet[PS].TEC << "\t";
            ofile << FinalAfterRepParetoSet[PS].MS << "  "<< FinalAfterRepParetoSet[PS].TEC << ","<< "\t";
            TSCCEAFinalAfterRepParetoSet[ins].push_back(FinalAfterRepParetoSet[PS]);
        }
        ofile << endl;  // 每个实例结果结束后换行
        cout << endl;   // 控制台输出换行，仅仅是在实例结束时换行
    }

    ofile.close();
}
