#include "GCCEA.h"
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

GCCEA::GCCEA()
{
}

GCCEA::~GCCEA()
{
}

void GCCEA::SetParameters(int RefSize, int PS1, int PS2, int AgeLimit, int DesLen, long TimeLimit) //PS1>RefSize; PS2>RefSize;
{
    m_RefSize = RefSize;
    m_PS1 = PS1;
    m_PS2 = PS2;
    m_AgeLimit = AgeLimit;
    m_DesLen = DesLen;
    m_TimeLimit = TimeLimit;

    GCCEAPopulation.clear();
    GCCEAPopulation.resize(m_RefSize);

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

void GCCEA::InitPop()
{
    //Initialize Reference set and Set flags 档案集和参考集
    for (int PS = 0; PS < m_RefSize; PS++)
    {
        //变速
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

        if (PS == 0)
        {
            //JPA
            JPA_G(newFacFamSeq, JobSeqinFam);
        }

        else
        {
            if (PS == 1) //第一个参考解用LPT，其它随机生成 Generate job sequence in each family and family sequence using LPT
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
            GetSol_Include(FamPrmu, JobSeqinFam, FacFamSeq, FacSpan);

            // 对每个工厂里的组进行NEH初始化
            for (int Fac = 0; Fac < m_Factories; Fac++)
                NEHFam(FacFamSeq[Fac], JobSeqinFam, newFacFamSeq[Fac], FacSpan[Fac]);
        }

        int MS = 0; float TEC = 0;
        //cout << endl << "种群中第" << PS << "个个体" << endl;
        GetMSandTECForPerandToalFac(newFacFamSeq, JobSeqinFam, FacSpan, FacEC, MS, TEC);
        //cout << "Makespan：" << MS << "\t" << "TEC：" << TEC << endl;


        // 初始化档案集（赋值） Initialize a reference
        m_RefJobSeqinFamArray[PS] = JobSeqinFam;  //工件序列
        m_RefFacFamSeqArray[PS] = newFacFamSeq;   //组序列
        //m_ArchiveFacSpanArray[PS] = FacSpan;  //所有工厂的完工时间
        m_RefSpanArray[PS] = MS;  // 最大完工时间
        m_RefTECArray[PS] = TEC;  //总能耗
        m_RefSpeedVector[PS] = m_SpeedMatrix;
        m_bFlag1[PS] = false;
        m_bFlag2[PS] = false;

        //更新参考集
        GCCEAPopulation[PS].m_FacFamSeqArray = newFacFamSeq;
        GCCEAPopulation[PS].m_JobSeqInFamArray = JobSeqinFam;
        GCCEAPopulation[PS].MS = MS;
        GCCEAPopulation[PS].TEC = TEC;
        GCCEAPopulation[PS].m_SpeedVector = m_SpeedMatrix;  //速度矩阵

        //检查
        //this->CheckSol(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], m_ArchiveSpanArray[PS]);
        //this->CheckSolTEC(m_ArchiveFacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[PS], m_ArchiveTECArray[PS]);
        //CheckSol(GCCEAPopulation[PS].m_FacFamSeqArray, GCCEAPopulation[PS].m_JobSeqInFamArray, GCCEAPopulation[PS].MS);

    }

    //组初始化 Initilize Familiy-sequence population, i.e., PS1
    for (int PS = 0; PS < min(m_PS1, m_RefSize); PS++)
    {
        m_FacFamSeqArray[PS] = m_RefFacFamSeqArray[PS];
        m_SpanArray1[PS] = m_RefSpanArray[PS];
        m_TECArray1[PS] = m_RefTECArray[PS];
        m_Map1[PS] = PS;  //标记档案集中对应的jobpop序号
    }

    for (int PS = min(m_PS1, m_RefSize); PS < m_PS1; PS++)
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
        //ShiftFam(m_FacFamSeqArray[PS]);
        SwapFam(m_FacFamSeqArray[PS]);
       // SwapFam(m_FacFamSeqArray[PS]);

        //随机选一个jobpop从档案集中 randomly select a job seqence from reference set
        m_Map1[PS] = rand() % m_RefSize;  //随机选一个jobpop与当前fampop组成一对
        vector<int> FacSpan(m_Factories, 0);
        vector<float> FacEC(m_Factories, 0);
        GetMSandTECForPerandToalFac(m_FacFamSeqArray[PS], m_RefJobSeqinFamArray[m_Map1[PS]], FacSpan, FacEC, m_SpanArray1[PS], m_TECArray1[PS]);

        //检查
        //this->CheckSol(m_FacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[m_Map1[PS]], m_SpanArray1[PS]);
        //this->CheckSolTEC(m_FacFamSeqArray[PS], m_ArchiveJobSeqinFamArray[m_Map1[PS]], m_TECArray1[PS]);
    }

    //Initialize Job-Sequence population, i.e., PS2
    for (int PS = 0; PS < min(m_PS2, m_RefSize); PS++)
    {
        m_JobSeqinFamArray[PS] = m_RefJobSeqinFamArray[PS];
        m_SpanArray2[PS] = m_RefSpanArray[PS];
        m_TECArray2[PS] = m_RefTECArray[PS];
        m_Map2[PS] = PS;
    }
    for (int PS = min(m_PS2, m_RefSize); PS < m_PS2; PS++)
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
        SwapJob(m_JobSeqinFamArray[PS]);
        //randomly select a Family seqence from reference set
        m_Map2[PS] = rand() % m_RefSize;

        vector<int> FacSpan(m_Factories);
        vector<float> FacEC(m_Factories, 0);
        GetMSandTECForPerandToalFac(m_RefFacFamSeqArray[m_Map2[PS]], m_JobSeqinFamArray[PS], FacSpan, FacEC, m_SpanArray2[PS], m_TECArray2[PS]);

        //检查
        //this->CheckSol(m_ArchiveFacFamSeqArray[m_Map2[PS]], m_JobSeqinFamArray[PS], m_SpanArray2[PS]);
        //this->CheckSolTEC(m_ArchiveFacFamSeqArray[m_Map2[PS]], m_JobSeqinFamArray[PS], m_TECArray2[PS]);
    }

}

int GCCEA::Destruction_Construction(int Len, vector<vector<int>>& FacFamSeq, vector<vector<int>> JobSeqinFam, vector<int>& FacSpan)
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
        if (FacFamSeq[Fac].size())
        {
            int Pos = rand() % FacFamSeq[Fac].size();
            ExtractFamSeq.push_back(FacFamSeq[Fac][Pos]);  //将挑选出的组添加到ExtractFamSeq
            FacFamSeq[Fac].erase(FacFamSeq[Fac].begin() + Pos);  //从原序列中删除pos位置的组
        }
    } while (ExtractFamSeq.size() < Len);

    GetSpan(FacFamSeq, JobSeqinFam, FacSpan);
    // 将挑选出来的d个组插入到最好的位置 Insert the extracted Families into the best Positions
    return NEHFam(ExtractFamSeq, JobSeqinFam, FacFamSeq, FacSpan);
}

void GCCEA::GCCEAReInitPop0() //reproduce a solution if age Limit is met
{
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
            for (int i = 0; i < GCCEAPopulation.size(); i++)
            {
                if ((MS == GCCEAPopulation[i].MS) && (TEC == GCCEAPopulation[i].TEC))
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
                GCCEAPopulation.push_back(tempIndi);
            }
        }

    }

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
            for (int i = 0; i < GCCEAPopulation.size(); i++)
            {
                if ((MS == GCCEAPopulation[i].MS) && (TEC == GCCEAPopulation[i].TEC))
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
                GCCEAPopulation.push_back(tempIndi);
            }
        }
    }
}

void GCCEA::Evolution()//randomly setup the mapping relationship, compare with the original relationship, considering span and the number of critical pathes
{
    m_InitTime = Base::GetElapsedProcessTime();
   // int count = 0;
    while (true)
    {
        //count++;
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

            SwapFam(FacFamSeq);

            for (int fac = 0; fac < m_Factories; fac++)
            {
                FamInsert(FacFamSeq[fac], m_RefJobSeqinFamArray[Map1], FacSpan[fac]);
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
                    for (int i = 0; i < GCCEAPopulation.size(); i++)
                    {
                        if ((MS == GCCEAPopulation[i].MS) && (TEC == GCCEAPopulation[i].TEC))
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
                        GCCEAPopulation.push_back(tempIndi);
                    }
                }
            }
            m_Age1[PS]++;
        }//end Fam Evolution

        //Reset m_bFlag2;
        for (int PS = 0; PS < m_RefSize; PS++)
            m_bFlag2[PS] = false;

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

            for (int fac = 0; fac < m_Factories; fac++)
            {
                JobInsert(m_RefFacFamSeqArray[Map2][fac], JobSeqinFam, FacSpan[fac]);
                JobSwap(m_RefFacFamSeqArray[Map2][fac], JobSeqinFam, FacSpan[fac]);
            }

            //int Span = GetSpan(m_RefFamSeqArray[Map2], JobSeqinFam);

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
                    for (int i = 0; i < GCCEAPopulation.size(); i++)
                    {
                        if ((MS == GCCEAPopulation[i].MS) && (TEC == GCCEAPopulation[i].TEC))
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
                        GCCEAPopulation.push_back(tempIndi);
                    }
                }
            }
            m_Age2[PS]++;
        }//end Job Evolution

        //Reset m_bFlag1;
        for (int PS = 0; PS < m_RefSize; PS++)
            m_bFlag1[PS] = false;

        //重新初始化
        GCCEAReInitPop0();

        long CurTime = Base::GetElapsedProcessTime();
        long ElapsedTime = CurTime - m_InitTime;
        if (ElapsedTime >= m_TimeLimit) break;
    }
    //cout << "Iter: " << count << endl;

    //检查
    //for (int PS = 0; PS < GCCEAPopulation.size(); PS++)
    //{
    //	//m_SpeedMatrix = m_ArchiveSpeedVector[PS];
    //	m_SpeedMatrix = GCCEAPopulation[PS].m_SpeedVector;
    //	//得到真正的处理时间 及单位加工能耗
    //	for (int j = 0; j < m_Jobs; j++)
    //	{
    //		for (int i = 0; i < m_Machines; i++)
    //		{
    //			m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
    //			UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
    //		}
    //	}
    //
    //	CheckSol(GCCEAPopulation[PS].m_FacFamSeqArray, GCCEAPopulation[PS].m_JobSeqInFamArray, GCCEAPopulation[PS].MS);
    //	CheckSolTEC(GCCEAPopulation[PS].m_FacFamSeqArray, GCCEAPopulation[PS].m_JobSeqInFamArray, GCCEAPopulation[PS].TEC);
    //}
}


void GCCEA::RunEvolution(int CPUFactor, int Rep, vector<vector<Individual>>& GCCEAFinalAfterRepParetoSet)
{
    ReadInstanceFileNameList("..\\Benchmark\\");
    int RefSize = 5; //Five paramters are Calibrated on 20190724
    int PS1 = 10;
    int PS2 = 10;
    int AgeLimit = 5;
    int DesLen = 2;
    int Instances = 405;
    vector<int> result;

    string FileDirectory = "..\\Result\\";
    ostringstream str;
    str << "GCCEA_" << CPUFactor << "_experiment" << ".txt"; //不同的算法
    ofstream ofile;
    ofile.open(FileDirectory + str.str());


    for (int ins = 0; ins < Instances; ins++)
    {
        srand((unsigned int)time(NULL));
        ReadInstance(ins);

        //XIN
        vector<Individual> FinalGCCEAParetoSet;
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

            SetParameters(RefSize, PS1, PS2, AgeLimit, DesLen, CPUFactor * m_Families * m_Machines);

            InitPop();

            Evolution();

            //非支配解
            GCCEAParetoSet.clear();

            Pareto_relation(GCCEAPopulation);

            //弱化支配解 de-emphasize dominated solutions
            for (int j = 0; j < GCCEAPopulation.size(); j++)
            {
                for (int i = 0; i < GCCEAPopulation.size(); i++)
                {
                    if (GCCEAPopulation[i].flag == 0)
                    {
                        if (GCCEAPopulation[i].pareto_rel[j] == 1)
                        {
                            GCCEAPopulation[i].flag = 999;
                        }
                    }
                }
            }

            for (int i = 0; i < GCCEAPopulation.size(); i++)
            {
                if (GCCEAPopulation[i].flag == 0)
                    GCCEAParetoSet.push_back(GCCEAPopulation[i]);
            }


            //变速
            vector<Individual> Temp1;
            Temp1.clear();
            for (int i = 0; i < GCCEAParetoSet.size(); i++)
            {
                Temp1.push_back(GCCEAParetoSet[i]);
            }

            Speed_mutation(GCCEAParetoSet, Temp1);

            //去除重复
            FinalGCCEAParetoSet.clear();
            for (int i = 0; i < GCCEAParetoSet.size(); i++)
            {
                bool fg = true;
                for (int j = 0; j < FinalGCCEAParetoSet.size(); j++)
                {
                    if ((FinalGCCEAParetoSet[j].MS == GCCEAParetoSet[i].MS) && (FinalGCCEAParetoSet[j].TEC == GCCEAParetoSet[i].TEC))
                    {
                        fg = false;
                        break;
                    }
                }
                if (fg)
                {
                    FinalGCCEAParetoSet.push_back(GCCEAParetoSet[i]);
                }

            }

            long EndTime_Ins = Base::GetElapsedProcessTime();

            //检查
            for (int PS = 0; PS < FinalGCCEAParetoSet.size(); PS++)
            {
                m_SpeedMatrix = FinalGCCEAParetoSet[PS].m_SpeedVector;
                //得到真正的处理时间 及单位加工能耗
                for (int j = 0; j < m_Jobs; j++)
                {
                    for (int i = 0; i < m_Machines; i++)
                    {
                        m_TureJobOpertime[j][i] = static_cast<int>(100 * (m_JobOperPTime[j][i] / m_Speed[m_SpeedMatrix[j][i]]));
                        UnitPEC[j][i] = 4 * m_Speed[m_SpeedMatrix[j][i]] * m_Speed[m_SpeedMatrix[j][i]];
                    }
                }

                //CheckSol(FinalGCCEAParetoSet[PS].m_FacFamSeqArray, FinalGCCEAParetoSet[PS].m_JobSeqInFamArray, FinalGCCEAParetoSet[PS].MS);
                //CheckSolTEC(FinalGCCEAParetoSet[PS].m_FacFamSeqArray, FinalGCCEAParetoSet[PS].m_JobSeqInFamArray, FinalGCCEAParetoSet[PS].TEC);
            }


            for (int PS = 0; PS < FinalGCCEAParetoSet.size(); PS++)
            {
               // cout << FinalGCCEAParetoSet[PS].MS << "\t" << FinalGCCEAParetoSet[PS].TEC << "\t";

                //ofile << FinalGCCEAParetoSet[PS].MS << "," << FinalGCCEAParetoSet[PS].TEC << ",";

                TempAfterRepParetoSet.push_back(FinalGCCEAParetoSet[PS]);
            }
            //cout << r + 1 << "\t" << EndTime_Ins - StartTime_Ins << endl;
            //ofile << endl;

        }//end rep

        //非支配解
        AfterRepParetoSet.clear();

        Pareto_relation(TempAfterRepParetoSet);

        //弱化支配解 de-emphasize dominated solutions
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

        //去除重复
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

        for (int PS = 0; PS < FinalAfterRepParetoSet.size(); PS++)
        {
            // 输出 MS 和 TEC，所有解都放在同一行
            cout << "MS:" << FinalAfterRepParetoSet[PS].MS
                 << " " << "TEC:" << FinalAfterRepParetoSet[PS].TEC << "\t";
            ofile << FinalAfterRepParetoSet[PS].MS << "  "<< FinalAfterRepParetoSet[PS].TEC << ","<< "\t";
            GCCEAFinalAfterRepParetoSet[ins].push_back(FinalAfterRepParetoSet[PS]);
        }
        ofile << endl;  // 每个实例结果结束后换行
        cout << endl;   // 控制台输出换行，仅仅是在实例结束时换行

    }//end ins

    ofile.close();
}