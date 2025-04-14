#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <valarray>
#include "Problem.h"
using namespace std;

Problem::Problem()
{
}

Problem::~Problem()
{
}


//得到实例文件名
void Problem::ReadInstanceFileNameList(string Dir)
{
    m_InstanceFileNameList.clear();
    m_InstanceNameList.clear();
    ifstream ifile;
    ifile.open(Dir + "\\" + "FileNameList.txt");
    while (true)
    {
        int x;
        string FName;
        ifile >> x >> FName;
        if (ifile.peek() != EOF)
        {
            m_InstanceFileNameList.push_back(Dir + "\\" + FName);
            m_InstanceNameList.push_back(FName);
        }
        else
            break;
    }
    ifile.close();
}

void Problem::ReadInstance(int InsNo)
{
    ifstream ifile;
    ifile.open(m_InstanceFileNameList[InsNo]);
    string str;
    int data;

    // Read Configuration
    ifile >> str >> m_Factories;
   // std::cout << "Factories: " << m_Factories << std::endl;

    ifile >> str >> m_Families;
   // std::cout << "Families: " << m_Families << std::endl;

    ifile >> str >> m_Machines;
  //  std::cout << "Machines: " << m_Machines << std::endl;

    ifile >> str >> data;  // Setup time type
  //  std::cout << "Setup Type: " << data << std::endl;

    // Read Number of Jobs in each Family
    vector<int> NumbofJobsinEachFamily(m_Families);  // Read the number of jobs in each family
    for (int i = 0; i < 6; i++)
        ifile >> str;  // Skip unnecessary lines if needed

    for (int Fam = 0; Fam < m_Families; Fam++) {
        ifile >> NumbofJobsinEachFamily[Fam];
      //  std::cout << "Number of Jobs in Family " << Fam << ": " << NumbofJobsinEachFamily[Fam] << std::endl;
    }

    // Read Jobs in each Family
    for (int i = 0; i < 4; i++)
        ifile >> str;  // Skip unnecessary lines

    m_JobsInEachFamily.clear();
    m_JobsInEachFamily.resize(m_Families);
    for (int Fam = 0; Fam < m_Families; Fam++)
        m_JobsInEachFamily[Fam].resize(NumbofJobsinEachFamily[Fam]);

    for (int Fam = 0; Fam < m_Families; Fam++)
        for (int i = 0; i < m_JobsInEachFamily[Fam].size(); i++) {
            ifile >> m_JobsInEachFamily[Fam][i];  // Read jobs for each family
          //  std::cout << "Job " << m_JobsInEachFamily[Fam][i] << " in Family " << Fam << std::endl;
        }

    // Read Total Number of Jobs
    for (int i = 0; i < 4; i++)
        ifile >> str;  // Skip unnecessary lines

    ifile >> m_Jobs;  // Read total number of jobs
   // std::cout << "Total Jobs: " << m_Jobs << std::endl;

    // Read Processing times of Jobs
    for (int i = 0; i < 4; i++)
        ifile >> str;  // Skip unnecessary lines

    m_JobOperPTime.clear();
    m_JobOperPTime.resize(m_Jobs);
    for (int j = 0; j < m_Jobs; j++)
        m_JobOperPTime[j].resize(m_Machines);

    for (int j = 0; j < m_Jobs; j++)
        for (int m = 0; m < m_Machines; m++) {
            ifile >> m_JobOperPTime[j][m];  // Read processing time of job j on machine m
//            std::cout << "Processing Time of Job " << j << " on Machine " << m << ": " << m_JobOperPTime[j][m] << std::endl;
        }

    // Read Setup times between Families
    for (int i = 0; i < 10; i++)
        ifile >> str;  // Skip unnecessary lines

    m_SetupTime.clear();
    m_SetupTime.resize(m_Machines);
    for (int mac = 0; mac < m_Machines; mac++) {
        m_SetupTime[mac].resize(m_Families);
        for (int Fam = 0; Fam < m_Families; Fam++)
            m_SetupTime[mac][Fam].resize(m_Families);
    }

    for (int mac = 0; mac < m_Machines; mac++) {
        for (int i = 0; i < 3; i++)
            ifile >> str;  // Skip unnecessary lines

        for (int Fam1 = 0; Fam1 < m_Families; Fam1++)
            for (int Fam2 = 0; Fam2 < m_Families; Fam2++) {
                ifile >> m_SetupTime[mac][Fam1][Fam2];  // Read setup times
//                std::cout << "Setup Time between Family " << Fam1 << " and Family " << Fam2
//                          << " on Machine " << mac << ": " << m_SetupTime[mac][Fam1][Fam2] << std::endl;
            }
    }

    m_Speed.clear();
    m_Speed.resize(3);
    m_Speed[0] = 1.0;
    m_Speed[1] = 1.5;
    m_Speed[2] = 2.0;

    UnitIdleEC = 1.0;
    UnitSetupEC = 0.50;
    ifile.close();
}


void Problem::GetJobTotalPTime()
{
    this->m_JobTotalPTime.clear();
    this->m_JobTotalPTime.resize(this->m_Jobs, 0);
    m_AllJobTotalPTime = 0;
    for (int j = 0; j < this->m_JobTotalPTime.size(); j++)
    {
        for (int m = 0; m < this->m_Machines; m++)
        {
            this->m_JobTotalPTime[j] += this->m_JobOperPTime[j][m];
            m_AllJobTotalPTime += this->m_JobOperPTime[j][m];
        }
    }
}


void Problem::GetJobWeightTotalPTime()
{
    this->m_JobWeightTotalPTime.clear();
    this->m_JobWeightTotalPTime.resize(this->m_Jobs, 0);
    for (int j = 0; j < this->m_JobWeightTotalPTime.size(); j++)
        for (int m = 0; m < this->m_Machines; m++)
            this->m_JobWeightTotalPTime[j] += this->m_JobOperPTime[j][m] * (this->m_Machines - m);//一个工件在前面机器上的PTime权重大。
}

//得到组内工件的加工时间之和
void Problem::GetFamTotalPTime()
{
    this->m_FamTotalPTime.clear();
    this->m_FamTotalPTime.resize(this->m_JobsInEachFamily.size(), 0);
    for (int Fam = 0; Fam < this->m_JobsInEachFamily.size(); Fam++)
        for (int j = 0; j < this->m_JobsInEachFamily[Fam].size(); j++)
            this->m_FamTotalPTime[Fam] += this->m_JobTotalPTime[this->m_JobsInEachFamily[Fam][j]];
}


void Problem::GetFamWeightTotalPTime()
{
    this->m_FamWeightTotalPTime.clear();
    this->m_FamWeightTotalPTime.resize(this->m_JobsInEachFamily.size(), 0);
    for (int Fam = 0; Fam < this->m_JobsInEachFamily.size(); Fam++)
        for (int j = 0; j < this->m_JobsInEachFamily[Fam].size(); j++)
            this->m_FamWeightTotalPTime[Fam] += this->m_JobWeightTotalPTime[this->m_JobsInEachFamily[Fam][j]];
}


void Problem::GetFamTotalPTimeOnFirstMachine()
{
    this->m_FamTotalPTimeOnFirstMachine.clear();
    this->m_FamTotalPTimeOnFirstMachine.resize(this->m_JobsInEachFamily.size(), 0);
    for (int Fam = 0; Fam < this->m_JobsInEachFamily.size(); Fam++)
        for (int j = 0; j < this->m_JobsInEachFamily[Fam].size(); j++)
            this->m_FamTotalPTimeOnFirstMachine[Fam] += this->m_JobOperPTime[this->m_JobsInEachFamily[Fam][j]][0];
}


void Problem::GetFamTotalPTimeOnLastMachine()
{
    this->m_FamTotalPTimeOnLastMachine.clear();
    this->m_FamTotalPTimeOnLastMachine.resize(this->m_JobsInEachFamily.size(), 0);
    for (int Fam = 0; Fam < this->m_JobsInEachFamily.size(); Fam++)
        for (int j = 0; j < this->m_JobsInEachFamily[Fam].size(); j++)
            this->m_FamTotalPTimeOnLastMachine[Fam] += this->m_JobOperPTime[this->m_JobsInEachFamily[Fam][j]][this->m_Machines - 1];
}


void Problem::GetFamSumSetupTime()
{
    this->m_FamSumSetupTime.clear();
    this->m_FamSumSetupTime.resize(this->m_Families);
    for (int Fam = 0; Fam < this->m_FamSumSetupTime.size(); Fam++)
        this->m_FamSumSetupTime[Fam].resize(this->m_Families);
    for (int CurFam = 0; CurFam < this->m_FamSumSetupTime.size(); CurFam++)
        for (int PreFam = 0; PreFam < this->m_FamSumSetupTime.size(); PreFam++)
        {
            if (CurFam == PreFam)
            {
                this->m_FamSumSetupTime[PreFam][CurFam] = RAND_MAX;
            }
            else {

                this->m_FamSumSetupTime[PreFam][CurFam] = 0;
                for (int m = 0; m < this->m_Machines; m++)
                    this->m_FamSumSetupTime[PreFam][CurFam] = this->m_FamSumSetupTime[PreFam][CurFam] + this->m_SetupTime[m][PreFam][CurFam];
            }
        }
}

void Problem::GetFamMaxSetupTime()
{
    this->m_FamMaxSetupTime.clear();
    this->m_FamMaxSetupTime.resize(this->m_Families);
    for (int Fam = 0; Fam < this->m_FamMaxSetupTime.size(); Fam++)
        this->m_FamMaxSetupTime[Fam].resize(this->m_Families);
    for (int CurFam = 0; CurFam < this->m_FamMaxSetupTime.size(); CurFam++)
        for (int PreFam = 0; PreFam < this->m_FamMaxSetupTime.size(); PreFam++)
        {
            this->m_FamMaxSetupTime[PreFam][CurFam] = this->m_SetupTime[0][PreFam][CurFam];
            for (int m = 1; m < this->m_Machines; m++)
                if (this->m_FamMaxSetupTime[PreFam][CurFam] < this->m_SetupTime[m][PreFam][CurFam])
                    this->m_FamMaxSetupTime[PreFam][CurFam] = this->m_SetupTime[m][PreFam][CurFam];
        }
}


void Problem::GetFamAvgMaxSetupTime()
{
    this->m_FamAvgMaxSetupTime.clear();
    this->m_FamAvgMaxSetupTime.resize(this->m_Families, 0);
    for (int CurFam = 0; CurFam < this->m_FamAvgMaxSetupTime.size(); CurFam++)
    {
        for (int PreFam = 0; PreFam < this->m_FamAvgMaxSetupTime.size(); PreFam++)
            this->m_FamAvgMaxSetupTime[CurFam] += this->m_FamMaxSetupTime[PreFam][CurFam];
        this->m_FamAvgMaxSetupTime[CurFam] /= this->m_FamAvgSetupTime.size();
    }
}

//得到组平均准备时间
void Problem::GetFamAvgSetupTime()
{
    this->m_FamAvgSetupTime.clear();
    this->m_FamAvgSetupTime.resize(this->m_Families, 0);
    //当前组与所有组在所有机器上准备时间的和 / 组数
    for (int CurFam = 0; CurFam < this->m_FamAvgSetupTime.size(); CurFam++)
    {
        for (int PreFam = 0; PreFam < this->m_FamAvgSetupTime.size(); PreFam++)
            for (int m = 0; m < this->m_Machines; m++)
                this->m_FamAvgSetupTime[CurFam] += this->m_SetupTime[m][PreFam][CurFam];
        this->m_FamAvgSetupTime[CurFam] /= this->m_FamAvgSetupTime.size();
    }
}



/****************************QCCEA************************************/


//得到组内工件的加工时间之和
void Problem::GetFamTotalPTime_QCCEA()
{
    this->m_FamTotalPTime.clear();
    this->m_FamTotalPTime.resize(this->m_JobsInEachFamily.size(), 0);
    for (int Fam = 0; Fam < this->m_JobsInEachFamily.size(); Fam++)
        for (int j = 0; j < this->m_JobsInEachFamily[Fam].size(); j++)
            this->m_FamTotalPTime[Fam] += this->m_JobTotalPTime[this->m_JobsInEachFamily[Fam][j]];

    // ========== 归一化 ==========
    double min_val = *min_element(this->m_FamTotalPTime.begin(), this->m_FamTotalPTime.end());
    double max_val = *max_element(this->m_FamTotalPTime.begin(), this->m_FamTotalPTime.end());
    for (double &val : this->m_FamTotalPTime) {
        val = (max_val > min_val) ? (val - min_val) / (max_val - min_val) : 0.0;
    }
}


void Problem::GetFamWeightTotalPTime_QCCEA()
{
    this->m_FamWeightTotalPTime.clear();
    this->m_FamWeightTotalPTime.resize(this->m_JobsInEachFamily.size(), 0);
    for (int Fam = 0; Fam < this->m_JobsInEachFamily.size(); Fam++)
        for (int j = 0; j < this->m_JobsInEachFamily[Fam].size(); j++)
            this->m_FamWeightTotalPTime[Fam] += this->m_JobWeightTotalPTime[this->m_JobsInEachFamily[Fam][j]];

    // ========== 归一化 ==========
    double min_val = *min_element(this->m_FamWeightTotalPTime.begin(), this->m_FamWeightTotalPTime.end());
    double max_val = *max_element(this->m_FamWeightTotalPTime.begin(), this->m_FamWeightTotalPTime.end());

    for (double &val : this->m_FamWeightTotalPTime) {
        val = (max_val > min_val) ? (val - min_val) / (max_val - min_val) : 0.0;
    }
}


void Problem::GetFamTotalPTimeOnFirstMachine_QCCEA()
{
    this->m_FamTotalPTimeOnFirstMachine.clear();
    this->m_FamTotalPTimeOnFirstMachine.resize(this->m_JobsInEachFamily.size(), 0);
    for (int Fam = 0; Fam < this->m_JobsInEachFamily.size(); Fam++)
        for (int j = 0; j < this->m_JobsInEachFamily[Fam].size(); j++)
            this->m_FamTotalPTimeOnFirstMachine[Fam] += this->m_JobOperPTime[this->m_JobsInEachFamily[Fam][j]][0];

    // ========== 归一化 ==========
    double min_val = *min_element(this->m_FamTotalPTimeOnFirstMachine.begin(), this->m_FamTotalPTimeOnFirstMachine.end());
    double max_val = *max_element(this->m_FamTotalPTimeOnFirstMachine.begin(), this->m_FamTotalPTimeOnFirstMachine.end());

    for (double &val : this->m_FamTotalPTimeOnFirstMachine) {
        val = (max_val > min_val) ? (val - min_val) / (max_val - min_val) : 0.0;
    }
}


void Problem::GetFamTotalPTimeOnLastMachine_QCCEA()
{
    this->m_FamTotalPTimeOnLastMachine.clear();
    this->m_FamTotalPTimeOnLastMachine.resize(this->m_JobsInEachFamily.size(), 0);
    for (int Fam = 0; Fam < this->m_JobsInEachFamily.size(); Fam++)
        for (int j = 0; j < this->m_JobsInEachFamily[Fam].size(); j++)
            this->m_FamTotalPTimeOnLastMachine[Fam] += this->m_JobOperPTime[this->m_JobsInEachFamily[Fam][j]][this->m_Machines - 1];
    // ========== 归一化 ==========
    double min_val = *min_element(this->m_FamTotalPTimeOnLastMachine.begin(), this->m_FamTotalPTimeOnLastMachine.end());
    double max_val = *max_element(this->m_FamTotalPTimeOnLastMachine.begin(), this->m_FamTotalPTimeOnLastMachine.end());

    for (double &val : this->m_FamTotalPTimeOnLastMachine) {
        val = (max_val > min_val) ? (val - min_val) / (max_val - min_val) : 0.0;
    }
}


//得到组平均准备时间
void Problem::GetFamAvgSetupTime_QCCEA()
{
    this->m_FamAvgSetupTime.clear();
    this->m_FamAvgSetupTime.resize(this->m_Families, 0);
    //当前组与所有组在所有机器上准备时间的和 / 组数
    for (int CurFam = 0; CurFam < this->m_FamAvgSetupTime.size(); CurFam++)
    {
        for (int PreFam = 0; PreFam < this->m_FamAvgSetupTime.size(); PreFam++)
            for (int m = 0; m < this->m_Machines; m++)
                this->m_FamAvgSetupTime[CurFam] += this->m_SetupTime[m][PreFam][CurFam];
        this->m_FamAvgSetupTime[CurFam] /= this->m_FamAvgSetupTime.size();
    }
    // ========== 归一化 ==========
    double min_val = *min_element(this->m_FamAvgSetupTime.begin(), this->m_FamAvgSetupTime.end());
    double max_val = *max_element(this->m_FamAvgSetupTime.begin(), this->m_FamAvgSetupTime.end());

    for (double &val : this->m_FamAvgSetupTime) {
        val = (max_val > min_val) ? (val - min_val) / (max_val - min_val) : 0.0;
    }
}

void Problem::GetFamTotalSkewness()
{
    this->m_FamTotalSkewness.clear();
    this->m_FamTotalSkewness.resize(this->m_JobsInEachFamily.size(), 0);

    // 计算每个组的偏度
    for (int Fam = 0; Fam < this->m_JobsInEachFamily.size(); Fam++)
    {
        vector<double> jobProcessingTimes;

        // 获取当前组的所有工件的加工时间
        for (int j = 0; j < this->m_JobsInEachFamily[Fam].size(); j++)
        {
            int jobId = this->m_JobsInEachFamily[Fam][j];
            jobProcessingTimes.push_back(this->m_JobTotalPTime[jobId]);
        }

        // 计算当前组的加工时间的偏度
        this->m_FamTotalSkewness[Fam] = CalculateSkewness(jobProcessingTimes);
    }

    // ========================= 归一化偏度 =========================
    double min_skew = *min_element(this->m_FamTotalSkewness.begin(), this->m_FamTotalSkewness.end());
    double max_skew = *max_element(this->m_FamTotalSkewness.begin(), this->m_FamTotalSkewness.end());

    for (double &val : this->m_FamTotalSkewness) {
        val = (max_skew > min_skew) ? (val - min_skew) / (max_skew - min_skew) : 0.0;
    }
}


// 计算偏度的函数
double Problem::CalculateSkewness(const vector<double>& data)
{
    int n = data.size();
    if (n < 3) return 0.0;  // 数据少于3个点时，偏度不可计算

    // 计算均值
    double mean = 0.0;
    double variance = 0.0;
    for (int i = 0; i < n; i++) {
        mean += data[i];
    }
    mean /= n;

    // 计算方差
    for (int i = 0; i < n; i++) {
        variance += (data[i] - mean) * (data[i] - mean);
    }
    variance /= n;
    double stddev = sqrt(variance);

    // 如果标准差为0，返回偏度为0
    if (stddev == 0) return 0.0;

    // 计算偏度
    double skewness = 0.0;
    for (int i = 0; i < n; i++) {
        skewness += pow((data[i] - mean) / stddev, 3);
    }

    skewness *= n / ((n - 1) * (n - 2));  // 偏度标准化
    return skewness;

}
