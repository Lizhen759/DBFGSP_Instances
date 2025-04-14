//Base.cpp
#include <windows.h>
#include <minwindef.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include "Base.h"
#include "rand.h"
using namespace std;
long Base::GetElapsedProcessTime()   //得到运行过程时间
{
    FILETIME createTime;
    FILETIME exitTime;
    FILETIME kernelTime;    //核心时间
    FILETIME userTime;

    long ElapsedTime;   //运行时间
    if (GetProcessTimes(GetCurrentProcess(), &createTime, &exitTime, &kernelTime, &userTime) != 0)
    {
        //  Returns total user time.
        SYSTEMTIME userSystemTime;
        if (FileTimeToSystemTime(&userTime, &userSystemTime) != -1)
            ElapsedTime =(userSystemTime.wDay - 1) * 24 * 3600 * 1000
                         + userSystemTime.wHour * 3600 * 1000 +
                         userSystemTime.wMinute * 60 * 1000 +
                         userSystemTime.wSecond * 1000 +
                         userSystemTime.wMilliseconds;
        else
            ElapsedTime = 0;
    }
    else
        ElapsedTime = 0;
    return ElapsedTime;
}

// 找到一个文件里的最佳完工时间 GetBestSpan in a file. bDetail=true, file format includes more information
vector<int> Base::GetBestSpan(string ResultsFileName, int Instances, int ExpRep) //Get Best makespan by an algorithm
{
    vector<int> BestSpan(Instances, INT_MAX);  //每个实例的完工时间
    ifstream ifile;
    ifile.open(ResultsFileName);  //打开文件
    if (ifile.is_open())
    {
        int Span;
        for (int i = 0; i < Instances * ExpRep; i++)
        {
            ifile >> Span;  //从文件读取span
            BestSpan[i / ExpRep] = min(BestSpan[i / ExpRep], Span);//更新bestspan
        }
        ifile.close();
    }
    else
    {
        cout << ResultsFileName << "\t is not open" << endl;
        exit(0);
    }

    return BestSpan;
}

// 找到多个文件中的最佳完工时间 GetBestSpan in a number of files
vector<int> Base::GetBestSpan(vector <string> ResultsFileNames, int Instances, int ExpRep)
{
    vector<vector<int>> SpanArray(ResultsFileNames.size());
    for (int i = 0; i < ResultsFileNames.size(); i++)
        SpanArray[i]=GetBestSpan(ResultsFileNames[i], Instances, ExpRep);  //spanarray存放每个文件最佳时间

    vector<int> BestSpan = SpanArray[0];
    for (int j = 0; j < Instances; j++)
    {
        for (int i = 1; i < ResultsFileNames.size(); i++)
            BestSpan[j] = min(BestSpan[j], SpanArray[i][j]);  //得到所有实例中最好的
    }
    return BestSpan;
}

// 得到一个文件中的平均RPI GetRPI for a file
double Base::GetRPI(string ResultFileName, vector<int> BestSpan, int ExpRep, string RPIFileName)
{
    double AvgRPI = 0;
    ifstream ifile;
    ofstream ofile;
    ifile.open(ResultFileName);
    ofile.open(RPIFileName);
    for (int i = 0; i < BestSpan.size(); i++)
    {
        for (int j = 0; j < ExpRep; j++)
        {
            int Span;
            ifile>> Span;
            double RPI = double(Span - BestSpan[i]) * 100 / BestSpan[i];   //计算每个实例的RPI
            ofile << i << "\t" << RPI << "\t" << endl;  //将RPI写入文件
            AvgRPI += RPI;
        }
    }
    ifile.close();
    ofile.close();
    AvgRPI /= BestSpan.size() * ExpRep;  //平均RPI
    return AvgRPI;
}

//将文件名合并
void Base::CombFiles(vector <string> FileNames, string MergeFileName)
{
    ofstream ofile;
    ofile.open(MergeFileName);
    ifstream ifile;
    vector <string> Temp;
    for (int i = 0; i < FileNames.size(); i++)
    {
        ifile.open(FileNames[i]);  //打开文件i
        while (true)
        {
            char buf[1000];
            ifile.getline(buf, 1000);  //读取文件写入buf
            string str = buf;
            ofile << str << endl;     //将str写入文件
            Temp.push_back(str);      //将str加入Temp
            if (ifile.peek() == EOF)  //若打开空文件则退出
                break;
        }
        ifile.close();
    }
    ofile.close();
}

//记录增加
void Base::RecordsAddConfs(vector<int> ParmLevels, int Instances, int Reps, string iFName, string oFName)
{
    int nConf = 1;
    for (int i = 0; i < ParmLevels.size(); i++)
        nConf *= ParmLevels[i];
    ifstream ifile;
    ofstream ofile;
    ifile.open(iFName);
    ofile.open(oFName);
    for (int No = 0; No < nConf; No++)
    {
        //得到每个参数的标准 Get Level for each parameter
        vector<int> Level(ParmLevels.size());
        int remainder = No;
        for (int i = 0; i < Level.size(); i++)
        {
            int divisor = 1;
            for (int j = i + 1; j < ParmLevels.size(); j++)
                divisor *= ParmLevels[j];
            Level[i] = remainder / divisor;
            remainder %= divisor;
        }

        // 将标准写入 add levels in the records
        for (int Record = 0; Record < Instances * Reps; Record++)
        {
            double RPI;
            int Ins;
            ifile >>Ins >> RPI;
            ofile << Ins << "\t";
            for (int i = 0; i < Level.size(); i++)
            {
                ofile << Level[i] << "\t";
            }
            ofile << RPI << endl;
        }
    }
    ifile.close();
    ofile.close();
}

void Base::WriteBatFiles(int Confs, int PCs, string ExeFileName)
{
    srand(1971);
    vector<int> RunSeq(Confs, 0);
    for (int i = 0; i < Confs; i++)
        RunSeq[i] = i;
    shuffle(RunSeq.begin(), RunSeq.end(), rand_generator());
    vector <vector<int>> ConfPerPC(PCs);
    int Index_PC = 0;
    for (int i = 0; i < RunSeq.size(); i++)
    {
        Index_PC %= PCs;
        ConfPerPC[Index_PC].push_back(RunSeq[i]);
        Index_PC++;
    }

    for (int index = 0; index < ConfPerPC.size(); index++)
    {
        ostringstream FName;
        FName << "Run_" << index << ".bat";
        ofstream ofile;
        ofile.open(FName.str());
        for (int i = 0; i < ConfPerPC[index].size(); i++)
            ofile << ExeFileName << "\t" << ConfPerPC[index][i] << endl;
        ofile.close();
    }
}

void Base::WriteBatFiles(int Confs, vector<string> IPArray, string ExeFileName)
{
    srand(1971);
    vector<int> RunSeq(Confs, 0);
    for (int i = 0; i < Confs; i++)
        RunSeq[i] = i;
    shuffle(RunSeq.begin(), RunSeq.end(), rand_generator());
    ofstream ofile;
    ofile.open("Run.bat");
    for (int i = 0; i < Confs; i++)
        ofile << ExeFileName << "\t" << RunSeq[i] << "\t" << IPArray[i%IPArray.size()] << endl;
    ofile.close();
}

void Base::WriteRunConfs(string FileName,int Algs,int Reps,int CPULevels)
{
    ofstream ofile;
    ofile.open(FileName);
    for (int a = 0; a < Algs; a++)
    {
        for (int rep = 0; rep < Reps; rep++)
        {
            for (int TimeLevel = 9; TimeLevel < CPULevels; TimeLevel++)
            {
                ofile << a << "\t" << (TimeLevel + 1) * 50 << "\t" << rep  << endl;
            }
        }
    }
    ofile << -1 << "\t" << -1 << "\t" << -1 << endl;
    ofile.close();
}

void Base::GenerateResultsFileNames()
{
    ofstream ofile;
    ofile.open("ResultFileName.txt");
    int Algs = 5;
    int Reps = 5;
    int CPULevels = 10;
    for (int a = 0; a < Algs; a++)
    {
        for (int rep = 0; rep < Reps; rep++)
        {
            for (int TimeLevel = 0; TimeLevel < CPULevels; TimeLevel++)
            {
                ostringstream ostr;
                ostr<<"CCEA"<<a<<"_"<<(TimeLevel + 1) * 50 <<"_"<< rep<<".txt";
                ofile << ostr.str() << endl;
            }
        }
    }
    ofile.close();
}

//void Base::FindIP(string &ip)
//{
//    WORD v = MAKEWORD(1, 1);
//    WSADATA wsaData;
//    WSAStartup(v, &wsaData); // 加载套接字库
//
//    struct hostent *phostinfo = gethostbyname("");
//    char *p = inet_ntoa(*((struct in_addr *)(*phostinfo->h_addr_list)));
//    int xx = sizeof(p);
//    ip = p;
//    WSACleanup();
//}

void Base::FindLostConfs(string Dir)
{
    int Confs = 300;
    ofstream ofile;
    ofile.open(Dir+"Confs.txt");
    for (int i = 0; i < Confs; i++)
    {
        ostringstream str;
        str << "EA_" << i << ".txt";
        ifstream ifile;
        ifile.open(Dir+str.str());
        if (ifile.is_open())
        {
            ifile.close();
        }
        else
        {
            cout << i << endl;
            ofile << i << endl;
        }
    }
    ofile << -1 << endl;//End of the congfig file
    ofile.close();
}

//Get minimum TT for each instance by an algorithm
vector<int> Base::GetMinTT_OneFile(string ResultsFileName, int Instances, int ExpRep)
{
    vector<int> MinTTForEachInstance_OneFile(Instances, INT_MAX);
    ifstream ifile;
    ifile.open(ResultsFileName);
    if (ifile.is_open())
    {
        int TT;
        for (int i = 0; i < Instances * ExpRep; i++)
        {
            ifile >> TT;
            MinTTForEachInstance_OneFile[i / ExpRep] = min(MinTTForEachInstance_OneFile[i / ExpRep], TT);
        }
        ifile.close();
    }
    else
    {
        cout << ResultsFileName << "\t is not open" << endl;
        exit(0);
    }
    return MinTTForEachInstance_OneFile;
}

//Get maximum TT for each instance by an algorithm
vector<int> Base::GetMaxTT_OneFile(string ResultsFileName, int Instances, int ExpRep)
{
    vector<int> MaxTTForEachInstance_OneFile(Instances, INT_MIN);
    ifstream ifile;
    ifile.open(ResultsFileName);
    if (ifile.is_open())
    {
        int TT;
        for (int i = 0; i < Instances * ExpRep; i++)
        {
            ifile >> TT;
            MaxTTForEachInstance_OneFile[i / ExpRep] = max(MaxTTForEachInstance_OneFile[i / ExpRep], TT);
        }
        ifile.close();
    }
    else
    {
        cout << ResultsFileName << "\t is not open" << endl;
        exit(0);
    }
    return MaxTTForEachInstance_OneFile;
}

//Get minimum TT for each instance by all algorithms
vector<int> Base::GetMinTT_MultipleFiles(vector <string> ResultsFileNames, int Instances, int ExpRep)
{
    vector<vector<int>> TTArrayForEachAlg(ResultsFileNames.size());
    for (int i = 0; i < ResultsFileNames.size(); i++)
        TTArrayForEachAlg[i] = GetMinTT_OneFile(ResultsFileNames[i], Instances, ExpRep);

    vector<int> BestTTForEachInstance_MultipleFile = TTArrayForEachAlg[0];
    for (int j = 0; j < Instances; j++)
    {
        for (int i = 1; i < ResultsFileNames.size(); i++)
            BestTTForEachInstance_MultipleFile[j] = min(BestTTForEachInstance_MultipleFile[j], TTArrayForEachAlg[i][j]);
    }
    return BestTTForEachInstance_MultipleFile;
}

//Get maximum TT for each instance by all algorithms
vector<int> Base::GetMaxTT_MultipleFiles(vector <string> ResultsFileNames, int Instances, int ExpRep)
{
    vector<vector<int>> TTArrayForEachAlg(ResultsFileNames.size());
    for (int i = 0; i < ResultsFileNames.size(); i++)
        TTArrayForEachAlg[i] = GetMaxTT_OneFile(ResultsFileNames[i], Instances, ExpRep);

    vector<int> MaxTTForEachInstance_MultipleFile = TTArrayForEachAlg[0];
    for (int j = 0; j < Instances; j++)
    {
        for (int i = 1; i < ResultsFileNames.size(); i++)
            MaxTTForEachInstance_MultipleFile[j] = max(MaxTTForEachInstance_MultipleFile[j], TTArrayForEachAlg[i][j]);
    }
    return MaxTTForEachInstance_MultipleFile;
}



void Base::GetRDI(string ResultFileName, vector<int> MaxTTForEachIns, vector<int> MinTTForEachIns, int ExpRep, vector<vector<double>>& RDI_Ins_Rep)
{
    ifstream ifile;
    ifile.open(ResultFileName);
    double RDI = 0;
    int TT = 0;

    for (int i = 0; i < MaxTTForEachIns.size(); i++)
    {
        for (int r = 0; r < ExpRep; r++)
        {
            ifile >> TT;
            if (MaxTTForEachIns[i] == MinTTForEachIns[i])
                RDI = 0;
            else
                RDI = double(TT - MinTTForEachIns[i]) / (MaxTTForEachIns[i] - MinTTForEachIns[i]);//RDI计算公式。
            //RDI = double (TT - MinTTForEachIns[i]) * 100 / (MaxTTForEachIns[i] - MinTTForEachIns[i]);//RDI计算公式。
            RDI_Ins_Rep[i][r] = RDI;
        }
    }
    ifile.close();
}