#ifndef DBFGSP_NEW_BASE_H
#define DBFGSP_NEW_BASE_H

#include <vector>
#include <string>
using namespace std;
namespace Base
{
    template<typename T>
    struct Pair
    {
        int dim;
        T value;
    };

    template<typename T>
    struct PairGreater
    {
        bool operator()(Pair<T> a, Pair<T> b)
        {
            return a.value > b.value;
        }
    };

    template<typename T>
    struct PairLess
    {
        bool operator()(Pair<T> a, Pair<T> b)
        {
            return a.value < b.value;
        }
    };

    long GetElapsedProcessTime();//Return Process times in millionSeconds;
    vector<int> GetBestSpan(string ResultsFileName, int Instances, int ExpRep);    //GetBestSpan in a file
    vector<int>
    GetBestSpan(vector<string> ResultsFileNames, int Instances, int ExpRep);//GetBestSpan in a number of files
    double GetRPI(string ResultFileName, vector<int> BestSpan, int ExpRep, string RPIFileName);//GetRPI for a file
    void CombFiles(vector<string> FileNames, string MergeFileName);

    void RecordsAddConfs(vector<int> ParmLevels, int Instances, int Reps, string iFName, string oFName);

    void WriteBatFiles(int Confs, int PCs, string ExeFileName);

    void WriteBatFiles(int Confs, vector<string> IPArray, string ExeFileName);

    void FindIP(string &ip);

    void FindLostConfs(string Dir);

    void WriteRunConfs(string FileName, int Algs = 5, int Reps = 5, int CPULevels = 10);

    void GenerateResultsFileNames();

    void GetRDI(string ResultFileName, vector<int> MaxTTForEachIns, vector<int> MinTTForEachIns, int ExpRep,
                vector<vector<double>> &RDI_Ins_Rep);

    vector<int> GetMinTT_MultipleFiles(vector<string> ResultsFileNames, int Instances, int ExpRep);

    vector<int> GetMaxTT_MultipleFiles(vector<string> ResultsFileNames, int Instances, int ExpRep);

    vector<int> GetMaxTT_OneFile(string ResultsFileName, int Instances, int ExpRep);

    vector<int> GetMinTT_OneFile(string ResultsFileName, int Instances, int ExpRep);
}

#endif //DBFGSP_NEW_BASE_H
