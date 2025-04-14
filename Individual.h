#pragma once
#include "Problem.h"
#include "stdlib.h"


using namespace std;

#define maxpop  10000

class Individual
{
public:
    //新
    vector<vector<int>> m_FacFamSeqArray;
    vector<vector<int>> m_JobSeqInFamArray;

    int MS;				   //the objective: makespan
    float TEC;			   //the objective: total energy cost objective

    float normalMS;       //归一化后的
    float normalTEC;

    float convergence_ind;//value of convergence indicator
    float distribution_ind[maxpop];//value of distribution indicator

    vector<vector<int>> m_SpeedVector;  //速度向量

    int	flag;  //判断标志

    int pareto_rel[maxpop];

    Individual& operator=(const Individual& other)
    {
        if (this != &other)
        {
            m_FacFamSeqArray = other.m_FacFamSeqArray;
            m_JobSeqInFamArray = other.m_JobSeqInFamArray;
            MS = other.MS;
            TEC = other.TEC;
            normalMS = other.normalMS;
            normalTEC = other.normalTEC;
            convergence_ind = other.convergence_ind;
            m_SpeedVector = other.m_SpeedVector;
            //flag = other.flag;
        }
        return *this;
    }

public:
    Individual();    //constructor
    //Individual(int k);
    ~Individual();    //deconstructor

    bool cmp_convergence(const Individual& x, const Individual& y);

    //归一化
    void Normalize(vector<Individual>& Pop, int& nadirpointMS, float& nadirpointTEC, int& idealpointMS, float& indealpointTEC);

    //计算收敛性指标
    void calc_convergence_ind(vector<Individual>& Pop);

    //计算分布性指标
    static void calc_distribution_ind(vector<Individual>& Pop);

    void getObjectives();

    static float calc_distribution_ind_QCCEA(vector<Individual> &CCEAPopulation, float normalMS, float normalTEC);

    static float calculate_norm(float normalMS, float normalTEC);
};