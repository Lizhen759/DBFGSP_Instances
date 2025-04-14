#include <unordered_set>
#include "Heu.h"

Heu::Heu()
{
}

Heu::~Heu()
{
}

/**
 * 对组内的工件根据某一属性排序
 * @param Factor
 * @param SortMethod
 * @param JobSeqinFam
 */
void Heu::SortJobsinFam(int Factor, int SortMethod, vector<vector<int>>& JobSeqinFam)
{
    JobSeqinFam = this->m_JobsInEachFamily; //initialize memeory
    for (int Fam = 0; Fam < JobSeqinFam.size(); Fam++)
    {
        Pair<int>* ch = new Pair <int>[JobSeqinFam[Fam].size()];
        for (int j = 0; j < JobSeqinFam[Fam].size(); j++)
        {
            ch[j].dim = JobSeqinFam[Fam][j];
            switch (Factor)
            {
                case 0:
                    ch[j].value = this->m_JobTotalPTime[JobSeqinFam[Fam][j]];
                break;
                case 1:
                    ch[j].value = this->m_JobOperPTime[JobSeqinFam[Fam][j]][0];
                    break;
                case 2:
                    ch[j].value = -this->m_JobOperPTime[JobSeqinFam[Fam][j]][this->m_Machines - 1];
                    break;
                case 3:
                    ch[j].value = this->m_JobOperPTime[JobSeqinFam[Fam][j]][0] - this->m_JobOperPTime[JobSeqinFam[Fam][j]][this->m_Machines - 1];
                    break;
                default:break;
            }
        }

        if (SortMethod == 0)
            sort(ch, ch + JobSeqinFam[Fam].size(), PairGreater<int>());
        else
            sort(ch, ch + JobSeqinFam[Fam].size(), PairLess<int>());
        for (int j = 0; j < JobSeqinFam[Fam].size(); j++)
            JobSeqinFam[Fam][j] = ch[j].dim;
        delete[] ch;
    }
}

void Heu::SortJobsInFam(int SortMethod, vector<vector<int>>& JobSeqInFam) //0:LPT,   1:SPT,  2:JobWeightTotalPTime,
{

    JobSeqInFam = this->m_JobsInEachFamily; //initialize memeory
    if (SortMethod == 0) //LPT
    {
        for (auto& JobSeq : JobSeqInFam)
        {
            sort(begin(JobSeq), end(JobSeq), [this](int j1, int j2)
            {
                return this->m_JobTotalPTime[j1] > this->m_JobTotalPTime[j2];
            });
        }
    }

    if (SortMethod == 1) //SPT
    {
        for (auto& JobSeq : JobSeqInFam)
        {
            sort(begin(JobSeq), end(JobSeq), [this](int j1, int j2)
            {
                return this->m_JobTotalPTime[j1] < this->m_JobTotalPTime[j2];
            });
        }
    }

    if (SortMethod == 2) //in non-increasing order JobWeightTotalPTime
    {
        for (auto& JobSeq : JobSeqInFam)
        {
            sort(begin(JobSeq), end(JobSeq), [this](int j1, int j2)
            {
                return this->m_JobWeightTotalPTime[j1] > this->m_JobWeightTotalPTime[j2];
            });
        }
    }

    if (SortMethod == 3) //in non-increasing order (job ptime on the first machine)
    {
        for (auto& JobSeq : JobSeqInFam)
        {
            sort(begin(JobSeq), end(JobSeq), [this](int j1, int j2)
            {
                return this->m_JobOperPTime[j1][0] > this->m_JobOperPTime[j2][0];
            });
        }
    }

    if (SortMethod == 4) //in non-increasing order (job ptime on the last machine)
    {
        for (auto& JobSeq : JobSeqInFam)
        {
            sort(begin(JobSeq), end(JobSeq), [this](int j1, int j2)
            {
                return this->m_JobOperPTime[j1][this->m_Machines - 1] > this->m_JobOperPTime[j2][this->m_Machines - 1];
            });
        }
    }

    if (SortMethod == 5) //in non-increasing order (job ptime on the first machine - job ptime on the last machine)
    {
        for (auto& JobSeq : JobSeqInFam)
        {
            sort(begin(JobSeq), end(JobSeq), [this](int j1, int j2)
            {
                return this->m_JobOperPTime[j1][0] - this->m_JobOperPTime[j1][this->m_Machines - 1] >
                       this->m_JobOperPTime[j2][0] - this->m_JobOperPTime[j2][this->m_Machines - 1];
            });
        }
    }

    if (SortMethod == 6) //黄颖颖师姐文章里面针对组成产品的工件的索引函数
    {
        for (int Fam = 0; Fam < JobSeqInFam.size(); Fam++)
        {
            Pair<int>* ch = new Pair<int>[JobSeqInFam[Fam].size()];
            vector<int> temp(m_Families, 0);

            for (int j = 0; j < JobSeqInFam[Fam].size(); j++)
            {
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


/**
 * 对单个组内的工件根据某一属性（所有场景下）排序
 * 0:LPT 1:SPT
 * @param Factor
 * @param SortMethod
 * @param JobSeqinFam
 */
void Heu::SortJobsinOtherFam(int Factor, int SortMethod, vector<vector<int>> &JobSeqinFam)
{
    for (int Fam = 0; Fam < JobSeqinFam.size(); Fam++)
    {
        Pair<int> *ch = new Pair<int>[JobSeqinFam[Fam].size()];
        for (int j = 0; j < JobSeqinFam[Fam].size(); j++)
        {
            ch[j].dim = JobSeqinFam[Fam][j];  //.dim表示实际的工件编号
            switch (Factor)
            {
                case 0:
                    ch[j].value = this->m_JobTotalPTime[JobSeqinFam[Fam][j]];
                    break; //工件总加工时间
                case 1:
                {
                    vector<int> JPT(m_Jobs, 0);
                        JPT[JobSeqinFam[Fam][j]] += m_JobOperPTime[JobSeqinFam[Fam][j]][this->m_Machines - 1];
                    ch[j].value = JPT[JobSeqinFam[Fam][j]];
                    break; //工件在最后一台机器上的加工时间
                }
                default:
                    break;
            }
        }
        if (SortMethod == 0)
            sort(ch, ch + JobSeqinFam[Fam].size(), PairGreater<int>());  //根据属性值大的排序
        else
            sort(ch, ch + JobSeqinFam[Fam].size(), PairLess<int>());     //根据属性值小的排序
        for (int j = 0; j < JobSeqinFam[Fam].size(); j++)
            JobSeqinFam[Fam][j] = ch[j].dim;  //将排好的工件赋给JobSeqinFam
        delete[] ch;
    }
}

/**
 * 对组根据某一属性排序
 * @param Factor
 * @param SortMethod
 * @param FamPrmu
 */
void Heu::SortFam(int Factor, int SortMethod, vector<int>& FamPrmu)
{
    FamPrmu.clear();
    FamPrmu.resize(this->m_Families);
    for (int Fam = 0; Fam < this->m_Families; Fam++)
        FamPrmu[Fam] = Fam;

    Pair<double>* ch = new Pair<double>[FamPrmu.size()];
    for (int Fam = 0; Fam < FamPrmu.size(); Fam++)
    {
        ch[Fam].dim = FamPrmu[Fam];
        switch (Factor)
        {
            case 0:
                ch[Fam].value = this->m_FamTotalPTime[FamPrmu[Fam]];
                break;
            case 1:
                ch[Fam].value = this->m_FamTotalPTime[FamPrmu[Fam]] + this->m_FamAvgSetupTime[FamPrmu[Fam]];
            break;
            case 2:
                ch[Fam].value = this->m_FamTotalPTime[FamPrmu[Fam]] + this->m_FamAvgMaxSetupTime[FamPrmu[Fam]];
                break;
            case 3:
                ch[Fam].value = this->m_FamWeightTotalPTime[FamPrmu[Fam]];
                break;
            case 4:
                ch[Fam].value = this->m_FamTotalPTimeOnFirstMachine[FamPrmu[Fam]];
                break;
            case 5:
                ch[Fam].value = this->m_FamTotalPTimeOnLastMachine[FamPrmu[Fam]];
                break;
            case 6:
                ch[Fam].value = this->m_FamTotalPTimeOnLastMachine[FamPrmu[Fam]] + this->m_FamTotalPTimeOnFirstMachine[FamPrmu[Fam]];
                break;
            case 7:
                ch[Fam].value = this->m_FamTotalPTimeOnLastMachine[FamPrmu[Fam]] - this->m_FamTotalPTimeOnFirstMachine[FamPrmu[Fam]];
                break;

            default:cout << "Factor is out of range"; break;
        }
    }
    if (SortMethod == 0)
        sort(ch, ch + this->m_Families, PairGreater<double>());
    else
        sort(ch, ch + this->m_Families, PairLess<double>());
    for (int Fam = 0; Fam < this->m_Families; Fam++)
        FamPrmu[Fam] = ch[Fam].dim;
    delete[]ch;
}

/**
 * 对组根据某一属性排序
 * @param SortMethod
 * @param FamPermu
 */
void Heu::SortFam(int SortMethod, vector<int>& FamPermu)
{
    FamPermu.clear();
    FamPermu.resize(this->m_Families);
    for (int Fam = 0; Fam < this->m_Families; Fam++)
        FamPermu[Fam] = Fam;
    Pair<double>* ch = new Pair<double>[FamPermu.size()];

    if (SortMethod == 0) //LPT
    {
        for (int Fam = 0; Fam < FamPermu.size(); Fam++)
        {
            ch[Fam].dim = FamPermu[Fam];
            ch[Fam].value = this->m_FamTotalPTime[FamPermu[Fam]];
        }
        sort(ch, ch + this->m_Families, PairGreater<double>());
    }
    if (SortMethod == 1) //SPT
    {
        for (int Fam = 0; Fam < FamPermu.size(); Fam++)
        {
            ch[Fam].dim = FamPermu[Fam];
            ch[Fam].value = this->m_FamTotalPTime[FamPermu[Fam]];
        }
        sort(ch, ch + this->m_Families, PairLess<double>());
    }

    for (int Fam = 0; Fam < this->m_Families; Fam++)
        FamPermu[Fam] = ch[Fam].dim;
    delete[]ch;
}
/**
 * 对不完整的组序列进行排序
 * @param Factor
 * @param SortMethod
 * @param FamSeq
 */
void Heu::SortOtherFam(int Factor, int SortMethod, vector<int> &FamSeq)
{
    Pair<double> *ch = new Pair<double>[FamSeq.size()];
    for (int Fam = 0; Fam < FamSeq.size(); ++Fam)
    {
        ch[Fam].dim = FamSeq[Fam];
        switch (Factor)
        {
            case 0:
                ch[Fam].value = this->m_FamTotalPTime[FamSeq[Fam]];
                //根据总加工时间排序
                break;
            case 1:
                ch[Fam].value = this->m_FamTotalPTime[FamSeq[Fam]] + this->m_FamAvgSetupTime[FamSeq[Fam]];
                //根据总加工时间+组平均设置时间排序
                break;
            case 2:
                ch[Fam].value = this->m_FamTotalPTimeOnLastMachine[FamSeq[Fam]];
                //根据最后一台机器上的总加工时间排序
                break;
            default:
                cout << "Factor is out of range";
                break;
        }
    }
    if (SortMethod == 0)
        sort(ch, ch + FamSeq.size(), PairGreater<double>());  //根据属性降序排列
    else
        sort(ch, ch + FamSeq.size(), PairLess<double>());     //升序
    for (int Fam = 0; Fam < FamSeq.size(); Fam++)
        FamSeq[Fam] = ch[Fam].dim;   //将排序完的重新赋给FamSeq
    delete[]ch;

}


//Get a FamSeq in a factory by NEHinsert to Fam 20190420
int Heu::NEHFam(vector<int> FamPrmuInFac, vector<vector<int>> JobSeqInFam, vector<int>& FamSeqInFac, int& SpanInFac)
{
    vector<vector<int>> JFDTime(this->m_Jobs), JBDTime;
    for (int j = 0; j < this->m_Jobs; j++)
    {
        JFDTime[j].resize(this->m_Machines, 0);
    }
    JBDTime = JFDTime;
    for (int Fam = 0; Fam < FamPrmuInFac.size(); Fam++)
    {
        this->GetJFDTime_Forward_InFactory(FamSeqInFac, JobSeqInFam, JFDTime);
        this->GetJBDTime_Backward_InFactory(FamSeqInFac, JobSeqInFam, JBDTime);

        int CurFam = FamPrmuInFac[Fam];
        int bestFac = -1, bestPos = -1;
        SpanInFac = this->FindBestPosToInsertFam_InFactory(FamSeqInFac, JobSeqInFam, JFDTime, JBDTime ,CurFam, bestPos);

        //Insert CurFam to bestPos at bestFac
        FamSeqInFac.insert(FamSeqInFac.begin() + bestPos, CurFam);

    }

    return SpanInFac;
}

//xin重载
int Heu::NEHFam(vector<int> FamPrmu, vector<vector<int>> JobSeqinFam, vector<vector<int>>& FacFamSeq, vector<int>& FacSpan)
{
    // 检查 FacFamSeq and FacSpan
    if (FacFamSeq.size() != FacSpan.size())
    {
        cout << "Error in NEHJob: FacFamSeq.size() != FacSpan.size()" << endl;
        exit(0);
    }
    if (!FacFamSeq.size())
    {
        cout << "Error in NEHJob: FacFamSeq.size()=0" << endl;
        exit(0);
    }

    vector<vector<int>> JFDTime(this->m_Jobs), JBDTime;
    for (int j = 0; j < this->m_Jobs; j++)
    {
        JFDTime[j].resize(this->m_Machines, 0);
    }

    JBDTime = JFDTime;
    vector<bool> Is_Calculated(FacFamSeq.size(), false); //是否计算得到jCTime，jSTime在给定的工厂中；初始化false  whether jCTime and jSTime are calculated in the given Factory

    for (int Fam = 0; Fam < FamPrmu.size(); Fam++)
    {
        for (int Fac = 0; Fac < FacFamSeq.size(); Fac++)
        {
            if (!Is_Calculated[Fac])
            {
                this->GetJFDTime_Forward_InFactory(FacFamSeq[Fac], JobSeqinFam, JFDTime);  //前向计算
                this->GetJBDTime_Backward_InFactory(FacFamSeq[Fac], JobSeqinFam,  JBDTime); //后向计算
                Is_Calculated[Fac] = true; //设为true
            }
        }

        int CurFam = FamPrmu[Fam];  //当前组
        int bestFac = -1, bestPos = -1;  //最好的工厂，位置
        int Span = this->FindBestPosToInsertFam(FacFamSeq, JobSeqinFam, JBDTime, JFDTime, CurFam,  bestFac,bestPos);  //在所有工厂中找到最好的位置插入组得到完工时间

        // 在最好的工厂的最好位置中插入当前组 Insert CurFam to bestPos at bestFac
        FacFamSeq[bestFac].insert(FacFamSeq[bestFac].begin() + bestPos, CurFam);
        FacSpan[bestFac] = Span; //记录span
        Is_Calculated[bestFac] = false;  //最好工厂设为false
    }
    return *max_element(FacSpan.begin(), FacSpan.end());
}


void Heu::setDifference(vector<int> OriginalSequence, vector<int> RemovedSequence, vector<int>& RemainingSequence) {
    // 用于计算 OriginalSequence 和 RemovedSequence 的差集，结果存储在 RemainingSequence 中
    vector<int> DifferenceResult;  // 存储差集结果
    DifferenceResult.resize(max(OriginalSequence.size(), RemovedSequence.size()));  // 初始化 DifferenceResult 数组的大小

    // 对原始序列和移除序列进行排序
    sort(OriginalSequence.begin(), OriginalSequence.end());
    sort(RemovedSequence.begin(), RemovedSequence.end());

    // 计算差集，set_difference 会将差集存入 DifferenceResult 中
    vector<int>::iterator diffIterator = set_difference(OriginalSequence.begin(), OriginalSequence.end(),
                                                        RemovedSequence.begin(), RemovedSequence.end(),
                                                        DifferenceResult.begin());

    // 计算差集的大小
    int differenceLength = static_cast<int>(diffIterator - DifferenceResult.begin());

    // 清空 RemainingSequence，并将差集结果存入 RemainingSequence
    RemainingSequence.clear();
    for (int i = 0; i < differenceLength; i++) {
        RemainingSequence.push_back(DifferenceResult[i]);
    }
}

void Heu::JPA_TS(vector<vector<int>>& FactoryFamilySequence) {
    // 计算每两个组之间的 setup time 总和
    GetFamSumSetupTime();
    vector<vector<int>> SetupTimesBetweenFamilies = this->m_FamSumSetupTime;  // 存储组之间的 setup time

    // 初始化组序列，包含所有的组编号（0 到 m_Families-1）
    vector<int> AllFamilyIndices(this->m_Families);
    for (int i = 0; i < this->m_Families; i++) {
        AllFamilyIndices[i] = i;  // 将组编号存入 AllFamilyIndices 中
    }

    // 初始化工厂组访问记录，每个工厂的组访问序列
    vector<vector<int>> FactoryVisitedGroups;
    FactoryVisitedGroups.resize(m_Factories);  // 根据工厂数量初始化 FactoryVisitedGroups

    // 临时存储已分配的组，初始化为空
    vector<int> AllocatedFamilies;
    AllocatedFamilies.clear();

    // 为每个工厂分配一个初始组（按顺序分配）
    for (int factoryIndex = 0; factoryIndex < m_Factories; factoryIndex++) {
        FactoryVisitedGroups[factoryIndex].push_back(AllFamilyIndices[factoryIndex]);  // 每个工厂初始分配一个组
        AllocatedFamilies.push_back(AllFamilyIndices[factoryIndex]);  // 已分配的组加入 AllocatedFamilies
    }

    // 计算剩余未分配的组
    vector<int> RemainingFamilies;
    setDifference(AllFamilyIndices, AllocatedFamilies, RemainingFamilies);  // 使用 setDifference 来计算剩余组

    // 循环直到所有工厂都有完整的组序列
    while (true) {
        for (int factoryIndex = 0; factoryIndex < m_Factories; factoryIndex++) {
            // 获取当前工厂 FactoryVisitedGroups[factoryIndex] 中的最后一个组（尾部组）
            int LastAssignedFamily = FactoryVisitedGroups[factoryIndex].back();
            int MinimumSetupTime = INT_MAX;  // 最小切换设置时间初始化为最大值
            int BestFamilyToAssign = -1;  // 待分配的最优组初始化为 -1

            // 在剩余的组中选择与当前尾部组切换设置时间最小的组
            for (int remainingFamilyIndex = 0; remainingFamilyIndex < RemainingFamilies.size(); remainingFamilyIndex++) {
                int candidateFamily = RemainingFamilies[remainingFamilyIndex];
                int setupTime = SetupTimesBetweenFamilies[LastAssignedFamily][candidateFamily];

                if (setupTime < MinimumSetupTime) {
                    MinimumSetupTime = setupTime;  // 更新最小切换时间
                    BestFamilyToAssign = candidateFamily;  // 更新最优分配的组
                }
            }

            // 将找到的最优组分配给当前工厂
            FactoryVisitedGroups[factoryIndex].push_back(BestFamilyToAssign);
            AllocatedFamilies.push_back(BestFamilyToAssign);  // 将已分配的组加入 AllocatedFamilies

            // 更新剩余未分配的组
            setDifference(AllFamilyIndices, AllocatedFamilies, RemainingFamilies);

            // 如果所有组都已分配，则退出循环
            if (AllFamilyIndices.size() == AllocatedFamilies.size()) {
                break;
            }
        }

        // 如果所有组都已分配，退出主循环
        if (AllFamilyIndices.size() == AllocatedFamilies.size()) {
            break;
        }
    }

    // 将最终的工厂-组序列存储到 FactoryFamilySequence 中
    FactoryFamilySequence = FactoryVisitedGroups;
}

void Heu::JPA_G(vector<vector<int>>& FacFamSeq, vector<vector<int>>& InitialJobinFamSeq) {
    GetFamSumSetupTime();
    vector<vector<int>>D = this->m_FamSumSetupTime;  //计算两个组之间的setup times总和

    //JPA被执行产生famSeq
    vector<int>permutation(this->m_Families);
    for (int j = 0; j < this->m_Families; j++)
        permutation[j] = j;

    vector<vector<int>> Facvisited;
    Facvisited.resize(m_Factories);
    vector<int> temp;
    temp.clear();
    for (int fac = 0; fac < m_Factories; fac++)
    {
        Facvisited[fac].push_back(permutation[fac]);
        temp.push_back(permutation[fac]);
    }

    vector<int> Restseq;
    setDifference(permutation, temp, Restseq);

    while (true)
    {
        for (int fac = 0; fac < m_Factories; fac++)
        {
            int Tail = Facvisited[fac].back();
            int minDis = INT_MAX;
            int fam = -1;
            for (int p = 0; p < Restseq.size(); p++)
            {
                if (D[Tail][Restseq[p]] < minDis)
                {
                    minDis = D[Tail][Restseq[p]];
                    fam = Restseq[p];
                }
            }
            Facvisited[fac].push_back(fam);
            temp.push_back(fam);
            setDifference(permutation, temp, Restseq);
            if (permutation.size() == temp.size())
            {
                break;
            }

        }
        if (permutation.size() == temp.size())
            break;
    }

    vector<vector<int>> JobSeqinFam = m_JobsInEachFamily;
    LRX(Facvisited, JobSeqinFam);

    int Span;
    for (int fac = 0; fac < m_Factories; fac++)
    {
        FamInsert(Facvisited[fac], JobSeqinFam, Span);
    }

    FacFamSeq = Facvisited;
    InitialJobinFamSeq = JobSeqinFam;
}

void Heu::LRX(vector<vector<int>> FamPop, vector<vector<int>>& JobSeqinFamPop) {
    vector<int>famseq;
    vector<vector<int>>JobSeqinFam = JobSeqinFamPop;

    for (int fac = 0; fac < m_Factories; fac++)
    {
        vector<int> MachReadyTime(this->m_Machines, 0);
        famseq = FamPop[fac];

        for (int Fam = 0; Fam < famseq.size(); Fam++) {
            int CurFam = famseq[Fam];
            if (Fam == 0) {
                for (int m = 0; m < this->m_Machines; m++)
                    MachReadyTime[m] = this->m_SetupTime[m][CurFam][CurFam];
            }
            else //from the second group of jobs to the end;
            {
                int PreFam = famseq[Fam - 1];
                for (int m = 0; m < this->m_Machines; m++)
                    MachReadyTime[m] += this->m_SetupTime[m][PreFam][CurFam];
            }
            //cout << "1" << endl;
            vector<int>jobseq = JobSeqinFam[CurFam];
            if (jobseq.size() == 1)
            {
                int CurJob = jobseq.back();
                MachReadyTime[0] =  max(MachReadyTime[0] + this->m_TureJobOpertime[CurJob][0], MachReadyTime[1]);;//on the first machine

                //on the rest machine
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
                continue;
            }
            //cout << "2" << endl;
            schedulejob(jobseq, MachReadyTime);
            //cout << "3" << endl;
            JobSeqinFam[CurFam] = jobseq;
        }
    }
    JobSeqinFamPop = JobSeqinFam;
}

void Heu::schedulejob(vector<int>& jobseq, vector<int>& MachReadyTime) {

    vector<int>unschedulejob = jobseq;
    vector<int>schedulejob;
    schedulejob.clear();
    for (int j = 0; j < jobseq.size(); j++) {
        int minYPCLjk = INT_MAX;
        vector<int> bestvalue(this->m_Machines, 0);
        int Appendjob;
        for (int uj = 0; uj < unschedulejob.size(); uj++) {
            vector<int> TempMachReadyTime = MachReadyTime;
            int job = unschedulejob[uj];
           TempMachReadyTime[0]= max(TempMachReadyTime[0] + this->m_JobOperPTime[job][0], TempMachReadyTime[1]);

            //计算第一部分的和
            int ITij = 0;
            int IT = 0;
            for (int m = 1; m < this->m_Machines; m++)
            {
                if (jobseq.size() > 2) {
                    IT = (m_Machines * max(TempMachReadyTime[m - 1] - TempMachReadyTime[m], 0)) / (m + (j * (m_Machines - m)) / (jobseq.size()- 2));
                }
                else {

                    IT = 0;
                }
                ITij = ITij + IT;
                if (m == this->m_Machines - 1)
                {
                    TempMachReadyTime[m] =  TempMachReadyTime[m - 1] + this->m_JobOperPTime[job][m];
                }
                else
                {
                    TempMachReadyTime[m] = max( TempMachReadyTime[m - 1] + m_JobOperPTime[job][m],TempMachReadyTime[m + 1]);
                }
            }

            if (unschedulejob.size() == 1) {

                Appendjob = job;
                bestvalue = TempMachReadyTime;
                break;
            }
            //计算第二部分
            vector<int> RestJob;
            vector<int> RemoveJob;
            RemoveJob.clear();
            RemoveJob.push_back(job);
            setDifference(unschedulejob, RemoveJob, RestJob);
            vector<int> artjobtime(this->m_Machines, 0);
            for (int m = 0; m < this->m_Machines; m++) {
                for (int j = 0; j < RestJob.size(); j++)
                    artjobtime[m] = artjobtime[m] + m_JobOperPTime[RestJob[j]][m];
                artjobtime[m] = artjobtime[m] / RestJob.size();
            }

            vector<int> tempcopy = TempMachReadyTime;
            tempcopy[0] = max( tempcopy[0] + artjobtime[0], tempcopy[1]);
            for (int m = 1; m < this->m_Machines; m++)
            {
                if (m == this->m_Machines - 1)
                {
                    tempcopy[m]= tempcopy[m - 1] + this->m_JobOperPTime[job][m];
                }
                else
                {
                    tempcopy[m] = max( tempcopy[m - 1] + m_JobOperPTime[job][m], tempcopy[m + 1]);
                }
            }
            int ATjk = tempcopy[m_Machines - 1] + MachReadyTime[m_Machines - 1];
            int YPCLjk = ITij + ATjk;
            if (YPCLjk < minYPCLjk) {
                minYPCLjk = YPCLjk;
                Appendjob = job;
                bestvalue = TempMachReadyTime;
            }

        }
        MachReadyTime = bestvalue;
        schedulejob.push_back(Appendjob);
        setDifference(jobseq, schedulejob, unschedulejob);
    }
    jobseq = schedulejob;
}

