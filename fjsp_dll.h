#pragma once

#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#include<utility>
#include<vector>
#include<unordered_map>
#include<stdlib.h>
#include<time.h>
#include<memory>
#include<new>
#include<list>
#include<algorithm>
#include<queue>

using namespace std;
/*应当统一四个算法的使用流程，以GA为例：
GA ga();
ga.initialize();
ga.run_mode1();
ga.show_best();*/

//算例类
struct __declspec(dllexport) Instance
{
public:
	vector<vector<unordered_map<int, int>>> message;		//零件加工信息
	int work_num;											//工件数量
	int machine_num;										//机器数量
	int avg_ma;												//工序平均可加工机器

	//辅助甘特图输出
	vector<int> works;
	vector<double> start_times;
	vector<double> finish_times;
	vector<int> machines;

	Instance() :work_num(0), machine_num(0), avg_ma(0) {};
	Instance(const string& filename);
	void print_Instance();									//展示算例信息
};
//甘特节点类
struct __declspec(dllexport) Node
{
public:
	int work;
	int job;
	int machine;
	double start_time;
	double finish_time;
	std::weak_ptr<Node> jPre;
	std::weak_ptr<Node> mPre;
	std::weak_ptr<Node> jNext;
	std::weak_ptr<Node> mNext;

	Node() : work(0), job(0), machine(0), start_time(0), finish_time(0) {}
	Node(int work_n, int job_n, int machine_n) :
		work(work_n), job(job_n), machine(machine_n), start_time(0), finish_time(0) {}
	Node(int work_n, int job_n, int machine_n, double start, double finish) :
		work(work_n), job(job_n), machine(machine_n), start_time(start), finish_time(finish) {}
};
//随机数生成类
struct __declspec(dllexport) RandomNumber
{
public:
	RandomNumber();										//初始化时间种子
	double getRandomD(double a, double b);				//产生[a,b)间的随机浮点数
	int getRandomI(int a, int b);						//产生[a,b]间的随机整数
	vector<int> get_n_diff(int a, int b, int n);		//产生[a,b]间不同的n个整数
};

//解决方案类
class __declspec(dllexport) Solution
{
public:
	Solution() :makespan(0.0), fitness(0.0) {}
	Solution(vector<int>& job_s, vector<vector<int>>& machine_s);
	
	void set_makespan();
	void set_fitness();
	void set_job_gene(vector<int>& job_s) { job_gene = job_s; }
	void set_machine_gene(vector<vector<int>>& machine_s) { machine_gene = machine_s; }
	void set_instance(Instance i) { ins = i; }		//设置算例，仅初始化fjsp算法时调用一次
	void get_gantt(vector<int>& job_s, vector<vector<int>>& machine_s);

	void show_gantt();			//展示甘特图
	void show_job_gene();		//展示工序序列
	void show_machine_gene();	//展示机器序列
	void set_arg(vector<int>& job_s, vector<vector<int>>& machine_s);	//重新设置参数

	void greed_initialize(int stand_cycle, int job_num, int mac_num);	//贪婪方法处理初始序列,stand_cycle未解决方案完工时间停滞迭代次数,job_num与mac_num未每轮迭代产生的工序、机器邻域解数量
	vector<Solution> search_neighborhood(int job_num, int mac_num);		//产生job_num个工序邻域解和mac_num个机器邻域解
	Solution search_job_neigh(int work1, int job1, int work2, int job2);//产生一个工序邻域解,下标从1开始
	Solution search_mac_neigh(int work, int job, int machine);			//产生一个机器邻域解，下标从1开始

	const vector<int>& get_job_gene() const { return job_gene; }
	const vector<vector<int>>& get_machine_gene() const { return machine_gene; }
	double get_makespan()const { return makespan; }
	double get_fitness()const { return fitness; }

	static Instance ins;
	static RandomNumber r;
	vector<int> jobPos;							//工序突变位置记录
	vector<int> machinePos;						//机器突变位置记录
private:
	double makespan;							//最大完工时间
	double fitness;								//个体适应度
	vector<int> job_gene;						//工序基因
	vector<vector<int>> machine_gene;			//机器基因
	vector<list<shared_ptr<Node>>> gantt;		//甘特图
};
__declspec(dllexport) Instance Solution::ins = Instance();
__declspec(dllexport) RandomNumber Solution::r = RandomNumber();
/*对<,>和==进行重载，<与>用于比较两个解决方案的最大完工时间大小，==用于比较两个解决方案工序及机器序列是否相同*/
inline bool operator<(const Solution& left, const Solution& right)
{
	return left.get_makespan() < right.get_makespan();
}
inline bool operator>(const Solution& left, const Solution& right)
{
	return left.get_makespan() > right.get_makespan();
}
inline bool operator==(const Solution& left, const Solution& right)
{
	return (left.get_job_gene() == right.get_job_gene()) && (left.get_machine_gene() == right.get_machine_gene());
}

//遗传算法类
class __declspec(dllexport) GA
{
public:
	/*构造及算法运行函数*/
	GA() :PS(0), maxCycle(0), PB(0.0), PC(0.0), PM(0.0), PY(0.0), R(0) {}
	GA(Instance _Ins, int _PS, int _MaxCycle, double _PB, double _PC, double _PM, double _PY, int _R) :
		ins(_Ins), PS(_PS), maxCycle(_MaxCycle), PB(_PB), PC(_PC), PM(_PM), PY(_PY), R(_R) {}
	GA(Instance _Ins, int _MaxCycle);		//输入算例ins与最大迭代次数maxcycle，其他默认
	GA(string _Filename, int _MaxCycle);	//输入算例文件名filename与最大迭代次数maxcycle，其他默认

	void initialize_mode1();					//按默认方法随机初始化种群
	void initialize_mode2();					//按限制种群数量与迭代次数的方法随机初始化种群
	void initialize_mode3(int stand_cycle, int job_num, int mac_num);	//按贪婪方法随机初始化种群
	Solution run_mode1();						//迭代n次(n为默认次数)
	Solution run_mode1(int cycle);				//迭代cycle次
	Solution run_mode2(double time);			//迭代time秒
	double get_makespan() const { return best_Individual.get_makespan(); }//获取最优个体的makespan
	void show_best();							//展示最优解

	/*获取私有变量*/
	int get_PS() const { return PS; }
	int get_maxCycle() const { return maxCycle; }
	double get_PB() const { return PB; }
	double get_PC() const { return PC; }
	double get_PM() const { return PM; }
	int get_R() const { return R; }
	const vector<Solution>& get_population() const { return population; }
	const Solution& get_best_Individual() const { return best_Individual; }

private:
	/*不需要使用方调用，可以放入private的部分函数: */
	void searchBestSource();					//找到适应度最好的个体
	vector<Solution> get_best_n(int n);			//找到适应度最好的前n个个体
	vector<Solution> get_best_n(vector<Solution>& population, int n);
	vector<Solution> get_worst_n(int n);		//找到适应度最差的前n个个体
	void select();								//选择操作
	void cross(Solution& sol1, Solution& sol2);	//交叉操作
	void mutate(Solution& sol);					//变异操作
	void evolve();								//进化一代

	/*私有成员变量*/
	int PS;				//种群规模
	int maxCycle = 30;	//终止条件迭代次数
	double PB = 0.05;	//最优个体赋值到下一代的比例
	double PC = 0.75;	//交叉概率
	double PM = 0.1;	//变异概率
	double PY = 0.02;	//最好基因保留概率
	int R = 4;			//交叉对数

	Instance ins;
	RandomNumber r;

	vector<Solution> population;		//种群
	Solution best_Individual;			//适应度最好的个体
};

//粒子群算法类
class __declspec(dllexport) PSO
{
public:
	/*构造及算法运行函数*/
	PSO() {}
	PSO(double _PM, int _PN) :PM(_PM), PN(_PN) {}

	void initialize(Instance _ins);				//随机初始化种群
	void run_mode1(int cycle);					//按步数迭代
	void run_mode1();
	void run_mode2(int cycle);					//按最高点停滞步数迭代
	void run_mode2();
	void run_mode3(double time);				//按运行时间迭代
	void run_mode3();
	void show_best();							//展示最优解

	/*获取私有变量*/
	const Solution& get_best()const { return best1; }		//返回最优方案
	int get_ps()const { return PS; }
	double get_pm()const { return PM; }
	int get_pn()const { return PN; }
	int get_max_cycle()const { return maxCycle; }
	int get_max_stand_cycle()const { return maxStandCycle; }
	double get_max_run_time()const { return maxRunTime; }

private:
	/*可以不提供给使用方，放入private的部分: */
	void run_step(Solution& sol, int pn);		//随机移动一步方案sol
	void run_step(Solution& sol, Solution& best_one, int pn);//模仿最优方案移动一步方案sol
	void search_best();							//找到最优解

	/*私有成员变量*/
	int PS = 50;								//方案群数量
	double PM = 0.2;							//偏移概率
	int PN = 2;									//每次改变点位对数
	int maxCycle = 100;							//终止条件：代数
	int maxStandCycle = 25;						//终止条件：最优方案停滞代数
	double maxRunTime = 20.0;					//终止条件：程序运行时间
	vector<Solution> population;				//方案群
	Solution best1;								//最优方案
	Solution best2;								//次优方案
	Solution best3;								//次优方案
	Instance ins;
	RandomNumber r;
};

//模拟退火算法类
class __declspec(dllexport) SAA
{
public:
	/*构造及算法运行函数*/
	SAA() { NowT = StartT; }
	SAA(double _StartT, double _EndT, double _K, int _L) :StartT(_StartT), EndT(_EndT), K(_K), L(_L) { NowT = StartT; }

	void initialize(Instance ins);
	void initialize(Instance ins, double _K, int _L, double _k);
	void run();							//运行算法
	void show_best();					//展示最终解

	/*获取私有变量*/
	const Solution& getBest()const { return best; }
	double get_NowT()const { return NowT; }
	double get_EndT()const { return EndT; }
	double get_K()const { return K; }
	int get_L()const { return L; }
	double get_k()const { return k; }
	int get_MinStep()const { return MinStep; }

private:
	/*可以不提供给使用方，放入private的部分: */
	int stepLength();						//返回当前步长(变异点个数)
	void run_step(Solution& sol, int pn);	//随意移动一步方案sol，步长pn
	
	/*私有成员变量*/
	Solution best;					//储存最优解
	Solution s;						//当前解
	double NowT;					//当前温度
	double StartT = 10000;			//初始温度
	double EndT = 1;				//终止温度
	double K = 0.95;				//温度下降速率
	int L = 100;					//等温过程迭代次数（马尔可夫链长度）
	double k = 0.3;					//较劣方案接受概率
	int MinStep = 1;				//最低突变次数

	RandomNumber r;
	const double e = 2.7183;		//自然常数e
};

//禁忌搜索类
class __declspec(dllexport) TS
{
public:
	/*构造及算法运行函数*/
	TS() :tabuLength(10), job_num(10), mac_num(10) {}
	TS(int _tabuLength, int _job_num, int _mac_num, Instance ins) :
		tabuLength(_tabuLength), job_num(_job_num), mac_num(_mac_num) {
		//initialize(ins);
	}

	void initialize(Instance ins);
	void run_mode1(int maxCycle);		//按最大迭代次数为终止条件的方式进行禁忌搜索
	void run_mode1();
	void run_mode2(int maxStandCycle);	//按最优解停滞代数为终止条件进行禁忌搜索
	void run_mode2();
	void run_mode3(double maxRunTime);	//按运行时间为终止条件进行禁忌搜索
	void run_mode3();

	/*获取私有变量*/
	Solution get_best()const { return best; }
	Solution get_current()const { return current; }
	const vector<Solution>& get_tabuTable()const { return tabuTable; }
	int get_tabuLength()const { return tabuLength; }
	int get_job_num()const { return job_num; }
	int get_mac_num()const { return mac_num; }
	int get_max_cycle()const { return max_cycle; }
	int get_max_stand_cycle()const { return max_stand_cycle; }
	double get_max_run_time()const { return max_run_time; }

private:
	/*可以不提供给使用方，放入private的部分: */
	void show_cycle(int i);						//输出第i次迭代结果
	bool is_in_tabuTable(const Solution& sol);	//判断sol是否在禁忌表中

	/*私有成员变量*/
	vector<Solution> tabuTable;		//禁忌表
	Solution best;					//最优解
	Solution current;				//当前解
	int tabuLength;					//禁忌长度
	int job_num;		//单次生成工序邻域数量
	int mac_num;		//单次生成机器邻域数量

	int max_cycle = 100;				//（终止条件）最大迭代次数
	int max_stand_cycle = 20;			//（终止条件）最优方案停滞代数
	double max_run_time = 100.0;		//（终止条件）运行时间
	RandomNumber r;
};
