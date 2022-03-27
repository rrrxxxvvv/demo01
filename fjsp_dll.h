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
/*Ӧ��ͳһ�ĸ��㷨��ʹ�����̣���GAΪ����
GA ga();
ga.initialize();
ga.run_mode1();
ga.show_best();*/

//������
struct __declspec(dllexport) Instance
{
public:
	vector<vector<unordered_map<int, int>>> message;		//����ӹ���Ϣ
	int work_num;											//��������
	int machine_num;										//��������
	int avg_ma;												//����ƽ���ɼӹ�����

	//��������ͼ���
	vector<int> works;
	vector<double> start_times;
	vector<double> finish_times;
	vector<int> machines;

	Instance() :work_num(0), machine_num(0), avg_ma(0) {};
	Instance(const string& filename);
	void print_Instance();									//չʾ������Ϣ
};
//���ؽڵ���
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
//�����������
struct __declspec(dllexport) RandomNumber
{
public:
	RandomNumber();										//��ʼ��ʱ������
	double getRandomD(double a, double b);				//����[a,b)������������
	int getRandomI(int a, int b);						//����[a,b]����������
	vector<int> get_n_diff(int a, int b, int n);		//����[a,b]�䲻ͬ��n������
};

//���������
class __declspec(dllexport) Solution
{
public:
	Solution() :makespan(0.0), fitness(0.0) {}
	Solution(vector<int>& job_s, vector<vector<int>>& machine_s);
	
	void set_makespan();
	void set_fitness();
	void set_job_gene(vector<int>& job_s) { job_gene = job_s; }
	void set_machine_gene(vector<vector<int>>& machine_s) { machine_gene = machine_s; }
	void set_instance(Instance i) { ins = i; }		//��������������ʼ��fjsp�㷨ʱ����һ��
	void get_gantt(vector<int>& job_s, vector<vector<int>>& machine_s);

	void show_gantt();			//չʾ����ͼ
	void show_job_gene();		//չʾ��������
	void show_machine_gene();	//չʾ��������
	void set_arg(vector<int>& job_s, vector<vector<int>>& machine_s);	//�������ò���

	void greed_initialize(int stand_cycle, int job_num, int mac_num);	//̰�����������ʼ����,stand_cycleδ��������깤ʱ��ͣ�͵�������,job_num��mac_numδÿ�ֵ��������Ĺ��򡢻������������
	vector<Solution> search_neighborhood(int job_num, int mac_num);		//����job_num������������mac_num�����������
	Solution search_job_neigh(int work1, int job1, int work2, int job2);//����һ�����������,�±��1��ʼ
	Solution search_mac_neigh(int work, int job, int machine);			//����һ����������⣬�±��1��ʼ

	const vector<int>& get_job_gene() const { return job_gene; }
	const vector<vector<int>>& get_machine_gene() const { return machine_gene; }
	double get_makespan()const { return makespan; }
	double get_fitness()const { return fitness; }

	static Instance ins;
	static RandomNumber r;
	vector<int> jobPos;							//����ͻ��λ�ü�¼
	vector<int> machinePos;						//����ͻ��λ�ü�¼
private:
	double makespan;							//����깤ʱ��
	double fitness;								//������Ӧ��
	vector<int> job_gene;						//�������
	vector<vector<int>> machine_gene;			//��������
	vector<list<shared_ptr<Node>>> gantt;		//����ͼ
};
__declspec(dllexport) Instance Solution::ins = Instance();
__declspec(dllexport) RandomNumber Solution::r = RandomNumber();
/*��<,>��==�������أ�<��>���ڱȽ������������������깤ʱ���С��==���ڱȽ���������������򼰻��������Ƿ���ͬ*/
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

//�Ŵ��㷨��
class __declspec(dllexport) GA
{
public:
	/*���켰�㷨���к���*/
	GA() :PS(0), maxCycle(0), PB(0.0), PC(0.0), PM(0.0), PY(0.0), R(0) {}
	GA(Instance _Ins, int _PS, int _MaxCycle, double _PB, double _PC, double _PM, double _PY, int _R) :
		ins(_Ins), PS(_PS), maxCycle(_MaxCycle), PB(_PB), PC(_PC), PM(_PM), PY(_PY), R(_R) {}
	GA(Instance _Ins, int _MaxCycle);		//��������ins������������maxcycle������Ĭ��
	GA(string _Filename, int _MaxCycle);	//���������ļ���filename������������maxcycle������Ĭ��

	void initialize_mode1();					//��Ĭ�Ϸ��������ʼ����Ⱥ
	void initialize_mode2();					//��������Ⱥ��������������ķ��������ʼ����Ⱥ
	void initialize_mode3(int stand_cycle, int job_num, int mac_num);	//��̰�����������ʼ����Ⱥ
	Solution run_mode1();						//����n��(nΪĬ�ϴ���)
	Solution run_mode1(int cycle);				//����cycle��
	Solution run_mode2(double time);			//����time��
	double get_makespan() const { return best_Individual.get_makespan(); }//��ȡ���Ÿ����makespan
	void show_best();							//չʾ���Ž�

	/*��ȡ˽�б���*/
	int get_PS() const { return PS; }
	int get_maxCycle() const { return maxCycle; }
	double get_PB() const { return PB; }
	double get_PC() const { return PC; }
	double get_PM() const { return PM; }
	int get_R() const { return R; }
	const vector<Solution>& get_population() const { return population; }
	const Solution& get_best_Individual() const { return best_Individual; }

private:
	/*����Ҫʹ�÷����ã����Է���private�Ĳ��ֺ���: */
	void searchBestSource();					//�ҵ���Ӧ����õĸ���
	vector<Solution> get_best_n(int n);			//�ҵ���Ӧ����õ�ǰn������
	vector<Solution> get_best_n(vector<Solution>& population, int n);
	vector<Solution> get_worst_n(int n);		//�ҵ���Ӧ������ǰn������
	void select();								//ѡ�����
	void cross(Solution& sol1, Solution& sol2);	//�������
	void mutate(Solution& sol);					//�������
	void evolve();								//����һ��

	/*˽�г�Ա����*/
	int PS;				//��Ⱥ��ģ
	int maxCycle = 30;	//��ֹ������������
	double PB = 0.05;	//���Ÿ��帳ֵ����һ���ı���
	double PC = 0.75;	//�������
	double PM = 0.1;	//�������
	double PY = 0.02;	//��û���������
	int R = 4;			//�������

	Instance ins;
	RandomNumber r;

	vector<Solution> population;		//��Ⱥ
	Solution best_Individual;			//��Ӧ����õĸ���
};

//����Ⱥ�㷨��
class __declspec(dllexport) PSO
{
public:
	/*���켰�㷨���к���*/
	PSO() {}
	PSO(double _PM, int _PN) :PM(_PM), PN(_PN) {}

	void initialize(Instance _ins);				//�����ʼ����Ⱥ
	void run_mode1(int cycle);					//����������
	void run_mode1();
	void run_mode2(int cycle);					//����ߵ�ͣ�Ͳ�������
	void run_mode2();
	void run_mode3(double time);				//������ʱ�����
	void run_mode3();
	void show_best();							//չʾ���Ž�

	/*��ȡ˽�б���*/
	const Solution& get_best()const { return best1; }		//�������ŷ���
	int get_ps()const { return PS; }
	double get_pm()const { return PM; }
	int get_pn()const { return PN; }
	int get_max_cycle()const { return maxCycle; }
	int get_max_stand_cycle()const { return maxStandCycle; }
	double get_max_run_time()const { return maxRunTime; }

private:
	/*���Բ��ṩ��ʹ�÷�������private�Ĳ���: */
	void run_step(Solution& sol, int pn);		//����ƶ�һ������sol
	void run_step(Solution& sol, Solution& best_one, int pn);//ģ�����ŷ����ƶ�һ������sol
	void search_best();							//�ҵ����Ž�

	/*˽�г�Ա����*/
	int PS = 50;								//����Ⱥ����
	double PM = 0.2;							//ƫ�Ƹ���
	int PN = 2;									//ÿ�θı��λ����
	int maxCycle = 100;							//��ֹ����������
	int maxStandCycle = 25;						//��ֹ���������ŷ���ͣ�ʹ���
	double maxRunTime = 20.0;					//��ֹ��������������ʱ��
	vector<Solution> population;				//����Ⱥ
	Solution best1;								//���ŷ���
	Solution best2;								//���ŷ���
	Solution best3;								//���ŷ���
	Instance ins;
	RandomNumber r;
};

//ģ���˻��㷨��
class __declspec(dllexport) SAA
{
public:
	/*���켰�㷨���к���*/
	SAA() { NowT = StartT; }
	SAA(double _StartT, double _EndT, double _K, int _L) :StartT(_StartT), EndT(_EndT), K(_K), L(_L) { NowT = StartT; }

	void initialize(Instance ins);
	void initialize(Instance ins, double _K, int _L, double _k);
	void run();							//�����㷨
	void show_best();					//չʾ���ս�

	/*��ȡ˽�б���*/
	const Solution& getBest()const { return best; }
	double get_NowT()const { return NowT; }
	double get_EndT()const { return EndT; }
	double get_K()const { return K; }
	int get_L()const { return L; }
	double get_k()const { return k; }
	int get_MinStep()const { return MinStep; }

private:
	/*���Բ��ṩ��ʹ�÷�������private�Ĳ���: */
	int stepLength();						//���ص�ǰ����(��������)
	void run_step(Solution& sol, int pn);	//�����ƶ�һ������sol������pn
	
	/*˽�г�Ա����*/
	Solution best;					//�������Ž�
	Solution s;						//��ǰ��
	double NowT;					//��ǰ�¶�
	double StartT = 10000;			//��ʼ�¶�
	double EndT = 1;				//��ֹ�¶�
	double K = 0.95;				//�¶��½�����
	int L = 100;					//���¹��̵�������������ɷ������ȣ�
	double k = 0.3;					//���ӷ������ܸ���
	int MinStep = 1;				//���ͻ�����

	RandomNumber r;
	const double e = 2.7183;		//��Ȼ����e
};

//����������
class __declspec(dllexport) TS
{
public:
	/*���켰�㷨���к���*/
	TS() :tabuLength(10), job_num(10), mac_num(10) {}
	TS(int _tabuLength, int _job_num, int _mac_num, Instance ins) :
		tabuLength(_tabuLength), job_num(_job_num), mac_num(_mac_num) {
		//initialize(ins);
	}

	void initialize(Instance ins);
	void run_mode1(int maxCycle);		//������������Ϊ��ֹ�����ķ�ʽ���н�������
	void run_mode1();
	void run_mode2(int maxStandCycle);	//�����Ž�ͣ�ʹ���Ϊ��ֹ�������н�������
	void run_mode2();
	void run_mode3(double maxRunTime);	//������ʱ��Ϊ��ֹ�������н�������
	void run_mode3();

	/*��ȡ˽�б���*/
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
	/*���Բ��ṩ��ʹ�÷�������private�Ĳ���: */
	void show_cycle(int i);						//�����i�ε������
	bool is_in_tabuTable(const Solution& sol);	//�ж�sol�Ƿ��ڽ��ɱ���

	/*˽�г�Ա����*/
	vector<Solution> tabuTable;		//���ɱ�
	Solution best;					//���Ž�
	Solution current;				//��ǰ��
	int tabuLength;					//���ɳ���
	int job_num;		//�������ɹ�����������
	int mac_num;		//�������ɻ�����������

	int max_cycle = 100;				//����ֹ����������������
	int max_stand_cycle = 20;			//����ֹ���������ŷ���ͣ�ʹ���
	double max_run_time = 100.0;		//����ֹ����������ʱ��
	RandomNumber r;
};
