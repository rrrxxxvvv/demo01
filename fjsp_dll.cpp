#include "fjsp_dll.h"

Instance::Instance(const string & filename)
{
	ifstream ifs;
	ifs.open(filename, ios::in);
	string temp;
	getline(ifs, temp);
	istringstream line(temp);
	line >> work_num;
	message.resize(work_num + 1);

	line >> machine_num;
	line >> avg_ma;
	//gantt_message.resize(mNumber + 1);
	int work = 1;
	while (work <= work_num && ifs.eof() == 0)
	{
		getline(ifs, temp);
		istringstream line(temp);
		int job = 0;
		line >> job;
		message[work].resize(job + 1);
		for (int i = 1; i <= job; ++i)
		{
			int machine = 0;
			line >> machine;
			pair<int, int> p;
			while (machine--)
			{
				line >> p.first >> p.second;
				message[work][i].insert(p);//work号工件第i道工序，插入不同机器的加工时间
			}
		}
		++work;
	}
}

void Instance::print_Instance()
{
	cout << work_num << "\t" << machine_num << "\t" << avg_ma << endl;
	int i = 0;
	for (auto it = message.begin(); it != message.end(); it++)
	{
		cout << "工件" << ++i << ": " << endl;
		for (auto vit = (*it).begin(); vit != (*it).end(); vit++)
		{
			for (auto vvit = (*vit).begin(); vvit != (*vit).end(); vvit++)
				cout << "(" << vvit->first << "," << vvit->second << ")  ";
			cout << endl;
		}
	}
}

RandomNumber::RandomNumber()
{
	srand((unsigned)time(NULL));
}

double RandomNumber::getRandomD(double a, double b)
{
	return (rand() % 10000) / 10000.0 * (b - a) + a;
}

int RandomNumber::getRandomI(int a, int b)
{
	return rand() % (b - a + 1) + a;
}

vector<int> RandomNumber::get_n_diff(int a, int b, int n)
{
	vector<int> temp;
	for (int i = 0; i < n; ++i)
	{
		while (true)
		{
			int tmp = getRandomI(a, b);
			if (find(temp.begin(), temp.end(), tmp) == temp.end())
			{
				temp.push_back(tmp); break;
			}
		}
	}
	return temp;
}

Solution::Solution(vector<int>& job_s, vector<vector<int>>& machine_s)
{
	set_arg(job_s, machine_s);

	const int len = (int)job_s.size();
	jobPos = r.get_n_diff(0, len - 1, 2);
	machinePos = r.get_n_diff(1, len, 2);

}

void Solution::get_gantt(vector<int>& job_s, vector<vector<int>>& machine_s)
{
	auto mess = ins.message;
	//所有甘特图的第一个node设为NULL
	vector<vector<shared_ptr<Node>>> help(ins.work_num + 1, vector<shared_ptr<Node>>(1, shared_ptr<Node>(new Node())));	//记录每个工件加工到第几个工序；
	gantt.resize(ins.machine_num + 1);
	for (int i = 1; i <= ins.machine_num; ++i)gantt[i].clear();
	for (int i = 1; i <= ins.machine_num; ++i)gantt[i].push_back(shared_ptr<Node>(new Node()));

	for (int i = 0; i < (int)job_s.size(); ++i)
	{
		int work = job_s[i];							//工件号
		int job = help[work].size();					//该工件加工到第几个工序
		int machine = machine_s[work][job];				//该工件工序的加工机器
		double time = mess[work][job][machine];		//该工件工序的加工时间

		shared_ptr<Node> creat_node(new Node(work, job, machine));	//新工序节点
		double pre_finish_time = help[work].back()->finish_time;	//该工件上一道工序完成时间

		auto it_l = gantt[machine].begin();			//设置两个左右迭代器，判断之间是否可插入该工序
		auto it_r = it_l;
		++it_r;
		for (; it_r != gantt[machine].end(); ++it_r, ++it_l)
		{
			shared_ptr<Node> left = *it_l;
			shared_ptr<Node> right = *it_r;
			double start = left->finish_time;		//左边完成时间
			double finish = right->start_time;		//右边开始时间

			double beg = max(pre_finish_time, start);
			if (finish - beg >= time)
			{
				creat_node->start_time = beg;
				creat_node->finish_time = beg + time;
				(*it_l)->mNext = creat_node;
				creat_node->mPre = *it_l;
				creat_node->mNext = *it_r;
				(*it_r)->mPre = creat_node;
				break;
			}
		}
		if (it_r == gantt[machine].end())
		{
			creat_node->start_time = max((*it_l)->finish_time, pre_finish_time);
			creat_node->finish_time = creat_node->start_time + time;
			(*it_l)->mNext = creat_node;
			creat_node->mPre = *it_l;
		}
		gantt[machine].insert(it_r, creat_node);			//插入工序节点到甘特图

		shared_ptr<Node> temp = help[work].back();
		temp->jNext = creat_node;
		creat_node->jPre = temp;
		help[work].push_back(creat_node);					//插入工序节点到工件序列图

	}
}

void Solution::set_makespan()
{
	double t = 0.0;
	for (int i = 1; i <= ins.machine_num; ++i)
	{
		double temp = (gantt[i].back())->finish_time;
		if (t < temp)t = temp;
	}
	makespan = t;
}

void Solution::set_fitness()
{
	fitness = 1.0 / makespan;
}

void Solution::show_gantt()
{
	int j = 0;
	for (auto& i : gantt)
	{
		for (auto it = i.begin(); it != i.end(); ++it)
		{
			int work = (*it)->work;
			ins.works.push_back(work);
			ins.start_times.push_back((*it)->start_time);
			ins.finish_times.push_back((*it)->finish_time);
			ins.machines.push_back(j);
			double time = (*it)->finish_time - (*it)->start_time;
			while (time--)cout << work;
			shared_ptr<Node> shared_ptr_it = ((*it)->mNext).lock();
			if (shared_ptr_it)
			{
				double time = shared_ptr_it->start_time - (*it)->finish_time;
				while (time--)cout << " ";
			}
		}
		++j;
		cout << endl;
	}
}

void Solution::show_job_gene()
{
	for (int i : job_gene)cout << i;
	cout << endl;
}

void Solution::show_machine_gene()
{
	for (int i = 1; i < (int)machine_gene.size(); ++i)
	{
		for (int j = 1; j < (int)machine_gene[i].size(); ++j)
			cout << machine_gene[i][j];
		cout << " ";
	}
	cout << endl;
}

void Solution::set_arg(vector<int>& job_s, vector<vector<int>>& machine_s)
{
	job_gene = job_s;
	machine_gene = machine_s;
	get_gantt(job_s, machine_s);
	set_makespan();
	set_fitness();
}

void Solution::greed_initialize(int stand_cycle, int job_num, int mac_num)
{
	for (int i = 0; i < stand_cycle; i++)
	{
		Solution temp = *search_neighborhood(job_num, mac_num).begin();
		//cout << temp.get_makespan() << endl;
		if (temp.get_makespan() < this->get_makespan())
		{
			vector<int> job_s = temp.get_job_gene();
			vector<vector<int>> machine_s = temp.get_machine_gene();
			this->set_arg(job_s, machine_s);
			i = 0;
		}
	}
}

vector<Solution> Solution::search_neighborhood(int job_num, int mac_num)
{
	vector<Solution> neighborhood;
	int work = ins.work_num;
	//生成工序邻域解
	int work1, job1, work2, job2;//分别为交换的两个工序编号
	for (int i = 0; i < job_num; i++)
	{
		vector<int> tmp = r.get_n_diff(1, work, 2);
		work1 = tmp[0];
		work2 = tmp[1];
		job1 = r.getRandomI(1, machine_gene[work1].size() - 1);
		job2 = r.getRandomI(1, machine_gene[work2].size() - 1);
		Solution temp_sol = search_job_neigh(work1, job1, work2, job2);
		neighborhood.push_back(temp_sol);
	}
	//生成机器邻域解
	int work3, job3, machine3;//作为改变加工机器的工序编号和机器编号
	for (int i = 0; i < mac_num; i++)
	{
		work3 = r.getRandomI(1, work);
		job3 = r.getRandomI(1, machine_gene[work3].size() - 1);
		if (ins.message[work3][job3].size() == 1)
		{
			i--;
			continue;
		}
		vector<int> temp;//用来保存work3工件job3工序的可选加工机器号
		auto it = ins.message[work3][job3].begin();
		for (; it != ins.message[work3][job3].end(); it++)
		{
			if ((machine_gene[work3][job3]) != (*it).first)
				temp.push_back((*it).first);
		}
		machine3 = temp[r.getRandomI(0, temp.size() - 1)];
		Solution temp_sol = search_mac_neigh(work3, job3, machine3);
		neighborhood.push_back(temp_sol);
	}

	sort(neighborhood.begin(), neighborhood.end());			//此时neighborhood按makespan从低到高排列
	return neighborhood;
}

Solution Solution::search_job_neigh(int work1, int job1, int work2, int job2)
{
	vector<int> job_gene = (*this).job_gene;
	vector<vector<int>> machine_gene = (*this).machine_gene;

	int flag1 = 0, flag2 = 0;	//flag1和flag2分别记录遍历到两工件第几道工序
	for (int i = 0; i < job_gene.size(); i++)
	{
		if (job_gene[i] == work1)
			flag1++;
		if (job_gene[i] == work2)
			flag2++;
		if (job1 == flag1)
		{
			job_gene[i] = work2;
			flag1++;
		}
		if (job2 == flag2)
		{
			job_gene[i] = work1;
			flag2++;
		}
	}
	Solution jobNeigh(job_gene, machine_gene);
	return jobNeigh;
}

Solution Solution::search_mac_neigh(int work, int job, int machine)
{
	vector<int> job_gene = (*this).job_gene;
	vector<vector<int>> machine_gene = (*this).machine_gene;
	machine_gene[work][job] = machine;
	Solution macNeigh(job_gene, machine_gene);
	return macNeigh;
}

GA::GA(Instance _Ins, int _MaxCycle)
{
	this->ins = _Ins;
	PS = 10 * (ins.work_num + ins.machine_num);
	maxCycle = _MaxCycle;
}

GA::GA(string _Filename, int _MaxCycle)
{
	this->ins = Instance(_Filename);
	PS = 10 * (ins.work_num + ins.machine_num);
	maxCycle = _MaxCycle;
}

void GA::initialize_mode1()
{
	best_Individual.set_instance(ins);
	int work = ins.message.size();								//工件数
	vector<int> work_job(work);									//每个工件的工序数量
	for (int i = 1; i < work; ++i)
		work_job[i] = ins.message[i].size();

	vector<int> temp_job;										//产生工序序列
	for (int i = 1; i < work; ++i)
	{
		for (int j = 1; j < work_job[i]; ++j)
			temp_job.push_back(i);
	}

	vector<vector<int>> temp_machine(work);						//产生机器序列
	for (int i = 1; i < work; ++i)
		temp_machine[i].resize(work_job[i]);

	for (int i = 0; i < PS; ++i)								//产生PS个随机个体放入种群中
	{
		random_shuffle(temp_job.begin(), temp_job.end());
		for (int w = 1; w < work; ++w)
		{
			for (int j = 1; j < (int)ins.message[w].size(); ++j)
			{
				int machine = ins.message[w][j].size();
				vector<int> temp;
				auto b = ins.message[w][j].begin();
				auto e = ins.message[w][j].end();
				for (; b != e; ++b)temp.push_back(b->first);

				for (int m = 0; m < machine; ++m)
				{
					int ri = r.getRandomI(0, machine - 1);
					temp_machine[w][j] = temp[ri];
				}
			}
		}
		Solution ind(temp_job, temp_machine);
		population.push_back(ind);
	}
	searchBestSource();
}

void GA::initialize_mode2()
{
	best_Individual.set_instance(ins);
	PS = 6 * (ins.work_num + ins.machine_num);
	while (PS > 100)
		PS = 100;
	maxCycle = ins.avg_ma * (ins.work_num + ins.machine_num);
	while (maxCycle > 40)
		maxCycle = 40;

	int work = ins.message.size();								//工件数
	vector<int> work_job(work);									//每个工件的工序数量
	for (int i = 1; i < work; ++i)
		work_job[i] = ins.message[i].size();

	vector<int> temp_job;										//产生工序序列
	for (int i = 1; i < work; ++i)
	{
		for (int j = 1; j < work_job[i]; ++j)
			temp_job.push_back(i);
	}

	vector<vector<int>> temp_machine(work);						//产生机器序列
	for (int i = 1; i < work; ++i)
		temp_machine[i].resize(work_job[i]);

	for (int i = 0; i < PS; ++i)								//产生PS个随机个体放入种群中
	{
		random_shuffle(temp_job.begin(), temp_job.end());
		for (int w = 1; w < work; ++w)
		{
			for (int j = 1; j < (int)ins.message[w].size(); ++j)
			{
				int machine = ins.message[w][j].size();
				vector<int> temp;
				auto b = ins.message[w][j].begin();
				auto e = ins.message[w][j].end();
				for (; b != e; ++b)temp.push_back(b->first);

				for (int m = 0; m < machine; ++m)
				{
					int ri = r.getRandomI(0, machine - 1);
					temp_machine[w][j] = temp[ri];
				}
			}
		}
		Solution ind(temp_job, temp_machine);
		population.push_back(ind);
	}
	searchBestSource();
}

void GA::initialize_mode3(int stand_cycle, int job_num, int mac_num)
{
	best_Individual.set_instance(ins);
	PS = (ins.work_num + ins.machine_num) / 4 * 2;

	best_Individual.set_instance(ins);
	int work = ins.message.size();									//工件数
	vector<int> work_job(work);									//每个工件的工序数量
	for (int i = 1; i < work; ++i)
		work_job[i] = ins.message[i].size();

	vector<int> temp_job;										//产生工序序列
	for (int i = 1; i < work; ++i)
	{
		for (int j = 1; j < work_job[i]; ++j)
			temp_job.push_back(i);
	}

	vector<vector<int>> temp_machine(work);						//初始化机器序列容器尺寸
	for (int i = 1; i < work; ++i)
		temp_machine[i].resize(work_job[i]);

	cout << "开始贪婪方法的种群初始化:" << endl;
	for (int i = 0; i < PS; ++i)								//产生PS个随机个体放入种群中
	{
		random_shuffle(temp_job.begin(), temp_job.end());
		for (int w = 1; w < work; ++w)
		{
			for (int j = 1; j < (int)ins.message[w].size(); ++j)
			{
				int machine = ins.message[w][j].size();
				vector<int> temp;
				auto b = ins.message[w][j].begin();
				auto e = ins.message[w][j].end();
				for (; b != e; ++b)temp.push_back(b->first);

				for (int m = 0; m < machine; ++m)
				{
					int ri = r.getRandomI(0, machine - 1);
					temp_machine[w][j] = temp[ri];
				}
			}
		}
		Solution ind(temp_job, temp_machine);
		//cout << i + 1 << "初始个体: " << ind.get_makespan();
		ind.greed_initialize(stand_cycle, job_num, mac_num);
		population.push_back(ind);
		cout << i + 1 << ": " << ind.get_makespan() << "  ";
		//ind.show_job_gene();
		//ind.show_machine_gene();
	}
	searchBestSource();
}

void GA::searchBestSource()
{
	best_Individual = *min_element(population.begin(), population.end());
}

vector<Solution> GA::get_best_n(int n)
{
	partial_sort(population.begin(), population.begin() + n, population.end());
	return vector<Solution>(population.begin(), population.begin() + n);
}

vector<Solution> GA::get_best_n(vector<Solution>& population, int n)
{
	partial_sort(population.begin(), population.begin() + n, population.end());
	return vector<Solution>(population.begin(), population.begin() + n);
}

vector<Solution> GA::get_worst_n(int n)
{
	partial_sort(population.begin(), population.begin() + n, population.end(), std::greater<Solution>());
	return vector<Solution>(population.begin(), population.begin() + n);
}

void GA::select()
{
	vector<Solution> population_1(PS);
	int best_num = static_cast<int>(PS * PB);
	vector<Solution> best_n = get_best_n(best_num);
	for (int i = 0; i < best_num; ++i)
		population_1[i] = best_n[i];	//最优个体直接复制

	double total_fitness = 0.0;
	for (auto& it : population)
		total_fitness += it.get_fitness();
	vector<double> p(PS);
	for (int i = 0; i < PS; ++i)
		p[i] = population[i].get_fitness() / total_fitness;		//每个个体适应度占种群的比例
	for (int i = best_num; i < PS; ++i)
	{
		double r0_1 = r.getRandomD(0, 1);
		double cur = 0.0;
		for (int j = 0; j < PS; ++j)
		{
			cur += p[j];
			if (r0_1 < cur)
			{
				population_1[i] = population[j];				//轮盘赌选择剩余个体
				break;
			}
		}
	}
	population = population_1;
}

void GA::cross(Solution & sol_1, Solution & sol_2)
{
	//获取需要交叉的两个个体的工序基因和机器基因
	vector<int> job_gene_1 = sol_1.get_job_gene();
	vector<vector<int>> machine_gene_1 = sol_1.get_machine_gene();
	vector<int> job_gene_2 = sol_2.get_job_gene();
	vector<vector<int>> machine_gene_2 = sol_2.get_machine_gene();

	//现在进行工序基因交叉
	const int len = (int)job_gene_1.size();
	//随机从n个工件中选择R个进行交叉
	vector<int> choose_job = r.get_n_diff(1, ins.work_num, R);
	//记录工序基因中需要交叉的基因的位置
	vector<bool> choose_gene_1(len, 0);
	vector<bool> choose_gene_2(len, 0);
	for (int j = 0; j < len; ++j)
	{
		if (find(choose_job.begin(), choose_job.end(), job_gene_1[j]) != choose_job.end())
		{
			choose_gene_1[j] = 1;
		}
		if (find(choose_job.begin(), choose_job.end(), job_gene_2[j]) != choose_job.end())
		{
			choose_gene_2[j] = 1;
		}
	}
	//交叉
	queue<int> remain_gene_1;
	queue<int> remain_gene_2;
	for (int i = 0; i < len; ++i)
	{
		if (choose_gene_1[i] == 0)remain_gene_1.push(job_gene_1[i]);
		if (choose_gene_2[i] == 1)remain_gene_2.push(job_gene_2[i]);
	}
	for (int i = 0; i < len; ++i)
	{
		if (choose_gene_1[i] == 1)
		{
			job_gene_1[i] = remain_gene_2.front();
			remain_gene_2.pop();
		}
		if (choose_gene_2[i] == 0)
		{
			job_gene_2[i] = remain_gene_1.front();
			remain_gene_1.pop();
		}
	}

	//现在进行机器基因交叉
	const int work = (int)machine_gene_1.size();

	//随机产生0或1的机器数组
	vector<vector<int>> r01 = machine_gene_1;
	for (int i = 1; i < work; ++i)
	{
		for (int j = 1; j < (int)machine_gene_1[i].size(); ++j)
			r01[i][j] = r.getRandomI(0, 1);
	}
	//机器数组随机点交换
	for (int i = 1; i < work; ++i)
	{
		for (int j = 1; j < (int)machine_gene_1[i].size(); ++j)
		{
			if (r01[i][j] == 1)
				swap(machine_gene_1[i][j], machine_gene_2[i][j]);
		}
	}

	//重新设置两个个体的各项参数
	sol_1.set_arg(job_gene_1, machine_gene_1);
	sol_2.set_arg(job_gene_2, machine_gene_2);
}

void GA::mutate(Solution & sol)
{
	vector<int> job_gene = sol.get_job_gene();
	vector<vector<int>> machine_gene = sol.get_machine_gene();

	//对工序基因变异
	const int len = (int)job_gene.size();
	vector<int> pos = r.get_n_diff(0, len - 1, 2);
	int left = pos[0], right = pos[1];
	//double rate = r.getRandomD(0, 1);
	double rate = 0.66;
	if (rate > 0.5)
	{
		int num = job_gene[left];
		if (left > right)
		{
			job_gene.erase(job_gene.begin() + left);
			job_gene.insert(job_gene.begin() + right, num);
		}
		else
		{
			job_gene.insert(job_gene.begin() + right, num);
			job_gene.erase(job_gene.begin() + left);
		}
	}
	else
	{
		swap(job_gene[left], job_gene[right]);
	}

	//对机器基因变异
	//随机选择两个机器基因的位置
	vector<int> p = r.get_n_diff(1, len, 2);
	int p1 = p[0], p2 = p[1];
	int work1 = 1, work2 = 1, job1 = 1, job2 = 1;
	int count1 = 0, count2 = 0;
	while (work1 < ins.work_num)
	{
		count1 += (int)machine_gene[work1].size() - 1;
		if (count1 >= p1)
		{
			job1 = machine_gene[work1].size() - 1 + p1 - count1;
			break;
		}
		++work1;
	}
	while (work2 < ins.work_num)
	{
		count2 += (int)machine_gene[work2].size() - 1;
		if (count2 >= p2)
		{
			job2 = machine_gene[work2].size() - 1 + p2 - count2;
			break;
		}
		++work2;
	}

	int machine1 = ins.message[work1][job1].size();
	vector<int> temp1;
	int min = 70000;
	int ri1 = 1;
	auto b1 = ins.message[work1][job1].begin();
	auto e1 = ins.message[work1][job1].end();
	for (; b1 != e1; ++b1)temp1.push_back(b1->first);
	for (auto it = ins.message[work1][job1].begin(); it != e1; ++it) {

		if (min > it->second) {
			min = it->second;
			ri1 = it->first;
		}
	}
	machine_gene[work1][job1] = ri1;

	int machine2 = ins.message[work2][job2].size();
	vector<int> temp2;
	auto b2 = ins.message[work2][job2].begin();
	auto e2 = ins.message[work2][job2].end();
	for (; b2 != e2; ++b2)temp2.push_back(b2->first);
	int ri2 = r.getRandomI(0, machine2 - 1);
	machine_gene[work2][job2] = temp2[ri2];

	sol.set_arg(job_gene, machine_gene);
}

void GA::evolve()
{
	int n = static_cast<int>(PS * PY);
	vector<Solution> best_n = get_best_n(n);		//获取种群中最好的n个解
	select();
	vector<Solution> new_population;
	for (int i = 1; i <= PS / 2; ++i)
	{
		vector<int> get_2 = r.get_n_diff(0, PS - 1, 2);
		Solution ind_1 = population[get_2[0]];
		Solution ind_2 = population[get_2[1]];
		double rate_c = r.getRandomD(0, 1);
		if (rate_c < PC)					//如果小于交叉概率，进行交叉
		{
			cross(ind_1, ind_2);	//交叉
		}
		double rate_m_1 = r.getRandomD(0, 1);
		double rate_m_2 = r.getRandomD(0, 1);
		if (rate_m_1 < PM)
		{
			mutate(ind_1);					//如果小于变异概率，进行变异
		}
		if (rate_m_2 < PM)
		{
			mutate(ind_2);
		}

		new_population.push_back(ind_1);
		new_population.push_back(ind_2);
	}
	population = new_population;
	vector<Solution> new_best_n = get_best_n(n);
	vector<Solution> temp(2 * n);
	merge(best_n.begin(), best_n.end(), new_best_n.begin(), new_best_n.end(), temp.begin());
	best_n = get_best_n(temp, n);

	partial_sort(population.begin(), population.begin() + n, population.end(), greater<Solution>());
	for (int i = 0; i < n; ++i)
	{
		population[i] = best_n[i];
	}
	searchBestSource();
}

Solution GA::run_mode1()
{
	Solution best = best_Individual;
	cout << best.get_makespan() << endl;
	for (int i = 1; i < maxCycle; ++i)
	{
		evolve();
		double best_time = best_Individual.get_makespan();
		cout << best_time << endl;
		if (best_Individual < best)
		{
			best = best_Individual;
			cout << "迭代次数:" << i << " 当前最优解:" << best.get_makespan() << endl;
		}
	}
	return best;
}

Solution GA::run_mode1(int cycle)
{
	Solution best = best_Individual;
	cout << best.get_makespan() << endl;
	for (int i = 1; i < cycle; ++i)
	{
		evolve();
		double best_time = best_Individual.get_makespan();
		cout << best_time << endl;
		if (best_Individual < best)
		{
			best = best_Individual;
			cout << "迭代次数:" << i << " 当前最优解:" << best.get_makespan() << endl;
		}
	}
	return best;
}

Solution GA::run_mode2(double time)
{
	Solution best = best_Individual;
	cout << best.get_makespan() << endl;
	double start_time, current_time;
	start_time = current_time = clock();
	int i = 1;
	while ((double)(current_time - start_time) / CLOCKS_PER_SEC < time)
	{
		evolve();
		double best_time = best_Individual.get_makespan();
		cout << best_time << endl;
		if (best_Individual < best)
		{
			best = best_Individual;
			cout << "迭代次数:" << i++ << " 当前最优解:" << best.get_makespan() << endl;
		}
		current_time = clock();
	}
	return best;
}

void GA::show_best()
{
	cout << "最终解: " << best_Individual.get_makespan() << "\t最终解甘特图如下: ";
	best_Individual.show_gantt();
	cout << "最终解工序基因: " << endl;
	best_Individual.show_job_gene();
	cout << "最终解机器基因" << endl;
	best_Individual.show_machine_gene();
}

void PSO::initialize(Instance _ins)
{
	ins = _ins;
	best1.set_instance(ins);
	PS = 3.5 * (ins.work_num + ins.machine_num);
	maxCycle = 2 * (ins.work_num + ins.machine_num);
	maxStandCycle = ins.work_num + ins.machine_num;
	maxRunTime = ins.work_num + ins.machine_num;
	if (PS > 100)
		PS = 100;
	if (maxCycle > 80)
		maxCycle = 80;
	if (maxStandCycle > 25)
		maxStandCycle = 25;

	int work = ins.message.size();								//工件数
	vector<int> work_job(work);									//每个工件的工序数量
	for (int i = 1; i < work; ++i)
		work_job[i] = ins.message[i].size();

	vector<int> temp_job;										//产生工序序列
	for (int i = 1; i < work; ++i)
	{
		for (int j = 1; j < work_job[i]; ++j)
			temp_job.push_back(i);
	}

	vector<vector<int>> temp_machine(work);						//产生机器序列
	for (int i = 1; i < work; ++i)
		temp_machine[i].resize(work_job[i]);

	for (int i = 0; i < PS; ++i)								//产生PS个随机个体放入种群中
	{
		random_shuffle(temp_job.begin(), temp_job.end());
		for (int w = 1; w < work; ++w)
		{
			for (int j = 1; j < (int)ins.message[w].size(); ++j)
			{
				int machine = ins.message[w][j].size();
				vector<int> temp;
				auto b = ins.message[w][j].begin();
				auto e = ins.message[w][j].end();
				for (; b != e; ++b)temp.push_back(b->first);

				for (int m = 0; m < machine; ++m)
				{
					int ri = r.getRandomI(0, machine - 1);
					temp_machine[w][j] = temp[ri];
				}
			}
		}
		Solution ind(temp_job, temp_machine);
		population.push_back(ind);
	}
	sort(population.begin(), population.end());	//将population按makespan从小到大排序
	best1 = population[0];
	best2 = population[1];
	best3 = population[2];
}

void PSO::run_step(Solution & sol, int pn)
{
	vector<int> job_ge = sol.get_job_gene();
	vector<vector<int>>  machine_ge = sol.get_machine_gene();
	for (int i = 0; i < pn; i++)
	{
		//对工序基因变异
		const int len = (int)job_ge.size();
		vector<int> pos = r.get_n_diff(0, len - 1, 2);
		sol.jobPos = pos;
		int left = pos[0], right = pos[1];
		double rate = r.getRandomD(0, 1);
		if (rate > 0.5)
		{
			int num = job_ge[left];
			if (left > right)
			{
				job_ge.erase(job_ge.begin() + left);
				job_ge.insert(job_ge.begin() + right, num);
			}
			else
			{
				job_ge.insert(job_ge.begin() + right, num);
				job_ge.erase(job_ge.begin() + left);
			}
		}
		else
		{
			swap(job_ge[left], job_ge[right]);
		}

		//对机器基因变异
		vector<int> p = r.get_n_diff(1, len, 2);
		sol.machinePos = p;
		int p1 = p[0], p2 = p[1];
		int work1 = 1, work2 = 1, job1 = 1, job2 = 1;
		int count1 = 0, count2 = 0;
		while (work1 < ins.work_num)
		{
			count1 += (int)sol.get_machine_gene()[work1].size() - 1;
			if (count1 >= p1)
			{
				job1 = sol.get_machine_gene()[work1].size() - 1 + p1 - count1;
				break;
			}
			++work1;
		}
		while (work2 < ins.work_num)
		{
			count2 += (int)machine_ge[work2].size() - 1;
			if (count2 >= p2)
			{
				job2 = machine_ge[work2].size() - 1 + p2 - count2;
				break;
			}
			++work2;
		}

		int machine1 = ins.message[work1][job1].size();
		vector<int> temp1;
		int min = 70000;
		int ri1 = 1;
		auto b1 = ins.message[work1][job1].begin();
		auto e1 = ins.message[work1][job1].end();
		for (; b1 != e1; ++b1)temp1.push_back(b1->first);
		for (auto it = ins.message[work1][job1].begin(); it != e1; ++it) {

			if (min > it->second) {
				min = it->second;
				ri1 = it->first;
			}
		}
		machine_ge[work1][job1] = ri1;

		int machine2 = ins.message[work2][job2].size();
		vector<int> temp2;
		auto b2 = ins.message[work2][job2].begin();
		auto e2 = ins.message[work2][job2].end();
		for (; b2 != e2; ++b2)temp2.push_back(b2->first);
		int ri2 = r.getRandomI(0, machine2 - 1);
		machine_ge[work2][job2] = temp2[ri2];
	}

	sol.set_arg(job_ge, machine_ge);

}

void PSO::run_step(Solution & sol, Solution & best_one, int pn)
{
	vector<int> job_ge = sol.get_job_gene();
	vector<vector<int>> machine_ge = sol.get_machine_gene();

	if (r.getRandomD(0, 1) < 0.6)
	{
		//对工序基因变异
		const int len = (int)job_ge.size();
		vector<int> pos = best_one.jobPos;
		int left = pos[0], right = pos[1];
		double rate = r.getRandomD(0, 1);
		if (rate > 0.5)
		{
			int num = job_ge[left];
			if (left > right)
			{
				job_ge.erase(job_ge.begin() + left);
				job_ge.insert(job_ge.begin() + right, num);
			}
			else
			{
				job_ge.insert(job_ge.begin() + right, num);
				job_ge.erase(job_ge.begin() + left);
			}
		}
		else
		{
			swap(job_ge[left], job_ge[right]);
		}

		//对机器基因变异
		vector<int> p = best_one.machinePos;
		int p1 = p[0], p2 = p[1];
		int work1 = 1, work2 = 1, job1 = 1, job2 = 1;
		int count1 = 0, count2 = 0;
		while (work1 < ins.work_num)
		{
			count1 += (int)sol.get_machine_gene()[work1].size() - 1;
			if (count1 >= p1)
			{
				job1 = sol.get_machine_gene()[work1].size() - 1 + p1 - count1;
				break;
			}
			++work1;
		}
		while (work2 < ins.work_num)
		{
			count2 += (int)machine_ge[work2].size() - 1;
			if (count2 >= p2)
			{
				job2 = machine_ge[work2].size() - 1 + p2 - count2;
				break;
			}
			++work2;
		}

		int machine1 = ins.message[work1][job1].size();
		vector<int> temp1;
		int min = 70000;
		int ri1 = 1;
		auto b1 = ins.message[work1][job1].begin();
		auto e1 = ins.message[work1][job1].end();
		for (; b1 != e1; ++b1)temp1.push_back(b1->first);
		for (auto it = ins.message[work1][job1].begin(); it != e1; ++it) {

			if (min > it->second) {
				min = it->second;
				ri1 = it->first;
			}
		}
		machine_ge[work1][job1] = ri1;

		int machine2 = ins.message[work2][job2].size();
		vector<int> temp2;
		auto b2 = ins.message[work2][job2].begin();
		auto e2 = ins.message[work2][job2].end();
		for (; b2 != e2; ++b2)temp2.push_back(b2->first);
		int ri2 = r.getRandomI(0, machine2 - 1);
		machine_ge[work2][job2] = temp2[ri2];
	}
	else
	{
		//获取需要交叉的两个个体的工序基因和机器基因
		vector<int>& job_gene_1 = job_ge;
		vector<vector<int>>& machine_gene_1 = machine_ge;
		vector<int> job_gene_2 = best_one.get_job_gene();
		vector<vector<int>> machine_gene_2 = best_one.get_machine_gene();

		//现在进行工序基因交叉
		const int len = (int)job_gene_1.size();
		//随机从n个工件中选择R个进行交叉
		vector<int> choose_job = r.get_n_diff(1, ins.work_num, pn);
		//记录工序基因中需要交叉的基因的位置
		vector<bool> choose_gene_1(len, 0);
		vector<bool> choose_gene_2(len, 0);
		for (int j = 0; j < len; ++j)
		{
			if (find(choose_job.begin(), choose_job.end(), job_gene_1[j]) != choose_job.end())
			{
				choose_gene_1[j] = 1;
			}
			if (find(choose_job.begin(), choose_job.end(), job_gene_2[j]) != choose_job.end())
			{
				choose_gene_2[j] = 1;
			}
		}
		//交叉
		queue<int> remain_gene_1;
		queue<int> remain_gene_2;
		for (int i = 0; i < len; ++i)
		{
			if (choose_gene_1[i] == 0)remain_gene_1.push(job_gene_1[i]);
			if (choose_gene_2[i] == 1)remain_gene_2.push(job_gene_2[i]);
		}
		for (int i = 0; i < len; ++i)
		{
			if (choose_gene_1[i] == 1)
			{
				job_gene_1[i] = remain_gene_2.front();
				remain_gene_2.pop();
			}
			if (choose_gene_2[i] == 0)
			{
				job_gene_2[i] = remain_gene_1.front();
				remain_gene_1.pop();
			}
		}

		//现在进行机器基因交叉
		const int work = (int)machine_gene_1.size();

		//随机产生0或1的机器数组
		vector<vector<int>> r01 = machine_gene_1;
		for (int i = 1; i < work; ++i)
		{
			for (int j = 1; j < (int)machine_gene_1[i].size(); ++j)
				r01[i][j] = r.getRandomI(0, 1);
		}
		//机器数组随机点交换
		for (int i = 1; i < work; ++i)
		{
			for (int j = 1; j < (int)machine_gene_1[i].size(); ++j)
			{
				if (r01[i][j] == 1)
					swap(machine_gene_1[i][j], machine_gene_2[i][j]);
			}
		}
	}

	sol.set_arg(job_ge, machine_ge);
}

void PSO::search_best()
{
	sort(population.begin(), population.end());//按照makespan从小到大排列
	//reverse(population.begin(), population.end());

	vector<Solution>::iterator it = population.begin();
	vector<Solution> vtmp;
	vtmp.push_back(best1);
	vtmp.push_back(best2);
	vtmp.push_back(best3);
	vtmp.push_back(*it);
	vtmp.push_back(*(++it));
	vtmp.push_back(*(++it));
	sort(vtmp.begin(), vtmp.end());//按照makespan从小到大排列
	//reverse(vtmp.begin(), vtmp.end());

	/*it = vtmp.begin();
	best1 = (*it);
	best2 = (*(++it));
	best3 = (*(++it));*/
	best1 = vtmp[0];
	best2 = vtmp[1];
	best3 = vtmp[2];
}

void PSO::show_best()
{
	cout << "最终解: " << best1.get_makespan() << "\t最终解甘特图如下: ";
	best1.show_gantt();
	cout << "最终解工序基因: " << endl;
	best1.show_job_gene();
	cout << "最终解机器基因" << endl;
	best1.show_machine_gene();
}

void PSO::run_mode1(int cycle)
{
	search_best();
	for (int i = 0; i < cycle; i++)
	{
		vector<Solution>::iterator it = population.begin();
		for (; it != population.end(); it++)
		{
			double tmp = r.getRandomD(0, 1);
			if (tmp < PM)
				run_step(*it, PN);
			else
			{
				double rn = r.getRandomD(0, 1);
				if (rn < 0.5)
					run_step(*it,best1,PN);
				else if (rn > 0.8)
					run_step(*it, best3, PN);
				else
					run_step(*it, best2, PN);
			}
		}
		search_best();
		cout << " " << i + 1 << " " << best1.get_makespan() << endl;
	}
}

void PSO::run_mode1()
{
	search_best();
	for (int i = 0; i < maxCycle; i++)
	{
		vector<Solution>::iterator it = population.begin();
		for (; it != population.end(); it++)
		{
			double tmp = r.getRandomD(0, 1);
			if (tmp < PM)
				run_step(*it, PN);
			else
			{
				double rn = r.getRandomD(0, 1);
				if (rn < 0.5)
					run_step(*it, best1, PN);
				else if (rn > 0.8)
					run_step(*it, best3, PN);
				else
					run_step(*it, best2, PN);
			}
		}
		search_best();
		cout << " " << i + 1 << " " << best1.get_makespan() << endl;
	}
}

void PSO::run_mode2(int cycle)
{
	int i = 1, j = 0;
	bool is_true = true;
	search_best();
	Solution tmp = best1;
	while (is_true)
	{
		vector<Solution>::iterator it = population.begin();
		for (; it != population.end(); it++)
		{
			double tmp = 1.0 * rand() / RAND_MAX * 1.0;
			if (tmp < PM)
				run_step(*it, PN);
			else
			{
				double rn = r.getRandomD(0, 1);
				if (rn < 0.5)
					run_step(*it, best1, PN);
				else if (rn > 0.8)
					run_step(*it, best3, PN);
				else
					run_step(*it, best2, PN);
			}
		}
		search_best();
		cout << " " << i++ << " " << best1.get_makespan() << endl;
		if (tmp == best1)
		{
			j++;
			if (j == cycle)
				is_true = false;
		}
		else
			j = 0;
		tmp = best1;
	}
}

void PSO::run_mode2()
{
	int i = 1, j = 0;
	bool is_true = true;
	search_best();
	Solution tmp = best1;
	while (is_true)
	{
		vector<Solution>::iterator it = population.begin();
		for (; it != population.end(); it++)
		{
			double tmp = 1.0 * rand() / RAND_MAX * 1.0;
			if (tmp < PM)
				run_step(*it, PN);
			else
			{
				double rn = r.getRandomD(0, 1);
				if (rn < 0.5)
					run_step(*it, best1, PN);
				else if (rn > 0.8)
					run_step(*it, best3, PN);
				else
					run_step(*it, best2, PN);
			}
		}
		search_best();
		cout << " " << i++ << " " << best1.get_makespan() << endl;
		if (tmp == best1)
		{
			j++;
			if (j == maxStandCycle)
				is_true = false;
		}
		else
			j = 0;
		tmp = best1;
	}
}

void PSO::run_mode3(double time)
{
	clock_t startTime, nowTime; //开始时间与当前时间
	startTime = clock();
	nowTime = clock();
	search_best();
	int i = 1;
	while ((double)(nowTime - startTime) / CLOCKS_PER_SEC < time)
	{
		vector<Solution>::iterator it = population.begin();
		for (; it != population.end(); it++)
		{
			double tmp = r.getRandomD(0, 1);
			if (tmp < PM)
				run_step(*it, PN);
			else
			{
				double rn = r.getRandomD(0, 1);
				if (rn < 0.5)
					run_step(*it, best1, PN);
				else if (rn > 0.8)
					run_step(*it, best3, PN);
				else
					run_step(*it, best2, PN);
			}
		}
		search_best();
		cout << " " << i++ << " " << best1.get_makespan() << endl;
		nowTime = clock();
	}
}

void PSO::run_mode3()
{
	clock_t startTime, nowTime; //开始时间与当前时间
	startTime = clock();
	nowTime = clock();
	search_best();
	int i = 1;
	while ((double)(nowTime - startTime) / CLOCKS_PER_SEC < maxRunTime)
	{
		vector<Solution>::iterator it = population.begin();
		for (; it != population.end(); it++)
		{
			double tmp = r.getRandomD(0, 1);
			if (tmp < PM)
				run_step(*it, PN);
			else
			{
				double rn = r.getRandomD(0, 1);
				if (rn < 0.5)
					run_step(*it, best1, PN);
				else if (rn > 0.8)
					run_step(*it, best3, PN);
				else
					run_step(*it, best2, PN);
			}
		}
		search_best();
		cout << " " << i++ << " " << best1.get_makespan() << endl;
		nowTime = clock();
	}
}

void SAA::run_step(Solution& sol, int pn)
{
	vector<int> job_ge = sol.get_job_gene();
	vector<vector<int>> machine_ge = sol.get_machine_gene();
	const int len = (int)job_ge.size();
	for (int i = 0; i < pn; i++)
	{
		//对工序基因变异
		vector<int> pos = r.get_n_diff(0, len - 1, 2);

		int left = pos[0], right = pos[1];
		double rate = r.getRandomD(0, 1);
		if (rate > 0.5)
		{
			int num = job_ge[left];
			if (left > right)
			{
				job_ge.erase(job_ge.begin() + left);
				job_ge.insert(job_ge.begin() + right, num);
			}
			else
			{
				job_ge.insert(job_ge.begin() + right, num);
				job_ge.erase(job_ge.begin() + left);
			}
		}
		else
		{
			swap(job_ge[left], job_ge[right]);
		}

		//对机器基因变异
		vector<int> p = r.get_n_diff(1, len, 2);
		int p1 = p[0], p2 = p[1];
		int work1 = 1, work2 = 1, job1 = 1, job2 = 1;
		int count1 = 0, count2 = 0;
		while (work1 < best.ins.work_num)
		{
			count1 += (int)sol.get_machine_gene()[work1].size() - 1;
			if (count1 >= p1)
			{
				job1 = sol.get_machine_gene()[work1].size() - 1 + p1 - count1;
				break;
			}
			++work1;
		}
		while (work2 < best.ins.work_num)
		{
			count2 += (int)machine_ge[work2].size() - 1;
			if (count2 >= p2)
			{
				job2 = machine_ge[work2].size() - 1 + p2 - count2;
				break;
			}
			++work2;
		}
		int machine1 = best.ins.message[work1][job1].size();
		vector<int> temp1;
		int min = 70000;
		int ri1 = 1;
		auto b1 = best.ins.message[work1][job1].begin();
		auto e1 = best.ins.message[work1][job1].end();
		for (; b1 != e1; ++b1)temp1.push_back(b1->first);
		for (auto it = best.ins.message[work1][job1].begin(); it != e1; ++it) {

			if (min > it->second) {
				min = it->second;
				ri1 = it->first;
			}
		}
		machine_ge[work1][job1] = ri1;

		int machine2 = best.ins.message[work2][job2].size();
		vector<int> temp2;
		auto b2 = best.ins.message[work2][job2].begin();
		auto e2 = best.ins.message[work2][job2].end();
		for (; b2 != e2; ++b2)temp2.push_back(b2->first);
		int ri2 = r.getRandomI(0, machine2 - 1);
		machine_ge[work2][job2] = temp2[ri2];
	}

	sol.set_arg(job_ge, machine_ge);
}

int SAA::stepLength()
{
	int length = s.get_job_gene().size() - MinStep;
	return (length * pow(e, NowT / StartT - 1) + MinStep);
}

void SAA::run()
{
	while (NowT > EndT)
	{
		for (int i = 0; i < L; i++)
		{
			double l = stepLength();
			//cout << "步长: " << l << endl;
			Solution ds = s;
			run_step(ds, l);
			//对于Solution采用makespan作为关键参数进行大小的比较
			if (ds < s)
				s = ds;
			else if (r.getRandomD(0, 1) < k*NowT / StartT)
				s = ds;
			if (ds < best)
				best = ds;
		}
		cout << NowT << "\t: " << s.get_makespan() << endl;
		NowT *= K;
	}
}

void SAA::initialize(Instance ins)
{
	best.set_instance(ins);
	int work = ins.message.size();								//工件数
	vector<int> work_job(work);									//每个工件的工序数量
	for (int i = 1; i < work; ++i)
		work_job[i] = ins.message[i].size();

	vector<int> temp_job;										//产生工序序列
	for (int i = 1; i < work; ++i)
	{
		for (int j = 1; j < work_job[i]; ++j)
			temp_job.push_back(i);
	}

	vector<vector<int>> temp_machine(work);						//产生机器序列
	for (int i = 1; i < work; ++i)
		temp_machine[i].resize(work_job[i]);

	random_shuffle(temp_job.begin(), temp_job.end());
	for (int w = 1; w < work; ++w)
	{
		for (int j = 1; j < (int)ins.message[w].size(); ++j)
		{
			int machine = ins.message[w][j].size();
			vector<int> temp;
			auto b = ins.message[w][j].begin();
			auto c = ins.message[w][j].end();
			for (; b != c; ++b)temp.push_back(b->first);

			for (int m = 0; m < machine; ++m)
			{
				int ri = r.getRandomI(0, machine - 1);
				temp_machine[w][j] = temp[ri];
			}
		}
	}
	Solution ind(temp_job, temp_machine);
	best = ind;
	s = ind;
}

void SAA::initialize(Instance ins, double _K, int _L, double _k)
{
	K = _K;
	L = _L;
	k = _k;
	initialize(ins);
}

void SAA::show_best()
{
	cout << "最终解: " << best.get_makespan() << "\t最终解甘特图如下: ";
	best.show_gantt();
	cout << "最终解工序基因: " << endl;
	best.show_job_gene();
	cout << "最终解机器基因" << endl;
	best.show_machine_gene();
}

void TS::initialize(Instance ins)
{
	max_cycle = 3.5 * (ins.work_num + ins.machine_num);
	max_stand_cycle = ins.work_num + ins.machine_num;
	max_run_time = ins.work_num + ins.machine_num;
	best.set_instance(ins);
	int work = ins.message.size();								//工件数
	vector<int> work_job(work);									//每个工件的工序数量
	for (int i = 1; i < work; ++i)
		work_job[i] = ins.message[i].size();

	vector<int> temp_job;										//产生工序序列
	for (int i = 1; i < work; ++i)
	{
		for (int j = 1; j < work_job[i]; ++j)
			temp_job.push_back(i);
	}

	vector<vector<int>> temp_machine(work);						//产生机器序列
	for (int i = 1; i < work; ++i)
		temp_machine[i].resize(work_job[i]);

	random_shuffle(temp_job.begin(), temp_job.end());
	for (int w = 1; w < work; ++w)
	{
		for (int j = 1; j < (int)ins.message[w].size(); ++j)
		{
			int machine = ins.message[w][j].size();
			vector<int> temp;
			auto b = ins.message[w][j].begin();
			auto c = ins.message[w][j].end();
			for (; b != c; ++b)temp.push_back(b->first);

			for (int m = 0; m < machine; ++m)
			{
				int ri = r.getRandomI(0, machine - 1);
				temp_machine[w][j] = temp[ri];
			}
		}
	}
	Solution ind(temp_job, temp_machine);
	best = ind;
	current = ind;
}

bool TS::is_in_tabuTable(const Solution & sol)
{
	auto it = tabuTable.begin();
	bool flag = false;
	for (; (it != tabuTable.end()) && (flag == false); it++)
		if ((*it) == sol)
			flag = true;
	return flag;
}

void TS::run_mode1(int maxCycle)
{
	for (int i = 0; i < maxCycle; i++)
	{
		tabuTable.push_back(current);
		if (tabuTable.size() > tabuLength)
			tabuTable.erase(tabuTable.begin());
		vector<Solution> neighborhood = current.search_neighborhood(job_num, mac_num);
		for (auto it = neighborhood.begin(); it != neighborhood.end(); it++)
		{
			if (!is_in_tabuTable(*it))
			{
				current = *it;
				break;
			}
		}
		if (best > current)
			best = current;
		show_cycle(i + 1);
	}
}

void TS::run_mode1()
{
	for (int i = 0; i < max_cycle; i++)
	{
		tabuTable.push_back(current);
		if (tabuTable.size() > tabuLength)
			tabuTable.erase(tabuTable.begin());
		vector<Solution> neighborhood = current.search_neighborhood(job_num, mac_num);
		for (auto it = neighborhood.begin(); it != neighborhood.end(); it++)
		{
			if (!is_in_tabuTable(*it))
			{
				current = *it;
				break;
			}
		}
		if (best > current)
			best = current;
		show_cycle(i + 1);
	}
}

void TS::run_mode2(int maxStandCycle)
{
	int i = 0, j = 0;//i用以统计最优解停滞代数,j用以统计迭代次数
	while (true)
	{
		tabuTable.push_back(current);
		if (tabuTable.size() > tabuLength)
			tabuTable.erase(tabuTable.begin());
		vector<Solution> neighborhood = current.search_neighborhood(job_num, mac_num);
		for (auto it = neighborhood.begin(); it != neighborhood.end(); it++)
		{
			if (!is_in_tabuTable(*it))
			{
				current = *it;
				break;
			}
		}
		if (best > current)
		{
			best = current;
			i = 0;
		}
		else
			i++;
		show_cycle(++j);
		if (i == (maxStandCycle + tabuLength))
			break;
	}
}

void TS::run_mode2()
{
	int i = 0, j = 0;//i用以统计最优解停滞代数,j用以统计迭代次数
	while (true)
	{
		tabuTable.push_back(current);
		if (tabuTable.size() > tabuLength)
			tabuTable.erase(tabuTable.begin());
		vector<Solution> neighborhood = current.search_neighborhood(job_num, mac_num);
		for (auto it = neighborhood.begin(); it != neighborhood.end(); it++)
		{
			if (!is_in_tabuTable(*it))
			{
				current = *it;
				break;
			}
		}
		if (best > current)
		{
			best = current;
			i = 0;
		}
		else
			i++;
		show_cycle(++j);
		if (i == (max_stand_cycle + tabuLength))
			break;
	}
}

void TS::run_mode3(double maxRunTime)
{
	clock_t startTime, nowTime; //开始时间与当前时间
	startTime = clock();
	int i = 0;	//迭代次数
	while (true)
	{
		tabuTable.push_back(current);
		if (tabuTable.size() > tabuLength)
			tabuTable.erase(tabuTable.begin());
		vector<Solution> neighborhood = current.search_neighborhood(job_num, mac_num);
		for (auto it = neighborhood.begin(); it != neighborhood.end(); it++)
		{
			if (!is_in_tabuTable(*it))
			{
				current = *it;
				break;
			}
		}
		if (best > current)
			best = current;
		show_cycle(++i);
		nowTime = clock();
		if ((double)(nowTime - startTime) / CLOCKS_PER_SEC > maxRunTime)
			break;
	}
}

void TS::run_mode3()
{
	clock_t startTime, nowTime; //开始时间与当前时间
	startTime = clock();
	int i = 0;	//迭代次数
	while (true)
	{
		tabuTable.push_back(current);
		if (tabuTable.size() > tabuLength)
			tabuTable.erase(tabuTable.begin());
		vector<Solution> neighborhood = current.search_neighborhood(job_num, mac_num);
		for (auto it = neighborhood.begin(); it != neighborhood.end(); it++)
		{
			if (!is_in_tabuTable(*it))
			{
				current = *it;
				break;
			}
		}
		if (best > current)
			best = current;
		show_cycle(++i);
		nowTime = clock();
		if ((double)(nowTime - startTime) / CLOCKS_PER_SEC > max_run_time)
			break;
	}
}

void TS::show_cycle(int i)
{
	cout << "第 " << i << " 代 : " << endl;
	cout << "\t当前解 = " << current.get_makespan() << endl;
	//current.show_job_gene();
	cout << "\t最优解 = " << best.get_makespan() << endl;
	//best.show_job_gene();
}
