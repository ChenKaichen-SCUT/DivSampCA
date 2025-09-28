#include "cdclcasampler.h"
using std::cout;
CDCLCASampler::CDCLCASampler(string cnf_file_path, int seed): rnd_file_id_gen_(seed){
    cnf_file_path_ = cnf_file_path;
    seed_ = seed;
    SetDefaultPara();
}

CDCLCASampler::~CDCLCASampler(){
    delete cdcl_sampler;
}

void CDCLCASampler::SetTransferSize(int transfer_set_size){
    transfer_count = transfer_set_size;
}

void CDCLCASampler::SetCandidateSetSize(int candidate_set_size){
    candidate_set_size_ = candidate_set_size;
    candidate_sample_init_solution_set_.resize(candidate_set_size_);
    candidate_testcase_set_.resize(candidate_set_size_);
}

void CDCLCASampler::SetDefaultPara(){
    t_wise_ = 2;
    transfer_count = 100;
    candidate_set_size_ = 100;
    candidate_sample_init_solution_set_.resize(candidate_set_size_);
    candidate_testcase_set_.resize(candidate_set_size_);

    int pos = cnf_file_path_.find_last_of('/');
    cnf_instance_name_ = cnf_file_path_.substr(pos + 1);
    size_t suffix_pos = cnf_instance_name_.rfind(".cnf");
    if (suffix_pos != string::npos){
        cnf_instance_name_.replace(suffix_pos, 4, "");
    }

    reduced_cnf_file_path_ = "./tmp/" + get_random_prefix() + "_reduced.cnf";
    testcase_set_save_path_ = "./testcase_set/" + cnf_instance_name_ + "_testcase_set.txt";
}

void CDCLCASampler::Init(){
    read_cnf();
    original_num_var_ = num_var_;
    ReduceCNF();
    read_cnf();

    Init2TupleInfo();
    cdcl_sampler = new ExtMinisat::SamplingSolver(num_var_, clauses, seed_, true, false);
    GenerateInitTestcaseCDCL();
    InitSampleWeightByAppearance();
    sample_proportion.resize(num_var_);
}

void CDCLCASampler::ReduceCNF(){

    vector<pair<int,int>> to_erase;

    if(extract_backbones(cnf_file_path_,backbones)) processbackbones(clauses,backbones);
    for(auto clause : clauses){
        if(clause.size() == 2)
        {
            double_clause_Set.insert({clause[0],clause[1]});
        }
        else if (clause.size() > 2){
            tri_clauses.emplace_back(clause);
        } 
    }

    for (const auto& p : double_clause_Set) {
        int a = p.first;
        int b = p.second;
        if (a * b < 0) {
            auto k1 = make_pair(-a, -b);
            auto k2 = make_pair(-b, -a);

            if (double_clause_Set.count(k1)) {
                equivalence_set.emplace(abs(a), abs(b));
                to_erase.push_back(k1);
                to_erase.push_back(p); 
            } else if (double_clause_Set.count(k2)) {
                equivalence_set.emplace(abs(a), abs(b));
                to_erase.push_back(k2);
                to_erase.push_back(p);
            }
        }
    }

    for (const auto& k : to_erase)
        double_clause_Set.erase(k);
    
    
    for (auto& p : equivalence_set) {
        dsu.unionSets(p.first, p.second);
        equivalence_number_set.emplace(p.first);
        equivalence_number_set.emplace(p.second);   
    }

    for (int element : equivalence_number_set) {
        int representative = dsu.find(element);
        equivalence_map[representative].push_back(element);
    }

    unordered_set<string> no_header_clauses;
    vector<string> header_clauses;

    for (const auto& p : double_clause_Set) {

        string temp_string;
        auto first_abs = abs(p.first);
        auto second_abs = abs(p.second);

        auto it_first = equivalence_number_set.find(first_abs);
        auto it_second = equivalence_number_set.find(second_abs);

        if (it_first == equivalence_number_set.end()) {
            temp_string += to_string(p.first);
        } else {
            temp_string += to_string (  ((p.first< 0)? -1 : 1)   * dsu.find(first_abs)  );
        }
        temp_string += " ";

        if (it_second == equivalence_number_set.end()) {
            temp_string += to_string(p.second);
        } else {
            temp_string += to_string(((p.second< 0) ? -1 : 1) * dsu.find(second_abs));
        }

        temp_string += " 0";
        no_header_clauses.emplace(move(temp_string));
    }

    for (const auto& clause : tri_clauses) {
        string temp_string;

        unordered_set<int> already_string;
        already_string.reserve(clause.size()); 

        for (int literal : clause) {
            int abs_literal = abs(literal);
            auto it = equivalence_number_set.find(abs_literal);
            int num;

            if (it == equivalence_number_set.end()) {
                num = literal;
            } else {
                num = ((literal< 0) ? -1 : 1) * dsu.find(abs_literal);
            }

            if (already_string.emplace(num).second) { 
                temp_string += to_string(num);
                temp_string += " ";
            }
        }
        temp_string += "0";
        no_header_clauses.emplace(move(temp_string));
    }

    for (const auto& line : no_header_clauses) {
        istringstream iss(line);
        string clause;
        clause.reserve(line.size());
        int num;

        while (iss >> num) {
            if (num == 0) break; 

            int absNum = std::abs(num);
            int sign = (num < 0) ? -1 : 1;

            auto [it, inserted] = varSignMap.try_emplace(absNum, varSignMap.size() + 1);
            clause += std::to_string(sign * it->second) + ' ';
        }

        clause += '0'; 
        header_clauses.emplace_back(std::move(clause));
    }

    for(auto m : varSignMap){
        To_original_varSignMap.emplace(m.second,m.first);
    }

    string newLine = "p cnf "+ to_string(varSignMap.size()) + " "+ to_string(no_header_clauses.size());
    tmp_cnf_file = "./simplified_"+cnf_instance_name_+".cnf";
    ofstream outFile(tmp_cnf_file);
    outFile << newLine << '\n';
    for (const auto& l : header_clauses) {
        outFile << l << '\n';
    }
    outFile.close();

    string cmd = "./bin/coprocessor -enabled_cp3 -up -subsimp -no-bve -no-bce"
                 " -no-dense -dimacs=" + reduced_cnf_file_path_ + " " + tmp_cnf_file;

    int return_val = system(cmd.c_str());

    cnf_file_path_ = reduced_cnf_file_path_;
}

void CDCLCASampler::read_cnf_header(ifstream& ifs, int& nvar, int& nclauses){
    string line;
    istringstream iss;
    
    while (getline(ifs, line)){
        if (line.substr(0, 1) == "c")
            continue;
        else if (line.substr(0, 1) == "p"){
            string tempstr1, tempstr2;
            iss.clear();
            iss.str(line);
            iss.seekg(0, ios::beg);
            iss >> tempstr1 >> tempstr2 >> nvar >> nclauses;
            break;
        }
    }
}

bool CDCLCASampler::check_no_clauses(){
    ifstream fin(cnf_file_path_.c_str());
    if (!fin.is_open()) return true;

    int num_clauses_original;
    read_cnf_header(fin, num_var_, num_clauses_original);

    fin.close();
    return num_clauses_original <= 0;
}

bool CDCLCASampler::read_cnf(){

    clauses.clear();
    
    ifstream fin(cnf_file_path_.c_str());
    if (!fin.is_open()) return false;

    read_cnf_header(fin, num_var_, num_clauses_);

    if (num_var_ <= 0 || num_clauses_ < 0){
        fin.close();
        return false;
    }
    clauses.resize(num_clauses_);

    for (int c = 0; c < num_clauses_; ++c)
    {
        int cur_lit;
        fin >> cur_lit;
        while (cur_lit != 0)
        {
            clauses[c].emplace_back(cur_lit);
            fin >> cur_lit;
        }
    }

    fin.close();

    return true;
}

void CDCLCASampler::GenerateInitTestcaseCDCL(){
    vector<pair<int, int> > sample_prob = get_prob_in(true);
    cdcl_sampler->set_prob(sample_prob);
    vector<int> test_case_to_add_ = cdcl_sampler->get_solution();
    testcase_set_.emplace_back(test_case_to_add_);
    SimpleBitSet tmp(num_tuple_all_possible_);
    for (int i = 0; i < num_var_ - 1; i++)
    {
        for (int j = i + 1; j < num_var_; j++)
        {
            long long index_tuple = Get2TupleMapIndex(i, test_case_to_add_[i], j, test_case_to_add_[j]);
                tmp.set(index_tuple);
        }
    }
    testcase_SimpleBitset.emplace_back(move(tmp));
    num_generated_testcase_ = 1;
}

void CDCLCASampler::Init2TupleInfo(){
    num_combination_all_possible_ = (long long)num_var_ * (num_var_ - 1) / 2;
    num_tuple_all_possible_ = num_combination_all_possible_ * 4;
    already_t_wise = SimpleBitSet(num_tuple_all_possible_);
    num_tuple_ = 0;
}

void CDCLCASampler::Update2TupleInfo(){
    int index_testcase = num_generated_testcase_ - 1;
    already_t_wise.bitwise_or(testcase_SimpleBitset[index_testcase]);
    num_tuple_ = already_t_wise.count();
}

void CDCLCASampler::InitSampleWeightByAppearance(){
    var_positive_appearance_count_.resize(num_var_);
    var_positive_sample_weight_.resize(num_var_);
}

void CDCLCASampler::UpdateSampleWeightByAppearance(){
    int new_testcase_index = num_generated_testcase_ - 1;
    const vector<int>& new_testcase = testcase_set_[new_testcase_index];
    for (int v = 0; v < num_var_; v++)
    {
        var_positive_appearance_count_[v] += new_testcase[v];
        var_positive_sample_weight_[v] = 1. - double(var_positive_appearance_count_[v]) / num_generated_testcase_;
    }
}

void CDCLCASampler::processbackbones(vector<vector<int>>& clauses, const vector<int>& backbones) {
    unordered_set<int> backboneSet(backbones.begin(), backbones.end());

    auto it = clauses.begin();
    while (it != clauses.end()) {
        bool deleteClause = false;
        auto litIt = it->begin();
        while (litIt != it->end()) {
            int lit = *litIt;
            if (backboneSet.find(lit) != backboneSet.end()) {
                deleteClause = true;
                break;
            } else if (backboneSet.find(-lit) != backboneSet.end()) {
                litIt = it->erase(litIt);
            } else {
                ++litIt;
            }
        }
        if (deleteClause) {
            it = clauses.erase(it);
        } else {
            ++it;
        }
    }
}

void CDCLCASampler::GenerateCandidateTestcaseSet(){
    vector<pair<float, float>> sample_prob1 = new_get_prob(false);
    vector<pair<int, int>> sample_prob;

    for(int i =0;i<sample_prob1.size();i++)
    {
        sample_prob.emplace_back((int)sample_prob1[i].first,(int)sample_prob1[i].second);
    }
    
    for (int i = 0 ; i < candidate_set_size_; i++){
        cdcl_sampler->set_prob(sample_prob);
        cdcl_sampler->get_solution(candidate_testcase_set_[i]);
    }
}

int CDCLCASampler::SelectTestcaseFromCandidateSetByTupleNum(){
    for (int i = 0; i < candidate_set_size_; ++i){
        pq.emplace_back(candidate_testcase_set_[i], get_gain(candidate_testcase_set_[i]));
        pq_idx.push_back(0);
    }
    
    iota(pq_idx.begin(), pq_idx.end(), 0);
    sort(pq_idx.begin(), pq_idx.end(), [&](int i, int j){
        return pq[i].second.count() > pq[j].second.count();
    });

    return pq_idx[0];
}

void CDCLCASampler::GenerateTestcase(){
    GenerateCandidateTestcaseSet();
    selected_candidate_index_ = SelectTestcaseFromCandidateSetByTupleNum();

    int testcase_index = num_generated_testcase_;
    vector<int> test_case_to_add_  = pq[selected_candidate_index_].first;
    testcase_set_.emplace_back(test_case_to_add_);
    selected_candidate_bitset_ = pq[selected_candidate_index_].second;

    SimpleBitSet tmp(num_tuple_all_possible_);
    for (int i = 0; i < num_var_ - 1; i++)
    {
        for (int j = i + 1; j < num_var_; j++)
        {
            long long index_tuple = Get2TupleMapIndex(i, test_case_to_add_[i], j, test_case_to_add_[j]);
                tmp.set(index_tuple);
        }
    }
    testcase_SimpleBitset.emplace_back(move(tmp));
}

long long CDCLCASampler::Get2TupleMapIndex(long i, long v_i, long j, long v_j){
    long long base = (v_i << 1 | v_j) * num_combination_all_possible_;
    long long pos = (2ll * num_var_ - i - 1) * i / 2 + j - i - 1;
    return base + pos;
}

void CDCLCASampler::GenerateCoveringArray(){

    auto start_time = chrono::system_clock::now().time_since_epoch();
    Init();
    
    for (num_generated_testcase_ = 1; ; num_generated_testcase_++)
    {
        Update2TupleInfo();

        if (num_generated_testcase_ > 1) {
            clear_pq();
        }
        cout << num_generated_testcase_ << endl;
        UpdateSampleWeightByAppearance();
        GenerateTestcase();
        if (selected_candidate_bitset_.count() < transfer_count){
            testcase_set_.pop_back();
            testcase_SimpleBitset.pop_back();
            break;
        }
    }
     
    Reduce_redundancy(false);

    clear_final();
    Reduce_redundancy(true);
    vector<vector<int>> new_test_case_set_;
    for(auto test_case : testcase_set_){
        vector<int> tmp_test_case_set_(original_num_var_);
        for(int i = 1 ; i <= num_var_; i++)
        {
            int original = To_original_varSignMap[i];
            if(equivalence_map.find(original) != equivalence_map.end())
            {
                for(auto eq: equivalence_map[original])
                {
                    tmp_test_case_set_[eq -1] = test_case[i-1];
                }
            }
            else{
                tmp_test_case_set_[original - 1] = test_case[i-1];
            }
        }
        for(int i =0 ; i<backbones.size();i++)
        {
            if(backbones[i] >0)
            {
                tmp_test_case_set_[backbones[i] - 1 ] =1;
            }
            else{
                tmp_test_case_set_[abs(backbones[i])-1] =0;
            }
        }
        new_test_case_set_.emplace_back(tmp_test_case_set_);
    }
    
    SaveTestcaseSet(new_test_case_set_, testcase_set_save_path_);

    auto end_time = chrono::system_clock::now().time_since_epoch();
    cpu_time_ = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() / 1000.0;
    cout << "c Generate testcase set finished, containing " << testcase_set_.size() << " testcases!" << endl;
    cout << "c CPU time cost by generating testcase set: " << cpu_time_ << " seconds" << endl;

    remove_temp_files();

}

void CDCLCASampler::SaveTestcaseSet(vector<vector<int>> new_test_case_set_ ,string result_path){
    ofstream res_file(result_path);
    for (const vector<int>& testcase: new_test_case_set_)
    {
        for (int v = 0; v < original_num_var_; v++)
        {
            res_file << testcase[v] << " ";
        }
        res_file << endl;
    }
    res_file.close();
    cout << "c Testcase set saved in " << result_path << endl;
}

vector<pair<int, int> > CDCLCASampler::get_prob_in(bool init){
    vector<pair<int, int> > res;
    res.reserve(num_var_);
    if (!init){
        for (int i = 0; i < num_var_; ++i){
            int v1 = num_generated_testcase_ - var_positive_appearance_count_[i];
            int v2 = var_positive_appearance_count_[i];
            res.emplace_back(v1, v2);
        }
    }
    return res;
}

SimpleBitSet CDCLCASampler::get_gain(const vector<int>& asgn){
    SimpleBitSet res(num_tuple_all_possible_);
    for (int i = 0; i < num_var_ - 1; i++)
    {
        for (int j = i + 1; j < num_var_; j++)
        {
            long long index_tuple = Get2TupleMapIndex(i, asgn[i], j, asgn[j]);
            if (!already_t_wise.get(index_tuple)){
                res.set(index_tuple);
            }
        }
    }
    return res;
}

void CDCLCASampler::remove_temp_files(){
    string cmd;
    int ret;
    cmd = "rm " + reduced_cnf_file_path_;
    ret = system(cmd.c_str());
    cmd = "rm " +  tmp_cnf_file;
    ret = system(cmd.c_str());
}

string CDCLCASampler::get_random_prefix(){
    return cnf_instance_name_ + to_string(getpid()) + to_string(rnd_file_id_gen_());
}

void CDCLCASampler::clear_pq(){
    for (auto& bs: pq){
        bs.second.difference_(selected_candidate_bitset_);
    }

    sort(pq.begin(), pq.end(), [](const pair<vector<int>, SimpleBitSet>& si, const pair<vector<int>, SimpleBitSet>& sj){
        return si.second.count() > sj.second.count();
    });

    int cur_pqsize = (int) pq.size();
    while (cur_pqsize > candidate_set_size_ || (cur_pqsize > 0 && pq.back().second.count() == 0)){
        pq.pop_back();
        pq_idx.pop_back();
        --cur_pqsize;
    }
}

void CDCLCASampler::find_uncovered_tuples(bool simplify){
    for (int i = 0; i < num_var_ - 1; ++i){
        for (int v_i = 0; v_i < 2; ++v_i){
            for (int j = i + 1; j < num_var_; ++j){
                for (int v_j = 0; v_j < 2; ++v_j){
                    long long index_tuple = Get2TupleMapIndex(i, v_i, j, v_j);
                    if (!already_t_wise.get(index_tuple)){
                        bool flag = true;
                        if (simplify){
                            cdcl_solver->add_assumption(i, v_i);
                            cdcl_solver->add_assumption(j, v_j);
                            bool res = cdcl_solver->solve();
                            if (!res){
                                --num_tuple_all_exact_;
                                flag = false;
                            } else {
                                vector<int> tc(num_var_, 0);
                                get_cadical_solution(tc);
                            }
                            cdcl_solver->clear_assumptions();
                        }
                        if (flag){
                            uncovered_tuples.emplace_back(vector<int>{v_i ? (i+1):(-i-1), v_j ? (j+1):(-j-1)});
                        }
                    }
                }
            }
        }
    }
}

void CDCLCASampler::get_cadical_solution(vector<int>& tc){
    cdcl_solver->get_solution(tc);
}

void CDCLCASampler::new_choose_final_plain(){

    find_uncovered_tuples(true);

    random_device rd;  
	mt19937 g(rd());	
	shuffle(uncovered_tuples.begin(), uncovered_tuples.end(), g);  

    CDCLSolver::Solver* uncovered_filler = new CDCLSolver::Solver;
    uncovered_filler->read_clauses(num_var_, clauses);
	
	vector<int> satisifed_index;	
	unordered_map<int,int> assu;
	unordered_map<int,int> old_assu;

	while(uncovered_tuples.size()>0)
	{
		uncovered_filler->clear_assumptions();
		assu.clear();old_assu.clear();
		satisifed_index.clear();

		for(int j = 0;j<uncovered_tuples.size();j++)
		{
            old_assu = assu;

            if(assu.find(abs(uncovered_tuples[j][0]) - 1) != assu.end() && assu[abs(uncovered_tuples[j][0]) - 1] != (uncovered_tuples[j][0] < 0 ? 0: 1))
            {
                assu = old_assu;
                continue;
            }
            else{
                assu.insert({ abs(uncovered_tuples[j][0]) - 1, uncovered_tuples[j][0] < 0 ? 0: 1});
            }

            if(assu.find(abs(uncovered_tuples[j][1]) - 1) != assu.end() && assu[abs(uncovered_tuples[j][1]) - 1] != (uncovered_tuples[j][1] < 0 ? 0: 1))
            {
                assu = old_assu;
                continue;
            }
            else{
                assu.insert({ abs(uncovered_tuples[j][1]) - 1, uncovered_tuples[j][1] < 0 ? 0: 1});
            }
			        
			satisifed_index.emplace_back(j);

			if(assu.size() > old_assu.size())
			{
				unordered_map<int,int> difference = difference_(assu,old_assu);
				for(const auto pair: difference)
				{
					uncovered_filler->add_assumption(pair.first ,pair.second);
				}
				if(!uncovered_filler->solve1(num_var_, clauses))
				{
					for(int m = 0; m < difference.size();m++)
					{
						uncovered_filler->assumptions.pop_back();
					}
					satisifed_index.pop_back();
					assu = old_assu;
				}
			}
        }
        
        uncovered_filler->clear_assumptions();
        for(auto pair:assu){uncovered_filler->add_assumption(pair.first ,pair.second);}
        bool res =uncovered_filler ->solve();
        if(!res)
        {
            continue;
        }
        
        vector<int> new_tc;
        uncovered_filler->get_solution(new_tc);
        testcase_set_.emplace_back(new_tc);

        SimpleBitSet tmp(num_tuple_all_possible_);
        for (int i = 0; i < num_var_ - 1; i++)
        {
            for (int j = i + 1; j < num_var_; j++)
            {
                long long index_tuple = Get2TupleMapIndex(i, new_tc[i], j, new_tc[j]);
                    tmp.set(index_tuple);
            }
        }
        testcase_SimpleBitset.emplace_back(move(tmp));


        ++num_generated_testcase_;
        Update2TupleInfo1();
        cout << num_generated_testcase_ << endl;

        for (auto it = satisifed_index.rbegin(); it != satisifed_index.rend(); ++it) {  
            int index = *it;  
            if (index >= 0 && index < uncovered_tuples.size()) {  
            uncovered_tuples.erase(uncovered_tuples.begin() + index);  
            }  
        }  
       
	}
    
    delete uncovered_filler;

}

void CDCLCASampler::Update2TupleInfo1(){
    int index_testcase = num_generated_testcase_ - 1;
    const vector<int>& testcase = testcase_set_[index_testcase];

    for (int i = 0; i < num_var_ - 1; i++)
    {
        for (int j = i + 1; j < num_var_; j++)
        {
            long long index_tuple = Get2TupleMapIndex(i, testcase[i], j, testcase[j]);
            bool res = already_t_wise.set(index_tuple);
            if (res)
            {
                num_tuple_++;
            }
        }
    }
}

bool CDCLCASampler::extract_backbones(const string& cnf_path, vector<int>& nums)   {
    nums.clear();
    string cmd = "./bin/cadiback -q " + cnf_path;
    char buf[512];

    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return false;

    while (fgets(buf, sizeof(buf), pipe)) {
        string line(buf);
        istringstream iss(line);
        string tag;
        int val;
        if (iss >> tag >> val && tag == "b" && val != 0)
            nums.push_back(val);
    }
    pclose(pipe);

    return !nums.empty(); 
}

void CDCLCASampler::clear_final(){
    cout << endl << "c Clear final: fix-up all remaining tuples ..." << endl;

    num_tuple_all_exact_ = num_tuple_all_possible_;
    cdcl_solver = new CDCLSolver::Solver;
    cdcl_solver->read_clauses(num_var_, clauses);

    new_choose_final_plain();
    
    delete cdcl_solver;
}

vector<pair<float, float> > CDCLCASampler::new_get_prob(bool init){

    if (!init){
        for (int i = 0; i < num_var_ - 1 ; ++i)
        {
            for(int j = i + 1 ; j < num_var_ ; ++j)
            {
                bool get00 = already_t_wise.get(Get2TupleMapIndex(i,0, j, 0));
                bool get01 = already_t_wise.get(Get2TupleMapIndex(i,0, j, 1));
                bool get10 = already_t_wise.get(Get2TupleMapIndex(i,1, j, 0));
                bool get11 = already_t_wise.get(Get2TupleMapIndex(i,1, j, 1));
                
                if(! (get00 && get01 && get10 && get11))
                {
                    int tuple0_ = !get00 + !get01;
                    int tuple1_ = !get10 + !get11; 
                    int tuple_0 = !get00 + !get10; 
                    int tuple_1 = !get01 + !get11; 

                    float v1 = static_cast<float>(tuple1_) / (tuple0_+tuple1_);
                    sample_proportion[i].first += v1;
                    sample_proportion[i].second += 1. - v1;

                    float v2 = static_cast<float>(tuple_1) / (tuple_1+tuple_0);
                    sample_proportion[j].first += v2;
                    sample_proportion[j].second += 1. - v2;
                }
            }
        }
    }

    return sample_proportion;
}

void CDCLCASampler:: Reduce_redundancy(bool last){
    const size_t n = testcase_SimpleBitset.size();
    vector<SimpleBitSet> suffix(n, SimpleBitSet(num_tuple_all_possible_));
    suffix[n - 1] = testcase_SimpleBitset[n - 1];
    for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
        suffix[i] = suffix[i + 1];
        suffix[i].bitwise_or(testcase_SimpleBitset[i]);
    }

    SimpleBitSet prefix(num_tuple_all_possible_);
    vector<long long> uniq_cnt(n);

    for (size_t i = 0; i < n; ++i) {
        SimpleBitSet common(num_tuple_all_possible_);
        if (i > 0) common = prefix;               
        if (i < n - 1) common.bitwise_or(suffix[i + 1]);
        SimpleBitSet uniq = testcase_SimpleBitset[i];
        uniq.difference_(common);                 
        uniq_cnt[i] = uniq.count();
        prefix.bitwise_or(testcase_SimpleBitset[i]);
    }

    vector<size_t> to_del;
    if(last)
    {
        for (size_t i = 0; i < uniq_cnt.size(); ++i) {
            if (uniq_cnt[i] == 0) {
                to_del.push_back(i);
            }
        }
    }
    else{
        for (size_t i = 0; i < uniq_cnt.size(); ++i) {
            if (uniq_cnt[i] < reduce_size ) {
                to_del.push_back(i);
            }
        }
    }

    sort(to_del.rbegin(), to_del.rend());
    for (size_t idx : to_del) {
        testcase_set_.erase(testcase_set_.begin() + idx);
        testcase_SimpleBitset.erase(testcase_SimpleBitset.begin() + idx);
        num_generated_testcase_--;
    }

    already_t_wise.clear();
    for(int i = 0; i<testcase_SimpleBitset.size();i++)
    {
        already_t_wise.bitwise_or(testcase_SimpleBitset[i]);
    }
    num_tuple_ = already_t_wise.count();
}

unordered_map<int, int> CDCLCASampler::difference_(const unordered_map<int, int>& map1, const unordered_map<int, int>& map2) {  
    unordered_map<int, int> diff;  
    for (const auto& pair : map1) {  
        if (map2.find(pair.first) == map2.end()) {  
            diff[pair.first] = pair.second;  
        }  
    }  
    return diff;  
}  
