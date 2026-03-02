#include "cdclcasampler.h"
#include <queue>
using std::cout;

inline int lit2idx(int lit) {
    int v = abs(lit) - 1;
    return (lit > 0) ? 2*v : 2*v + 1;
}

inline int idx2lit(int idx) {
    int v = idx / 2 + 1;
    return (idx % 2 == 0) ? v : -v;
}

inline int neg(int idx) {
    return idx ^ 1;
}

CDCLCASampler::CDCLCASampler(string cnf_file_path, int seed){
    cnf_file_path_ = cnf_file_path;
    seed_ = seed;
    SetDefaultPara();
}

CDCLCASampler::~CDCLCASampler(){
    delete cdcl_sampler;
}

void CDCLCASampler::SetTransferSize(int transfer_size){
    transfer_size_ = transfer_size;
}

void CDCLCASampler::SetCandidateSetSize(int candidate_set_size){
    candidate_set_size_ = candidate_set_size;
    candidate_testcase_set_.resize(candidate_set_size_);
}

void CDCLCASampler::SetDefaultPara(){
    transfer_size_ = 150;
    candidate_set_size_ = 150;
    candidate_testcase_set_.resize(candidate_set_size_);

    int pos = cnf_file_path_.find_last_of('/');
    cnf_instance_name_ = cnf_file_path_.substr(pos + 1);
    size_t suffix_pos = cnf_instance_name_.rfind(".cnf");
    if (suffix_pos != string::npos){
        cnf_instance_name_.replace(suffix_pos, 4, "");
    }

    reduced_cnf_file_path_ = "./reduced_cnf_file/" + cnf_instance_name_ + "_reduced.cnf";
    testcase_set_save_path_ = "./testcase_set/" + cnf_instance_name_ + "_testcase_set.txt";
}

void CDCLCASampler::GenerateCoveringArray(){
    cout << "c Sampling started" << endl;
    auto start_time = chrono::system_clock::now().time_since_epoch();
    Init();
    for (num_generated_testcase_ = 1; ; num_generated_testcase_++){
        cout<<num_generated_testcase_<<endl;
        Update_already_t_wise();
        if (num_generated_testcase_ > 1) {
            clear_pq();
        }
        GenerateTestcase();
        if (selected_candidate_bitset_.count() <transfer_size_){
            testcase_set_.pop_back();
            break;
        }
    }
    
    clear_final();
    SaveTestcaseSet(final_testcase_set_, testcase_set_save_path_);
    remove_temp_files();

    auto end_time = chrono::system_clock::now().time_since_epoch();
    cpu_time_ = chrono::duration_cast<chrono::milliseconds>(end_time - start_time).count() / 1000.0;    
    cout << "c Generate testcase set finished, containing " << testcase_set_.size() << " testcases!" << endl;
    cout << "c CPU time cost by generating testcase set: " << cpu_time_ << " seconds" << endl;
}

void CDCLCASampler::Init(){
    read_cnf();
    original_num_var_ = num_var_;
    ReduceCNF();
    read_cnf();

    Init2TupleInfo();
    cdcl_sampler = new ExtMinisat::SamplingSolver(num_var_, clauses, seed_, true, false);
    GenerateInitTestcaseCDCL();
}

bool CDCLCASampler::read_cnf(){
    num_var_ = 0;
    num_clauses_ = 0;
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

void CDCLCASampler::Init2TupleInfo(){
    num_combination_all_possible_ = (long long)num_var_ * (num_var_ - 1) / 2;
    num_tuple_all_possible_ = num_combination_all_possible_ * 4;
    already_t_wise = SimpleBitSet(num_tuple_all_possible_);
    sample_proportion.resize(num_var_);
    InitCombinationOffset();
    num_tuple_ = 0;
}

void CDCLCASampler::GenerateInitTestcaseCDCL(){
    vector<pair<int, int>> res;
    res.reserve(num_var_);
    cdcl_sampler->set_prob(res);

    vector<int> init_test_case_ = cdcl_sampler->get_solution();
    testcase_set_.emplace_back(init_test_case_);
    num_generated_testcase_ = 1;
}

bool CDCLCASampler::extract_backbones(const string& cnf_path, vector<int>& nums){
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

void CDCLCASampler::processbackbones(vector<vector<int>>& clauses, const vector<int>& backbones){
    vector<int8_t> backboneSign(num_var_ + 1, 0);

    for (int b : backbones) {
        int v = std::abs(b);
        backboneSign[v] = (b > 0) ? 1 : -1;
    }

    vector<vector<int>> newClauses;
    newClauses.reserve(clauses.size());

    for (const auto& clause : clauses) {
        bool satisfied = false;
        vector<int> newClause;
        newClause.reserve(clause.size());

        for (int lit : clause) {
            int var  = std::abs(lit);
            int sign = (lit > 0) ? 1 : -1;

            int8_t bsign = backboneSign[var];

            if (bsign == sign) {
                satisfied = true;
                break;
            }

            if (bsign == 0) {
                newClause.push_back(lit);
            }
        }

        if (!satisfied) {
            newClauses.emplace_back(std::move(newClause));
        }
    }

    clauses.swap(newClauses);
}

void CDCLCASampler::ReduceCNF() {
    if(extract_backbones(cnf_file_path_,backbones)) processbackbones(clauses,backbones);
    vector<pair<int,int>> binary;
    vector<vector<int>> longc;

    for (auto& c : clauses) {
        if (c.size() == 2)
            binary.emplace_back(c[0], c[1]);
        else if (c.size() > 2)
            longc.emplace_back(c);
    }

    int N = num_var_ * 2;
    vector<vector<int>> G(N);

    for (auto& [a,b] : binary) {
        G[lit2idx(-a)].push_back(lit2idx(b));
        G[lit2idx(-b)].push_back(lit2idx(a));
    }

    vector<int> dfn(N, -1), low(N), comp(N);
    vector<bool> in(N, false);
    stack<int> st;
    int ts = 0, scc_cnt = 0;

    function<void(int)> dfs = [&](int u) {
        dfn[u] = low[u] = ts++;
        st.push(u);
        in[u] = true;

        for (int v : G[u]) {
            if (dfn[v] == -1) {
                dfs(v);
                low[u] = min(low[u], low[v]);
            } else if (in[v]) {
                low[u] = min(low[u], dfn[v]);
            }
        }

        if (dfn[u] == low[u]) {
            while (true) {
                int x = st.top(); st.pop();
                in[x] = false;
                comp[x] = scc_cnt;
                if (x == u) break;
            }
            scc_cnt++;
        }
    };

    for (int i = 0; i < N; ++i)
        if (dfn[i] == -1) dfs(i);

    for (int v = 1; v <= num_var_; ++v) {
        if (comp[lit2idx(v)] == comp[lit2idx(-v)]) {
            cerr << "UNSAT during preprocessing\n";
            exit(1);
        }
    }

    vector<vector<int>> bucket(scc_cnt);
    for (int i = 0; i < N; ++i)
        bucket[comp[i]].push_back(i);

    unordered_set<int> used_vars;

    for (auto& b : bucket) {
        if (b.size() < 2) continue;

        int base = b[0];
        int base_lit = idx2lit(base);

        for (size_t i = 1; i < b.size(); ++i) {
            int cur_lit = idx2lit(b[i]);

            int x = abs(base_lit);
            int y = abs(cur_lit);

            bool same = ((base_lit > 0) == (cur_lit > 0));
            dsu.unite(x, y, same);

            used_vars.insert(x);
            used_vars.insert(y);
        }
    }

    unordered_set<string> uniq;

    auto rewrite = [&](int lit) {
        int v = abs(lit);
        int r = dsu.find(v);
        int s = dsu.sign(v);
        int final_sign = (lit > 0) ? s : -s;
        return final_sign * r;
    };

    for (auto& [a,b] : binary) {
        int x = rewrite(a), y = rewrite(b);
        if (x == y) continue;
        uniq.insert(to_string(x) + " " + to_string(y) + " 0");
    }

    for (auto& c : longc) {
        unordered_set<int> seen;
        string s;
        for (int lit : c) {
            int r = rewrite(lit);
            if (seen.insert(r).second)
                s += to_string(r) + " ";
        }
        s += "0";
        uniq.insert(move(s));
    }

    unordered_map<int,int> mp;
    vector<string> finalc;

    for (auto& cl : uniq) {
        istringstream iss(cl);
        int x;
        string out;

        while (iss >> x && x != 0) {
            int v = abs(x), sg = (x > 0 ? 1 : -1);
            if (!mp.count(v))
                mp[v] = mp.size() + 1;
            out += to_string(sg * mp[v]) + " ";
        }
        out += "0";
        finalc.push_back(move(out));
    }

    new2old.assign(mp.size() + 1, 0);

    for (auto& [old_rep, new_id] : mp) {
        new2old[new_id] = old_rep;
    }

    tmp_cnf_file = "./simplified_" + cnf_instance_name_ + ".cnf";
    ofstream out(tmp_cnf_file);
    out << "p cnf " << mp.size() << " " << finalc.size() << "\n";
    for (auto& l : finalc) out << l << "\n";
    out.close();

    string cmd = string("./bin/coprocessor -verb=0 -enabled_cp3 -up -subsimp "
                        "-no-bve -no-bce -no-dense "
                        "-dimacs=") + reduced_cnf_file_path_ + " " + tmp_cnf_file;
    system(cmd.c_str());

    cnf_file_path_ = reduced_cnf_file_path_;
}

void CDCLCASampler::RecoverSolutions() {

    final_testcase_set_.reserve(testcase_set_.size());

    for (const auto& reduced_model : testcase_set_) {

        vector<int> full_model(original_num_var_, -1);

        for (int new_id = 0; new_id < num_var_; ++new_id) {
            int old_rep = new2old[new_id + 1];
            int val = reduced_model[new_id];

            full_model[old_rep - 1] = val;
        }

        for (int v = 1; v <= original_num_var_; ++v) {
            int idx = v - 1;
            if (full_model[idx] != -1)
                continue;

            int rep = dsu.find(v);
            int rep_idx = rep - 1;

            if (full_model[rep_idx] == -1)
                continue;

            int s = dsu.sign(v);
            full_model[idx] = (s == 1) ? full_model[rep_idx] : (1 - full_model[rep_idx]);
        }

        for (int lit : backbones) {
            int v = abs(lit);
            int idx = v - 1;

            full_model[idx] = (lit > 0) ? 1 : 0;
        }

        final_testcase_set_.push_back(move(full_model));
    }
}

void CDCLCASampler::InitCombinationOffset() {
    combination_offset.resize(num_var_);
    for (long i = 0; i < num_var_; ++i)
        combination_offset[i] = (2ll * num_var_ - i - 1) * i / 2;
}

void CDCLCASampler::clear_final(){
    Reduce_redundancy();
    already_t_wise.clear();
    for(int idx = 0; idx < testcase_set_.size(); idx++)
    {
        const vector<int>& testcase = testcase_set_[idx];
        for (int i = 0; i < num_var_ - 1; i++)
        {
            for (int j = i + 1; j < num_var_; j++)
            {
                long long index_tuple = Get2TupleMapIndex(i, testcase[i], j, testcase[j]);
                bool res = already_t_wise.set(index_tuple);
            }
        }
    }
    num_tuple_all_exact_ = num_tuple_all_possible_;
    cdcl_solver = new CDCLSolver::Solver;
    cdcl_solver->read_clauses(num_var_, clauses);
    choose_final_plain();
    delete cdcl_solver;
}

void CDCLCASampler::Update_already_t_wise(){
    int index_testcase = num_generated_testcase_ - 1;
    const vector<int>& testcase = testcase_set_[index_testcase];
    for (int i = 0; i < num_var_ - 1; i++)
    {
        for (int j = i + 1; j < num_var_; j++)
        {
            long long index_tuple = Get2TupleMapIndex(i, testcase[i], j, testcase[j]);
            bool res = already_t_wise.set(index_tuple);
        }
    }
    num_tuple_ = already_t_wise.count();
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

void CDCLCASampler::GenerateTestcase(){
    GenerateCandidateTestcaseSet();
    selected_candidate_index_ = SelectTestcaseFromCandidateSetByTupleNum();

    int testcase_index = num_generated_testcase_;
    vector<int> test_case_to_add_  = pq[selected_candidate_index_].first;
    testcase_set_.emplace_back(test_case_to_add_);
    selected_candidate_bitset_ = pq[selected_candidate_index_].second;
}

void CDCLCASampler::GenerateCandidateTestcaseSet(){
    auto prob_vec = get_prob_in();
    vector<pair<int, int>> sample_prob;
    sample_prob.reserve(prob_vec.size());
    for (const auto& p : prob_vec) {
        sample_prob.emplace_back(static_cast<int>(p.first), static_cast<int>(p.second));
    }

    for (int i = 0 ; i < candidate_set_size_; i++){
        cdcl_sampler->set_prob(sample_prob);
        cdcl_sampler->get_solution(candidate_testcase_set_[i]);
    }    
}

vector<pair<float, float> > CDCLCASampler::get_prob_in(){
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
    return sample_proportion;
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

void CDCLCASampler::Reduce_redundancy() {
    const int THRESHOLD = 200;
    const int N = (int)testcase_set_.size();
    if (N == 0) return;

    std::vector<int> tuple_owner(num_tuple_all_possible_, -1);
    std::vector<int> contribution(N, 0);

    for (int t = 0; t < N; ++t) {
        const auto& tc = testcase_set_[t];
        for (long i = 0; i < num_var_ - 1; ++i) {
            long v_i = tc[i];

            for (long j = i + 1; j < num_var_; ++j) {
                long v_j = tc[j];
                long long tid = Get2TupleMapIndex(i, v_i, j, v_j);
                int& owner = tuple_owner[tid];

                if (owner == -1) {
                    owner = t;
                    contribution[t]++;
                }
                else if (owner >= 0) {
                    contribution[owner]--;
                    owner = -2;
                }
            }
        }
    }

    std::vector<std::vector<int>> new_set;
    new_set.reserve(N);

    for (int t = 0; t < N; ++t) {
        if (contribution[t] >= THRESHOLD) {
            new_set.push_back(std::move(testcase_set_[t]));
        }
    }

    testcase_set_.swap(new_set);
}

void CDCLCASampler::remove_temp_files(){
    string cmd;
    int ret;
    cmd = "rm " +  tmp_cnf_file;
    ret = system(cmd.c_str());
}

void CDCLCASampler::find_uncovered_tuples(){
    uncovered_tuples.clear();
    for (int i = 0; i < num_var_ - 1; ++i) {
        for (int j = i + 1; j < num_var_; ++j) {
            int seen = 0;
            for (const auto& tc : testcase_set_) {
                int vi = tc[i] > 0 ? 1 : 0;
                int vj = tc[j] > 0 ? 1 : 0;
                int bit_index = (vi << 1) | vj;
                seen |= 1 << bit_index;
            }

            for (int k = 0; k < 4; ++k) {
                if (!(seen & (1 << k))) {
                    int vi = (k >> 1) & 1;
                    int vj = k & 1;

                    cdcl_solver->add_assumption(i, vi);
                    cdcl_solver->add_assumption(j, vj);
                    bool flag = cdcl_solver->solve();
                    cdcl_solver->clear_assumptions();

                    if (flag) {
                        uncovered_tuples.emplace_back(vector<int>{
                            vi ? (i + 1) : -(i + 1),
                            vj ? (j + 1) : -(j + 1)
                        });
                    } else {
                        --num_tuple_all_exact_;
                    }
                }
            }
        }
    }
}

void CDCLCASampler::Update2TupleInfo(){
    const vector<int>& testcase = testcase_set_.back();
    
    for (int i = 0; i < num_var_ - 1; i++)
    {
        for (int j = i + 1; j < num_var_; j++)
        {
            long long index_tuple = Get2TupleMapIndex(i, testcase[i], j, testcase[j]);
            bool res = already_t_wise.set(index_tuple);
            if (res){num_tuple_++;}
        }
    }
}

void CDCLCASampler::choose_final_plain()
{
    find_uncovered_tuples();
    shuffle(uncovered_tuples.begin(), uncovered_tuples.end(),
            mt19937(random_device{}()));

    CDCLSolver::Solver* solver = new CDCLSolver::Solver;
    solver->read_clauses(num_var_, clauses);

    vector<int8_t> assu(num_var_, -1);
    vector<int> assu_stack;
    vector<int> satisfied_indices;

    assu_stack.reserve(num_var_);
    satisfied_indices.reserve(uncovered_tuples.size());

    while (!uncovered_tuples.empty())
    {
        solver->clear_assumptions();
        fill(assu.begin(), assu.end(), -1);
        assu_stack.clear();
        satisfied_indices.clear();

        for (int j = 0; j < (int)uncovered_tuples.size(); ++j)
        {
            const auto& t = uncovered_tuples[j];

            int v1   = abs(t[0]) - 1;
            int val1 = (t[0] > 0);
            int v2   = abs(t[1]) - 1;
            int val2 = (t[1] > 0);

            bool conflict = false;
            int old_stack_size = assu_stack.size();
            int old_assu_size  = solver->assumptions.size();

            if (assu[v1] == -1) {
                assu[v1] = val1;
                assu_stack.push_back(v1);
                solver->add_assumption(v1, val1);
            } else if (assu[v1] != val1) {
                conflict = true;
            }

            if (!conflict) {
                if (assu[v2] == -1) {
                    assu[v2] = val2;
                    assu_stack.push_back(v2);
                    solver->add_assumption(v2, val2);
                } else if (assu[v2] != val2) {
                    conflict = true;
                }
            }

            if (conflict || !solver->solveLimited()) {
                while ((int)assu_stack.size() > old_stack_size) {
                    int v = assu_stack.back();
                    assu_stack.pop_back();
                    assu[v] = -1;
                }
                solver->assumptions.resize(old_assu_size);
                continue;
            }

            satisfied_indices.push_back(j);
        }

        solver->clear_assumptions();
        for (int v : assu_stack) {
            solver->add_assumption(v, assu[v]);
        }

        if (!solver->solve()) {
            continue;
        }

        vector<int> new_tc;
        solver->get_solution(new_tc);
        testcase_set_.emplace_back(new_tc);

        SimpleBitSet tmp(num_tuple_all_possible_);
        for (int i = 0; i < num_var_ - 1; ++i) {
            for (int j = i + 1; j < num_var_; ++j) {
                long long idx = Get2TupleMapIndex(
                    i, new_tc[i], j, new_tc[j]);
                tmp.set(idx);
            }
        }
        ++num_generated_testcase_;
        Update2TupleInfo();

        vector<char> covered(uncovered_tuples.size(), 0);
        for (int idx : satisfied_indices)
            covered[idx] = 1;

        vector<vector<int>> new_uncovered;
        new_uncovered.reserve(uncovered_tuples.size());

        for (int i = 0; i < (int)uncovered_tuples.size(); ++i) {
            if (!covered[i])
                new_uncovered.emplace_back(uncovered_tuples[i]);
        }

        uncovered_tuples.swap(new_uncovered);
    }

    delete solver;
}


