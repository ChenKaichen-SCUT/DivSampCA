#include "cdclcasampler.h"
#include <signal.h>

using namespace std;

void HandleInterrupt(int sig) {
	cout << "c" << endl;
	cout << "c caught signal... exiting" << endl;
	exit(-1);
}

void SetupSignalHandler() {
	signal(SIGTERM, HandleInterrupt);
	signal(SIGINT, HandleInterrupt);
	signal(SIGQUIT, HandleInterrupt);
	signal(SIGKILL, HandleInterrupt);
}

struct Argument{
	int seed;
	int transfer_size;
	int candidate_set_size;
	string input_cnf_path;
	string output_testcase_path;

	bool flag_input_cnf_path;
	bool flag_output_testcase_path;
	bool flag_candidate_set_size;
	bool flag_transfer_size;

};

bool ParseArgument(int argc, char **argv, Argument &argu){	
	argu.seed = 1;

	argu.flag_input_cnf_path = false;
	argu.flag_output_testcase_path = false;
	argu.flag_candidate_set_size = false;
	argu.flag_transfer_size = false;
	
	if (argc < 2) return false;

	for (int i = 1; i < argc; ++i){
        string arg = argv[i];

        if (arg == "-input_cnf_path"){
            i++;
            if (i >= argc) return false;
            argu.input_cnf_path = argv[i];
            argu.flag_input_cnf_path = true;
        }
        else if (arg == "-output_testcase_path"){
            i++;
            if (i >= argc) return false;
            argu.output_testcase_path = argv[i];
            argu.flag_output_testcase_path = true;
        }
        else if (arg == "-seed"){
            i++;
            if (i >= argc) return false;
            argu.seed = stoi(argv[i]); 
        }
        else if (arg == "-c"){
            i++;
            if (i >= argc) return false;
            argu.candidate_set_size = stoi(argv[i]);
            argu.flag_candidate_set_size = true;
        }
        else if (arg == "-t"){
            i++;
            if (i >= argc) return false;
            argu.transfer_size = stoi(argv[i]); 
            argu.flag_transfer_size = true;
        }
        else{
            return false;
        }
    }

	if(argu.flag_input_cnf_path) return true;
	else return false;
}

int main(int argc, char **argv)
{
	SetupSignalHandler();

	Argument argu;

	if (!ParseArgument(argc, argv, argu)){
		cout << "c Argument Error!" << endl;
		return -1;
	}

	CDCLCASampler ca_sampler(argu.input_cnf_path, argu.seed);

	if (argu.flag_output_testcase_path)
		ca_sampler.SetTestcaseSetSavePath(argu.output_testcase_path);
	if (argu.flag_candidate_set_size)
		ca_sampler.SetCandidateSetSize(argu.candidate_set_size);
	if (argu.flag_transfer_size)
		ca_sampler.SetTransferSize(argu.transfer_size);		
	ca_sampler.GenerateCoveringArray();

	return 0;
}
