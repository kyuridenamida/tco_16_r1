#include "classes_for_tco.cpp"
#include <fstream>
#include <cstdio>

#ifdef LOCAL


int main(int argc,char *argv[]) {
	if( argc != 3 ){
		cerr << "Something is wrong! ./exec [input] [output]" << endl;
		return 1;
	}
	
	ifstream input(argv[1]);
	ifstream output(argv[2]);
    int NP;
    input >> NP;
    int Npoints;
    input >> Npoints;
    vector<int> points(Npoints);
    for(int i = 0 ; i < points.size() ; i++) input >> points[i];
	
    int Nroots;
    input >> Nroots;
    vector<int> roots(Nroots);
	for(int i = 0 ; i < roots.size() ; i++) input >> roots[i];

	Problem problem(NP,points,roots);
	
	
	
	
	int ansN;
	output >> ansN;
	vector<int> ans(ansN);
	for(int i = 0 ; i < ansN ; i++)
		output >> ans[i];
	Answer answer = Answer(ans);
	printf("%.10lf\n",NaiveScoring::overall_score(problem,answer)/NaiveScoring::overall_score(problem,{})*1000000);
    cout.flush();
}
#endif
