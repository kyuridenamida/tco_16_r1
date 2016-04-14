#include "classes_for_tco.cpp"
#include <cstdio>

const double TIME_LIMIT = 9.8;
namespace MyTimer{
	unsigned long long int cycle_per_sec = 2500000000;
	unsigned long long int beginCycle;
	unsigned long long int getCycle()
	{
	  unsigned int low, high;
	  __asm__ volatile ("rdtsc" : "=a" (low), "=d" (high));
	  return ((unsigned long long int)low) | ((unsigned long long int)high << 32);
	}
	double nowTime()
	{
	  return (double)(getCycle() - beginCycle) / cycle_per_sec;
	}
	void resetTimer(){
		beginCycle = getCycle();
	}
	bool isTLE(double limit){
		return nowTime() >= limit;
	}
};
using namespace MyTimer;

struct Configuration{
	bool visualize_mode;
	Configuration(){
		visualize_mode = false;
	}
} configuration;
class CutTheRoots {
public:
	
    vector<int> makeCuts(int NP, vector<int> points, vector<int> roots) {
		resetTimer();
		
		Problem problem(NP,points,roots);
		double current = -1e9;
		Answer res_ans;
		while( nowTime() < TIME_LIMIT ){
			//cerr << NaiveScoring::overall_score(problem,{}) << endl;
			ExtendedAnswer answer = greedy1(problem);
			return answer.to_vector();
			if( nowTime() < TIME_LIMIT ){
				if( answer.overall_score() > current ){
					res_ans = answer;
					current = answer.overall_score();
				}
			}
			break;
		}
		// cerr << NaiveScoring::overall_score(problem,answer) << endl;
        return res_ans.to_vector();
    }
	ExtendedAnswer greedy1(Problem &problem){
		vector<RGB> colors;
		for(int i = 0 ; i < problem.trees.size() ; i++) colors.push_back(RGB::random());
		int N = problem.trees.size();

		// cerr << answer.overall_score() << endl;
		
	
		// if( configuration.visualize_mode ) answer.draw(problem,colors);
		vector<pair<double,pair<int,int>>> pairs;
		for(int i = 0 ; i < N ; i++){
			for(int j = i+1; j < N ; j++){
				P p1 = problem.trees[i].position[problem.trees[i].root];
				P p2 = problem.trees[j].position[problem.trees[j].root];
				pairs.push_back({norm(p1-p2),{i,j}});
			}
		}
		
		sort(pairs.begin(),pairs.end());
		
		
		ExtendedAnswer final_answer(&problem);
		double currentScore = -1;
		
		while(true){
			ExtendedAnswer answer(&problem);
			random_shuffle(final_answer.lines.begin(),final_answer.lines.end());
			for(int i = 0 ; i < final_answer.lines.size() * 8 / 10 ; i++){
				answer.add_line(final_answer.lines[i]);
			}
			
			for(int li = 0 ; li < pairs.size(); li++){
				int i = pairs[li].second.first;
				int j = pairs[li].second.second;
				P p1 = problem.trees[i].position[problem.trees[i].root];
				P p2 = problem.trees[j].position[problem.trees[j].root];
				bool separated = false;
				for( auto l : answer.lines ){
					if( GeomUtils::is_separating(l,p1,p2) ){
						separated = true;
						break;
					}
				}
				
				if( !separated ){
					vector<pair<double,L>> cand;
					double cur = answer.overall_score();
					for(int m = 0 ; m < 80000 / problem.trees.size() ; m++){
						if( nowTime() > TIME_LIMIT ){
							
							return final_answer;
						}
						double A = 1. * xor64() / (unsigned int)(-1);
						double B = 1. * xor64() / (unsigned int)(-1);
						P mp = p1 + (p2-p1) * A;
					// cerr << cur << " " << NaiveScoring::overall_score_fast(problem,answer) << "|" <<  NaiveScoring::overall_score(problem,answer) << endl;
						P vec = (p1-p2) * exp(P(0,PI*B));
						L l = L(mp,mp+vec);
						if( l[0] != P(-1,-1) and GeomUtils::is_separating(l,p1,p2) ){
							double loss = cur - answer.overall_score(l);
							//cerr << loss << endl;
							cand.emplace_back(loss,l);
						}
						// cerr << endl;
					}
					// cerr << endl;
					sort(cand.begin(),cand.end(),[&](const pair<double,L> &a,const pair<double,L> &b){
						return a.first < b.first;
						
					});
					for(int i = 0 ; i < cand.size() ; i++){
						L fix_l = GeomUtils::convert_to_integer_line(cand[i].second,p1,p2);
						if( fix_l[0] != null_point and GeomUtils::is_separating(fix_l,p1,p2) ){
								// double loss = answer.overall_score() - answer.overall_score(fix_l);
								// cerr << loss -  cand[0].first << endl;
								answer.add_line(fix_l);
								
								break;
						}
					}
					
					if( configuration.visualize_mode ) answer.draw(problem,colors);
					// cerr << answer.overall_score() << " "  << NaiveScoring::overall_score(problem,answer) << endl;
				}
			}
			
			// cerr << answer.overall_score() << endl;
			for(int i = 0 ; i < answer.lines.size() ; ){
				vector<L> ls;
				for(int j = 0 ; j < answer.lines.size() ; j++)
					if( i != j ) ls.push_back(answer.lines[j]);
				bool f = true;
				for(int j = 0 ; j < N ; j++){
					for(int k = j+1 ; k < N ; k++){
						P p1 = problem.trees[j].position[problem.trees[j].root];
						P p2 = problem.trees[k].position[problem.trees[k].root];
						bool separated = false;
						for( auto l : ls ){
							if( GeomUtils::is_separating(l,p1,p2) ){
								separated = true;
								break;
							}
						}
						if( !separated ) {
							f = false;
							break;
						}
					}
					if( !f ) break;
				}
				if( f ) {
					answer.lines.erase(answer.lines.begin()+i);
					
					//cerr << "OK" << endl;
				}else i++;
			}
			if( currentScore < answer.overall_score() ){
				ExtendedAnswer hogehoge(&problem);
				for(int i = 0 ; i < answer.lines.size() ; i++){
					hogehoge.add_line(answer.lines[i]);
				}
				final_answer = hogehoge;
				currentScore = hogehoge.overall_score();
				cerr << currentScore << endl;
			}
			cerr << nowTime() << endl;
			random_shuffle(pairs.begin(),pairs.end());
		}
	}
};
#ifdef LOCAL
template<class T> void getVector(vector<T>& v) {
    for (int i = 0; i < v.size(); ++i)
        cin >> v[i];
}



int main(int argc,char *argv[]) {
	if( argc >= 2 && string(argv[1]) == "-vis" ){
		cerr << "visualize mode on" << endl;
		configuration.visualize_mode = true;
	}
    int NP;
    cin >> NP;
    int Npoints;
    cin >> Npoints;
    vector<int> points(Npoints);
    getVector(points);

    int Nroots;
    cin >> Nroots;
    vector<int> roots(Nroots);
    getVector(roots);

    CutTheRoots cr;
    vector<int> ret = cr.makeCuts(NP, points, roots);
	
	if( !configuration.visualize_mode ){
		cout << ret.size() << endl;
		for (int i = 0; i < ret.size(); ++i) {
			cout << ret[i] << endl;
		}
		cout.flush();
	}
}
#endif
