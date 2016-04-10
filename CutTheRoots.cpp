#include "classes_for_tco.cpp"
#include <cstdio>
struct Configuration{
	bool visualize_mode;
	Configuration(){
		visualize_mode = false;
	}
} configuration;
class CutTheRoots {
public:
	
    vector<int> makeCuts(int NP, vector<int> points, vector<int> roots) {
		srand(time(NULL));
		//cout << intersectLS(L(P(0,0),P(100,0)),L(P(1,0),P(2,0))) << endl;
		Problem problem(NP,points,roots);
		//cerr << NaiveScoring::overall_score(problem,{}) << endl;
		Answer answer = greedy1(problem);
        return answer.to_vector();
    }
	Answer greedy1(const Problem &problem){
		vector<RGB> colors;
		for(int i = 0 ; i < problem.trees.size() ; i++) colors.push_back(RGB::random());
		int N = problem.trees.size();
		ExtendedAnswer answer;
		answer.init(problem);
		if( configuration.visualize_mode ) answer.draw(problem,colors);
		
		vector<pair<double,pair<int,int>>> pairs;
		for(int i = 0 ; i < N ; i++){
			for(int j = i+1; j < N ; j++){
				P p1 = problem.trees[i].position[problem.trees[i].root];
				P p2 = problem.trees[j].position[problem.trees[j].root];
				pairs.push_back({norm(p1-p2),{i,j}});
			}
		}
		
		sort(pairs.begin(),pairs.end());
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
				double cur = NaiveScoring::overall_score_fast_nonline(problem,answer);
				for(int m = 1 ; m < 10 ; m++){
					double lll = INF;
					P mp = p1 + (p2-p1) * (double)(m/10.);
					// cerr << cur << " " << NaiveScoring::overall_score_fast(problem,answer) << "|" <<  NaiveScoring::overall_score(problem,answer) << endl;
					for(int k = 0 ; k < 19 ; k++){
						P vec = (p1-p2) * exp(P(0,PI*k/19));
						L l = L(mp,mp+vec);
						if( l[0] != P(-1,-1) and GeomUtils::is_separating(l,p1,p2) ){						
							ExtendedAnswer copied_answer = answer;
							double loss = cur - NaiveScoring::overall_score_fast_differ_ver(problem,copied_answer,l);
							// printf("%5.0lf ",loss);
							lll = min(lll,loss);
							//cerr << loss << endl;
							cand.push_back({loss,l});
						}
					}
					// cerr << endl;
				}
				// cerr << endl;
				sort(cand.begin(),cand.end(),[&](const pair<double,L> &a,const pair<double,L> &b){
					return a.first < b.first;
					
				});
				for(int i = 0 ; i < cand.size() ; i++){
					L fix_l = GeomUtils::convert_to_integer_line(cand[i].second,p1,p2);
					if( fix_l[0] != P(-1,-1) and GeomUtils::is_separating(fix_l,p1,p2) ){
							NaiveScoring::overall_score_fast_differ_ver(problem,answer,fix_l);
							break;
					}
				}
				
				if( configuration.visualize_mode ) answer.draw(problem,colors);
			}
		}
		return answer;
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
