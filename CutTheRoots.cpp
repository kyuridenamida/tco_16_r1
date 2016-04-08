#include "classes_for_tco.cpp"

class CutTheRoots {
public:
	
    vector<int> makeCuts(int NP, vector<int> points, vector<int> roots) {
		//cout << intersectLS(L(P(0,0),P(100,0)),L(P(1,0),P(2,0))) << endl;
		Problem problem(NP,points,roots);
		//cerr << NaiveScoring::overall_score(problem,{}) << endl;
		
		Answer answer = greedy1(problem);
        return answer.to_vector();
    }
	Answer greedy1(const Problem &problem){
		int N = problem.trees.size();
		ExtendedAnswer answer;
		answer.init(problem);
		
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
				
				P mp = 0.5 * (p1 + p2);
				double cur = NaiveScoring::overall_score_fast_nonline(problem,answer);
				// cerr << cur << " " << NaiveScoring::overall_score_fast(problem,answer) << "|" <<  NaiveScoring::overall_score(problem,answer) << endl;
				vector<pair<double,ExtendedAnswer>> cand;
				for(int k = 0 ; k < 69 ; k++){
					P vec = (p1-p2) * exp(P(0,2*PI*k/69));
					L l = L(mp,mp+vec);
					L fix_l = GeomUtils::convert_to_integer_line(l,p1,p2);
					if( fix_l[0] != P(-1,-1) and GeomUtils::is_separating(fix_l,p1,p2) ){						
						ExtendedAnswer copied_answer = answer;
						double loss = cur - NaiveScoring::overall_score_fast_differ_ver(problem,copied_answer,fix_l);
						cand.push_back({loss,copied_answer});
					}
				}
				sort(cand.begin(),cand.end(),[&](const pair<double,ExtendedAnswer> &a,const pair<double,ExtendedAnswer> &b){
					return a.first < b.first;
					
				});
				// answer.add_line(cand[0].second);
				// cerr << answer.lines.size() << endl;	
				answer = cand[0].second;
				// cerr << fix_l[0] << " " << fix_l[1] << endl;
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


int main() {
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
	
    cout << ret.size() << endl;
    for (int i = 0; i < ret.size(); ++i) {
        cout << ret[i] << endl;
    }
    cout.flush();
}
#endif
