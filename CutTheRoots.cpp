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
		Answer answer;
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
				P vec = (p1-p2) * P(0,1);
				L l = L(mp,mp+vec);
				// cerr << l[0] << " " << l[1] << " -> ";
				L fix_l = GeomUtils::convert_to_integer_line(l,p1,p2);
				if( GeomUtils::is_separating(fix_l,p1,p2) ){
					answer.add_line(fix_l);
				}else{
					cerr << "oops" << endl;
				}
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
