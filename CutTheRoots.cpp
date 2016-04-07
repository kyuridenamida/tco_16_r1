#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stack>
#include <cassert>
#include <sstream>
#include <vector>
#include "geom.cpp"
using namespace std;


struct Edge{
	int src,dst;
	Edge(int src,int dst) : src(src), dst(dst) {}
};


//直線


class Tree{
public:
	int n,root;
	vector<P> position;
	vector<vector<int>> child;
	vector<int> parent;
	vector<int> depth;
	vector<double> tot_sum_of_tree;
	vector<Edge> es;
	Tree(int n,const vector<Edge> &es,vector<P> position,int root) : n(n), position(position), es(es), root(root){
		parent = vector<int>(n,-1);
		depth = vector<int>(n);
		child = vector<vector<int>>(n);
		
		vector<vector<int>> g(n);
		for( auto e : es ){
			g[e.src].push_back(e.dst);
			g[e.dst].push_back(e.src);
		}
		
		stack< array<int,3> > S;
		S.push(array<int,3>{root,-1,0});
		while( S.size() ){
			int x = S.top()[0];
			int p = S.top()[1];
			int d = S.top()[2];
			S.pop();
			parent[x] = p;
			depth[x] = d;
			for( auto e : g[x] ){
				if( e != p ){
					S.push(array<int,3>{e,x,d+1});
					child[x].push_back(e);
				}
			}
		}
		
		tot_sum_of_tree = vector<double>(n,0);
		init_dfs(root);
	}
private:
	double init_dfs(int x){
		for( auto c : child[x] ){
			init_dfs(c);
			tot_sum_of_tree[x] += tot_sum_of_tree[c] + abs(position[x]-position[c]);
		}
	}
};


struct Answer{
	vector<L> lines; 
	void add_line(const L &line){
		lines.push_back(line);
	}
	
	vector<int> to_vector(){
		vector<int> vs;
		for( auto line : lines ){
			vs.push_back(line[0].real()+0.5);
			vs.push_back(line[0].imag()+0.5);
			vs.push_back(line[1].real()+0.5);
			vs.push_back(line[1].imag()+0.5);
		}
		return vs;
	}
};
class Problem{
public:
	vector<Tree> trees;
	Problem(int NP,vector<int> points,vector<int> roots){
		int n = points.size() / 2;
		
		map<P,vector<P>> graph;
		map<int,P> id_to_pos;	
		for(int i = 0 ; i < n ; i++)
			id_to_pos[i] = P(points[2*i],points[2*i+1]);
		
		vector<vector<int>> g(id_to_pos.size());
		for(int i = 0 ; i < roots.size() ; i += 2 ){
			g[roots[i]].push_back(roots[i+1]);
			g[roots[i+1]].push_back(roots[i]);
		}
		vector<bool> done(n,false);
		for(int root = 0 ; root < n ; root++){
			if( !done[root] ){
				map<int,int> relabel;
				vector<Edge> es;
				inner_listup_component(root,g,done,relabel,es);
				vector<P> position(relabel.size());
				for( auto r : relabel )
					position[r.second] = id_to_pos[r.first];
				trees.push_back(Tree(points.size(),es,position,0));
			}
		}
		assert( trees.size() == NP );
	}
private:
	void inner_listup_component(int x,const vector<vector<int>> &g,vector<bool> &done, map<int,int> &relabel,vector<Edge> &es){
		assert( !done[x] );
		done[x] = true;
		int k = relabel.size();
		relabel[x] = k;
		for( auto to : g[x] ){
			if( !done[to] ){
				inner_listup_component(to,g,done,relabel,es);
				es.push_back(Edge(relabel[x],relabel[to]));
			}
		}
		
		
	
	}
};

class NaiveScoring{
public:
	static double score_of_tree(const Tree &tree,const Answer &answer){
		return inner_dfs_score_of_tree(tree,answer,tree.root);
	}
	static double overall_score(const Problem &problem,const Answer &answer){
		double sum = 0;
		for( const auto &tree : problem.trees ){
			sum += score_of_tree(tree,answer);
		}
		return sum;
	}
private:
	static double inner_dfs_score_of_tree(const Tree &tree,const Answer &answer,int x){
		double sum = 0;
		for( auto c : tree.child[x] ){
			// ここ前計算で114514倍くらい速くなると思う
			double cut_dist = INF;
			for( const auto& line : answer.lines ){
				//assert(!(intersectLP(line,tree.position[x]) and !intersectLP(line,tree.position[c])));
				if( intersectLS(line,L(tree.position[x],tree.position[c])) ){
					P cp = crosspoint(line,L(tree.position[x],tree.position[c]));
					//assert(abs(cp-tree.position[x]) > EPS and abs(cp-tree.position[c]) > EPS);
					cut_dist = min(cut_dist,abs(cp-tree.position[x]));
				}
			}
			if( cut_dist != INF ){
				// cerr << cut_dist << " " << abs(tree.position[x]-tree.position[c]) << endl;
				sum += cut_dist; 
			}else{
				sum += inner_dfs_score_of_tree(tree,answer,c);
				sum += abs(tree.position[x]-tree.position[c]);
			}
		}
		return sum;
	}

};

class GeomUtils{
public:
	static bool is_separating(L l,P p1,P p2){
		int r1 = ccw(l[0],l[1],p1);
		int r2 = ccw(l[0],l[1],p2);
		return abs(r1) == 1 && abs(r2) == 1 && r1 != r2;
	}
	static L convert_to_integer_line(L l){
		// cout << l[1] << " " << l[0] << endl;
		P vec = (l[1] - l[0]) / abs(l[1]-l[0]);
		// cout << vec << endl;
		if( abs(vec.real()) < EPS ) {
			int x = (l[0].real()+0.5);
			return L(P(x,0),P(x,1));
		}
		if( abs(vec.imag()) < EPS ){
			int y = (l[0].imag()+0.5);
			return L(P(0,y),P(1,y));
		}
		P vx = vec / vec.real();
		//cout << vx << endl;
		vector< pair<double,P> > ps;
		for(int i = 0 ; i <= 1024 ; i++){
			P t = l[0] + (i-l[0].real()) * vx;
			if( -EPS <= t.real() && t.real() <= 1024 + EPS && 
				-EPS <= t.imag() && t.imag() <= 1024 + EPS ){
				int X = t.real() + 0.5;
				int Y = t.imag() + 0.5;
				ps.push_back({abs(t.real()-X)+abs(t.imag()-Y),P(X,Y)});
			}
		}
		P vy = vec / vec.imag();
		// cout << l[0] << " " << vy << endl;
		for(int i = 0 ; i <= 1024 ; i++){
			P t = l[0] + (i-l[0].imag()) * vy;
			if( -EPS <= t.real() && t.real() <= 1024 + EPS && 
				-EPS <= t.imag() && t.imag() <= 1024 + EPS ){
				int X = t.real() + 0.5;
				int Y = t.imag() + 0.5;
				ps.push_back({abs(t.real()-X)+abs(t.imag()-Y),P(X,Y)});
			}
		}
		sort(ps.begin(),ps.end());
		for(int i = 1 ; i < ps.size() ; i++){
			// cout << ps[i].second << "(" << vx << " " << " " << vy << endl;
			if( abs(ps[i].second-ps[0].second) > EPS ){
				//cout << ps[0].second << " " << ps[i].second << "|" << ps[i].first << endl;
				return L(ps[0].second,ps[i].second);
			}
		}
		assert(false and "ps!!");
		
	}
};

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
				L fix_l = GeomUtils::convert_to_integer_line(l);
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
