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
	double score_of_tree(const Tree &tree,const Answer &answer){
		return inner_dfs_score_of_tree(tree,answer,tree.root);
	}
	double overall_score(const Problem &problem,const Answer &answer){
		double sum = 0;
		for( const auto &tree : problem.trees ){
			sum += score_of_tree(tree,answer);
		}
		return sum;
	}
private:
	double inner_dfs_score_of_tree(const Tree &tree,const Answer &answer,int x){
		double sum = 0;
		for( auto c : tree.child[x] ){
			// ここ前計算で114514倍くらい速くなると思う
			double cut_dist = INF;
			for( const auto& line : answer.lines ){
				if( intersectLS(line,L(tree.position[x],tree.position[c])) ){
					P cp = crosspoint(line,L(tree.position[x],tree.position[c]));
					assert(abs(cp-tree.position[x]) > EPS and abs(cp-tree.position[c]) > EPS);
					cut_dist = min(cut_dist,abs(cp-tree.position[x]));
				}
			}
			if( cut_dist != INF ){
				sum += cut_dist; 
			}else{
				sum += inner_dfs_score_of_tree(tree,answer,c);
				sum += abs(tree.position[x]-tree.position[c]);
			}
		}
		return sum;
	}

};

class CutTheRoots {
public:
	
    vector<int> makeCuts(int NP, vector<int> points, vector<int> roots) {
		Problem problem(NP,points,roots);
		for( auto tree : problem.trees ){
			cout << tree.position[tree.root] << endl;
		}
        // first NP points give coordinates of the bases of the plants, written as x, y
        // get x-coordinates of the bases, sort them, make cuts along y-axis to separate x-coordinates
        vector<int> xs(NP);
        for (int i = 0; i < NP; ++i) {
            xs[i] = points[2 * i];
        }
        sort(xs.begin(), xs.end());
        vector<int> ret(4 * (NP - 1));
        for (int i = 0; i < NP - 1; ++i) {
            int x = (xs[i] + xs[i + 1]) / 2;
            ret[4 * i] = x;
            ret[4 * i + 1] = 1;
            ret[4 * i + 2] = x;
            ret[4 * i + 3] = 2;
        }
        return ret;
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
