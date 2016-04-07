#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stack>
#include <cassert>
#include <sstream>
#include <vector>
using namespace std;
#include "geom.cpp"

struct Edge{
	int src,dst;
	Edge(int src,int dst) : src(src), dst(dst) {}
};


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
	Answer(){}
	Answer(vector<int> ps){
		assert(ps.size() % 2 == 0);
		for(int i = 0 ; i < ps.size() ; i += 4){
			lines.push_back(L(P(ps[i],ps[i+1]),P(ps[i+2],ps[i+3])));
		}
	}
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
	static L convert_to_integer_line(L l,P p1,P p2){
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
		
		for(int i = 0 ; i < ps.size() ; i++){
			for(int j = i+1 ; j < ps.size() ; j++){
				if( abs(ps[i].second-ps[j].second) < EPS ) continue;
				if( GeomUtils::is_separating(L(ps[i].second,ps[j].second),p1,p2) ){
					return L(ps[i].second,ps[j].second);
				}
			}
		}
		assert(false and "ps!!");
		
	}
};
