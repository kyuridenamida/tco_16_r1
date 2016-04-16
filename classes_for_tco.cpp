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

namespace Utils{
	void reflesh_max(vector<double> &a,const vector<double> &b){
		for(int i = 0 ; i < a.size() ; i++)
			a[i] = max(a[i],b[i]);
	}
	int to_int(double x){
		return x + 0.5;
	}
};
using namespace Utils;
const P null_point = P(-1,-1);
const L null_line = L(null_point,null_point);

uint32_t xor64(void) {
  static uint64_t x = 88172645463325252ULL;
  x = x ^ (x << 13); x = x ^ (x >> 7);
  return x = x ^ (x << 17);
}

struct Edge{
	int src,dst;
	Edge(int src,int dst) : src(src), dst(dst) {}
};


class Tree{
public:
	int id;
	int n;
	vector<Edge> es;
	vector<P> position;
	int root;
	vector<Circle> mec;
	vector<P> convex_polygon;
	vector<vector<int>> child;
	vector<vector<L>> child_line;
	vector<int> parent;
	vector<int> depth;
	vector<double> length_between_parent;
	vector<double> tot_sum_of_tree;
	Tree(int id,int n,const vector<Edge> es,vector<P> position,int root) : id(id), n(n), es(es), position(position), root(root){
		//convex_polygon = convex_hull(position);
		
		
		
		parent = vector<int>(n,-1);
		depth = vector<int>(n);
		child = vector<vector<int>>(n);
		child_line = vector<vector<L>>(n);
		tot_sum_of_tree = length_between_parent = vector<double>(n,0);
		mec =  vector<Circle>(n,Circle(P(0,0),-1));
		
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
					
					length_between_parent[e] = abs(position[e]-position[x]);
					S.push(array<int,3>{e,x,d+1});
					child[x].push_back(e);
					child_line[x].push_back(L(position[x],position[e]));
				}
			}
		}
		init_dfs(root);
	}
	vector<P> init_dfs(int x){
		tot_sum_of_tree[x] = 0.0;
		vector<P> ps;
		ps.push_back(position[x]);
		
		for( auto c : child[x] ){
			vector<P> X = init_dfs(c);
			tot_sum_of_tree[x] += tot_sum_of_tree[c];// + length_between_parent[c];
			ps.insert(ps.end(),X.begin(),X.end());
		}
	
		
		ps = convex_hull2(ps);
		mec[x] = Circle::minEnclosingCircle(ps);
		tot_sum_of_tree[x] += length_between_parent[x];
		// cerr << tot_sum_of_tree[x] << endl;
		//cerr << ps.size() << endl;
		return ps;
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
				trees.push_back(Tree(trees.size(),position.size(),es,position,0));
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


class Answer{
public:
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
			vs.push_back(line.a.x+0.5);
			vs.push_back(line.a.y+0.5);
			vs.push_back(line.b.x+0.5);
			vs.push_back(line.b.y+0.5);
		}
		return vs;
	}
};

class RGB{
public:
	int r,g,b;
	RGB(int r=0,int g=0,int b=0) : r(r), g(g), b(b) {}
	static RGB random(){
		RGB res;
		res.r = rand() % 256;
		res.g = rand() % 256;
		res.b = rand() % 256;
		return res;
	}
};



class ExtendedAnswer : public Answer{
public:
	Problem* problem;
	// pre
	vector<G> convex_polygon;
	vector<double> original_area;

	//2387500
	vector<Circle> mecs;
	

	ExtendedAnswer(Problem *problem) : problem(problem){
		double all_total = 0;
		for(const auto &tree : problem->trees ){
			all_total += tree.tot_sum_of_tree[0];
			// cerr << tree.tot_sum_of_tree[0] << "<" << endl;
			original_area.push_back(area2(tree.convex_polygon));
			convex_polygon.push_back(tree.convex_polygon);
			mecs.push_back(Circle::minEnclosingCircle(tree.convex_polygon));
			MEMO_score_of_tree.push_back(tree.tot_sum_of_tree[0]);
			current_weight.push_back(vector<double>(tree.n,-1.0));
			already_cut.push_back(vector<bool>(tree.n,false));
			inner_init_dfs_score_of_tree(tree.root,tree);
			
		}
		// cerr << all_total << "<" << endl;
		MEMO_overall_score = all_total;
		
	}
	void add_line(const L &line){
		lines.push_back(line);
		// overall_score_rough(line,true);
		overall_score(line,true);
	}
	
	vector<double> MEMO_score_of_tree;

	double score_of_tree(const Tree &tree,const L &l = null_line,bool reflesh=true){
		if( l == null_line )
			return MEMO_score_of_tree[tree.id];
		double res = inner_dfs_score_of_tree(tree.root,tree,l,reflesh);
		if( reflesh ){
			MEMO_score_of_tree[tree.id] = res;
		}
		return res;
	}


	double MEMO_overall_score;
	double overall_score(const L &l = null_line,bool reflesh=false){
		if( l == null_line )
			return MEMO_overall_score;

		double sum = 0;
		for( const auto &tree : problem->trees ){
			sum += score_of_tree(tree,l,reflesh);
		}
		if( reflesh ){
			MEMO_overall_score = sum;
		}
		return sum;
	}

	vector< vector<double> > current_weight;
	vector< vector<bool> > already_cut;
	
	double inner_init_dfs_score_of_tree(int x,const Tree &tree){
		double sum = 0;
		for( auto c : tree.child[x] ){
			sum += inner_init_dfs_score_of_tree(c,tree);
		}
		return current_weight[tree.id][x] = sum + tree.length_between_parent[x];
	}
	
	double inner_dfs_score_of_tree(int x,const Tree &tree,const L &line,bool reflesh){
		
		if( distanceLP_check(line,tree.mec[x].p,tree.mec[x].r) ){
			return current_weight[tree.id][x];
		}
		if( already_cut[tree.id][x] ) return current_weight[tree.id][x];


		double sum = 0;
		for(int i = 0 ; i < tree.child[x].size() ; i++){
			const int &c = tree.child[x][i];
			const L &cline = tree.child_line[x][i];
			
			if( intersectLS(line,cline) ){
				P cp = crosspoint(line,cline);
				double cut_dist = min( current_weight[tree.id][c], abs(cp-tree.position[x]) );
				if( reflesh ){
					already_cut[tree.id][c] = true;
					current_weight[tree.id][c] = cut_dist ; // 枝がさらにカットされたとき対策
				}
				sum += cut_dist;
			}else{
				sum += inner_dfs_score_of_tree(c,tree,line,reflesh);			
			}
		}
		
		double res = sum + tree.length_between_parent[x];
		if( reflesh ){
			current_weight[tree.id][x] = res;
		}
		return res;
	}
	vector<double> get_ratio(){
		vector<double> res;
		for( const auto &tree : problem->trees ){
			res.push_back(score_of_tree(tree) / tree.tot_sum_of_tree[0]);
			//cerr << tree.id << "|" << score_of_tree(tree) / tree.tot_sum_of_tree[0] << "|" << score_of_tree(tree)  << "/" <<  tree.tot_sum_of_tree[0] << endl;
		}
		//cerr << "--------" << endl;
		return res;
	}

	void inner_draw_tree(int x,bool dead,const Tree &tree,const RGB &alive){
		if( already_cut[tree.id][x] ) dead = true;
		
		
		for( auto c : tree.child[x] ){
			cout << "S " << alive.r << " " << alive.g << " " << alive.b << " " << to_int(tree.position[x].x) << " " << to_int(tree.position[x].y) << " " << to_int(tree.position[c].x) << " " << to_int(tree.position[c].y) << endl;
			inner_draw_tree(c,dead,tree,alive);
		}
		
		return;
	}
	void draw_sol(){
		cout << "C " << 128 << " " << 128 << " " << 128 << " " << 512 << " " << 512 << " " << 512 << endl;	
		static vector<RGB> cs;
		while( cs.size() < 200 ) cs.push_back(RGB::random());
		
		for( auto& tree : problem->trees ){
			inner_draw_tree(0,false,tree,cs[tree.id]);
		}
		cout << "END" << endl;
	}

};



void inner_dfs_remaked_tree(int x,int bit,Tree &tree,const vector<ExtendedAnswer> &answers){
	for(int i = 0 ; i < answers.size() ; i++){
		if( answers[i].already_cut[tree.id][x] ) bit |= 1 << i;	
	}
	// cerr << bit << endl;
	if( bit == (1<<answers.size()) - 1 ){
		tree.child[x].clear();
		return;
	}
	
	for( auto c : tree.child[x] ){
		inner_dfs_remaked_tree(c,bit,tree,answers);
	}
	
	return;
}
void remake_tree(Tree &tree, const vector<ExtendedAnswer> &answers){
	assert( answers.size() <= 32 );
	inner_dfs_remaked_tree(0,0,tree,answers);
	tree.init_dfs(0);
	
}
Problem* remake_trees(Problem *problem,const vector<ExtendedAnswer> &answers){
	Problem* new_prob = new Problem(*problem);
	for(int i = 0 ; i < new_prob->trees.size() ; i++){
		remake_tree(new_prob->trees[i],answers);
	}
	return new_prob;
}


class GeomUtils{
public:
	static bool is_separating(L l,P p1,P p2){
		int r1 = ccw(l.a,l.b,p1);
		int r2 = ccw(l.a,l.b,p2);
		return abs(r1) == 1 && abs(r2) == 1 && r1 != r2;
	}
	static L convert_to_integer_line(L l,P p1,P p2){
		// cout << l.b << " " << l.a << endl;
		P vec = (l.b - l.a) / abs(l.b-l.a);
		// cout << vec << endl;
		if( abs(vec.x) < EPS ) {
			int x = (l.a.x+0.5);
			return L(P(x,0),P(x,1));
		}
		if( abs(vec.y) < EPS ){
			int y = (l.a.y+0.5);
			return L(P(0,y),P(1,y));
		}
		P vx = vec / vec.x;
		//cout << vx << endl;
		vector< pair<double,P> > ps;
		for(int i = 0 ; i <= 1024 ; i++){
			P t = l.a + (i-l.a.x) * vx;
			if( -EPS <= t.x && t.x <= 1024 + EPS && 
				-EPS <= t.y && t.y <= 1024 + EPS ){
				int X = t.x + 0.5;
				int Y = t.y + 0.5;
				ps.push_back({abs(t.x-X)+abs(t.y-Y),P(X,Y)});
			}
		}
		P vy = vec / vec.y;
		// cout << l.a << " " << vy << endl;
		for(int i = 0 ; i <= 1024 ; i++){
			P t = l.a + (i-l.a.y) * vy;
			if( -EPS <= t.x && t.x <= 1024 + EPS && 
				-EPS <= t.y && t.y <= 1024 + EPS ){
				int X = t.x + 0.5;
				int Y = t.y + 0.5;
				ps.push_back({abs(t.x-X)+abs(t.y-Y),P(X,Y)});
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
		return L(P(-1,-1),P(-1,-1));
	}
};
