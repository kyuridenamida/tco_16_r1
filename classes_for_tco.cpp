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

const P null_point = P(-1,-1);
const L null_line = L(null_point,null_point);


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
	vector<int> parent;
	vector<int> depth;
	vector<double> length_between_parent;
	vector<double> tot_sum_of_tree;
	Tree(int id,int n,const vector<Edge> es,vector<P> position,int root) : id(id), n(n), es(es), position(position), root(root){
		//convex_polygon = convex_hull(position);
		
		
		
		parent = vector<int>(n,-1);
		depth = vector<int>(n);
		child = vector<vector<int>>(n);
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
				}
			}
		}
		init_dfs(root);
	}
private:
	void init_dfs(int x){
		tot_sum_of_tree[x] = 0.0;
		vector<P> ps;
		ps.push_back(position[x]);
		for( auto c : child[x] ){
			init_dfs(c);
			tot_sum_of_tree[x] += tot_sum_of_tree[c];// + length_between_parent[c];
			ps.insert(ps.end(),mec[c].ps.begin(),mec[c].ps.end());
		}
		mec[x] = Circle::minEnclosingCircle(ps);
		tot_sum_of_tree[x] += length_between_parent[x];
		// cerr << tot_sum_of_tree[x] << endl;
		//cerr << ps.size() << endl;
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
				trees.push_back(Tree(trees.size(),points.size(),es,position,0));
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
			vs.push_back(line[0].real()+0.5);
			vs.push_back(line[0].imag()+0.5);
			vs.push_back(line[1].real()+0.5);
			vs.push_back(line[1].imag()+0.5);
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


	vector<Circle> mecs;
	ExtendedAnswer(Problem *problem) : problem(problem){
		double all_total = 0;
		for(const auto &tree : problem->trees ){
			all_total += tree.tot_sum_of_tree[0];
			// cerr << tree.tot_sum_of_tree[0] << "<" << endl;
			original_area.push_back(area2(tree.convex_polygon));
			convex_polygon.push_back(tree.convex_polygon);
			mecs.push_back(Circle::minEnclosingCircle(tree.convex_polygon));
			MEMO_score_of_tree_rough.push_back(tree.tot_sum_of_tree[0]);
			MEMO_score_of_tree.push_back(tree.tot_sum_of_tree[0]);
			current_weight.push_back(vector<double>(tree.n,-1.0));
			already_cut.push_back(vector<bool>(tree.n,false));
			
		}
		// cerr << all_total << "<" << endl;
		MEMO_overall_score_rough = MEMO_overall_score = all_total;
		
	}
	void add_line(const L &line){
		lines.push_back(line);
		// overall_score_rough(line,true);
		overall_score(line,true);
	}
	
	vector<double> MEMO_score_of_tree_rough;
	
	double score_of_tree_rough(const Tree &tree,const L &l,bool reflesh){
		if( distanceLP(l,mecs[tree.id].p) > mecs[tree.id].r + EPS ){
			return MEMO_score_of_tree_rough[tree.id];
		}
		double all_area = original_area[tree.id];
		if( all_area < EPS ) return MEMO_score_of_tree_rough[tree.id] = 0;
		G& cut_g = convex_polygon[tree.id];
		if( cut_g.size() == 0 ) return MEMO_score_of_tree_rough[tree.id] = 0;
		auto g1 = convex_cut(cut_g,l);
		G next_cut_g;
		if( convex_contains(g1,tree.position[tree.root]) != OUT ){
			next_cut_g = g1;
		}else{
			auto g2 = convex_cut(cut_g,L(l[1],l[0]));
			next_cut_g = g2;
		}
		double sub_area =  area2(next_cut_g);
		double res = (sub_area / all_area) * tree.tot_sum_of_tree[tree.root];
		
		if( reflesh ) {
			cut_g = next_cut_g;
			mecs[tree.id] = Circle::minEnclosingCircle(next_cut_g);
			MEMO_score_of_tree_rough[tree.id] = res;
		}
		return res;
	}
	
	double MEMO_overall_score_rough ;
	double overall_score_rough(const L &l = null_line,bool reflesh=false){
		if( l == null_line ){
			return MEMO_overall_score_rough;
		}
		double sum = 0;
		for(const auto &tree : problem->trees ){
			sum += score_of_tree_rough(tree,l,reflesh);
		}
		if( reflesh ){
			MEMO_overall_score_rough = sum;
		}
		return sum;
	}
	

	void draw_MEC(const Problem &problem, vector<RGB> colors){
		cout << "C " << 0 << " " << 0 << " " << 0 << " " << 512 << " " << 512 << " " << 512 << endl;
		for(int i = 0 ; i < convex_polygon.size() ; i++){
			auto &g = convex_polygon[i];
			if( g.size() >= 3 ){
				cout << "G ";
				cout << colors[i].r << " " << colors[i].g << " " << colors[i].b;
				for( auto p : g ){
					cout << " " << (int)(p.real()+0.5) << " " << (int)(p.imag()+0.5);
				}
				cout << endl;
			}
			cout << "C ";
			cout << colors[i].r << " " << colors[i].g << " " << colors[i].b;
			auto mec = Circle::minEnclosingCircle(g);
			// cerr << mec.p.real() << " " << mec.p.imag() << " " << mec.r << endl;
			cout << " " << (int)(mec.p.real()+0.5) << " " << (int)(mec.p.imag()+0.5) << " " << (int)(mec.r+0.5) << endl;
			cout << "C ";
			cout << 255-colors[i].r << " " << 255-colors[i].g << " " << 255-colors[i].b;
			cout << " " << (int)(problem.trees[i].position[0].real()+0.5) << " " << (int)(problem.trees[i].position[0].imag()+0.5) << " " << 2 << endl;
		}
		for( auto l : lines ){
			cout << "L " << 255 << " " << 255 << " " << 255 << " " << (int)(l[0].real()+0.5) << " " << (int)(l[0].imag()+0.5) << " " << (int)(l[1].real()+0.5) << " " << (int)(l[1].imag()+0.5) << endl;
		}
		cout << "END" << endl;
	
	}
	void draw(const Problem &problem, vector<RGB> colors){
		overall_score_rough();
		draw_MEC(problem,colors);
		return;
		cout << "C " << 0 << " " << 0 << " " << 0 << " " << 512 << " " << 512 << " " << 512 << endl;
		for(int i = 0 ; i < convex_polygon.size() ; i++){
			auto &g = convex_polygon[i];
			if( g.size() >= 3 ){
				cout << "G ";
				cout << colors[i].r << " " << colors[i].g << " " << colors[i].b;
				for( auto p : g ){
					cout << " " << (int)(p.real()+0.5) << " " << (int)(p.imag()+0.5);
				}
				cout << endl;
			}
			cout << "C ";
			cout << 255-colors[i].r << " " << 255-colors[i].g << " " << 255-colors[i].b;
			cout << " " << (int)(problem.trees[i].position[0].real()+0.5) << " " << (int)(problem.trees[i].position[0].imag()+0.5) << " " << 2 << endl;
		}
		for( auto l : lines ){
			cout << "L " << 255 << " " << 255 << " " << 255 << " " << (int)(l[0].real()+0.5) << " " << (int)(l[0].imag()+0.5) << " " << (int)(l[1].real()+0.5) << " " << (int)(l[1].imag()+0.5) << endl;
		}
		cout << "END" << endl;
	}


	vector<double> MEMO_score_of_tree;

	double score_of_tree(const Tree &tree,const L &l,bool reflesh=true){
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
	double inner_dfs_score_of_tree(int x,const Tree &tree,const L &line,bool reflesh){
		
		// if( distanceLP(line,mecs[tree.id].p) > mecs[tree.id].r + EPS ){
		if( distanceLP(line,tree.mec[x].p) > tree.mec[x].r + EPS ){
			if( current_weight[tree.id][x] != -1.0 ){
				return current_weight[tree.id][x];
			}
		}
		// cerr <<  already_cut[tree.id][x] << endl;
		if( already_cut[tree.id][x] ) return current_weight[tree.id][x];


		double sum = 0;
		for( auto c : tree.child[x] ){
			// ここ前計算で114514倍くらい速くなると思う
			if( intersectLS(line,L(tree.position[x],tree.position[c])) ){
				P cp = crosspoint(line,L(tree.position[x],tree.position[c]));
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
	
	static double score_of_tree_fast(const Tree &tree,const Answer &answer){
		double all_area = area2(tree.convex_polygon);
		if( all_area < EPS ) return 0;
		G cut_g = tree.convex_polygon;
		for(auto l : answer.lines ){
			if( cut_g.size() == 0 ) return 0;
			auto g1 = convex_cut(cut_g,l);
			if( convex_contains(g1,tree.position[tree.root]) != OUT ){
				cut_g = g1;
			}else{
				auto g2 = convex_cut(cut_g,L(l[1],l[0]));
				cut_g = g2;
			}	
		}
		double sub_area =  area2(cut_g);
		return (sub_area / all_area) * tree.tot_sum_of_tree[tree.root];
	}
	static double overall_score_fast(const Problem &problem,const Answer &answer){
		double sum = 0;
		for( const auto &tree : problem.trees ){
			sum += score_of_tree_fast(tree,answer);
		}
		return sum;
	}

	static double score_of_tree_fast_differ_ver(int tree_id,const Tree &tree,ExtendedAnswer &answer,const L &l){
		double all_area = answer.original_area[tree_id];
		if( all_area < EPS ) return 0;
		G& cut_g = answer.convex_polygon[tree_id];
		if( cut_g.size() == 0 ) return 0;
		auto g1 = convex_cut(cut_g,l);
		if( convex_contains(g1,tree.position[tree.root]) != OUT ){
			cut_g = g1;
		}else{
			auto g2 = convex_cut(cut_g,L(l[1],l[0]));
			cut_g = g2;
		}
		double sub_area =  area2(cut_g);
		return (sub_area / all_area) * tree.tot_sum_of_tree[tree.root];
	}
	static double overall_score_fast_differ_ver(const Problem &problem,ExtendedAnswer &answer,const L &l){
		double sum = 0;
		for(int i = 0 ; i < problem.trees.size() ; i++){
			sum += score_of_tree_fast_differ_ver(i,problem.trees[i],answer,l);
		}
		answer.lines.push_back(l);
		return sum;
	}
	static double overall_score_fast_nonline(const Problem &problem,ExtendedAnswer &answer){
		double sum = 0;
		for(int i = 0 ; i < problem.trees.size() ; i++){
			if( answer.original_area[i] > EPS ){
				sum += area2(answer.convex_polygon[i]) / answer.original_area[i] * problem.trees[i].tot_sum_of_tree[problem.trees[i].root];
			}
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
		return L(P(-1,-1),P(-1,-1));
	}
};
