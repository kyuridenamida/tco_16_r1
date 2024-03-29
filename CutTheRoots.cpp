﻿#include "classes_for_tco.cpp"
#include <cstdio>

const double TIME_LIMIT = 9.85;

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
		// srand(time(NULL));
		resetTimer();
		Problem problem(NP,points,roots);
		// if( configuration.visualize_mode ){
			// ExtendedAnswer answer_tmp(&problem);
			// answer_tmp.draw_sol();
		// }
		
		// if( configuration.visualize_mode ){
			// ExtendedAnswer answer_tmp(&problem);
			// answer_tmp.draw_sol();
			// return {};
		// }
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
	Problem *preprocess(Problem &problem){
		int N = problem.trees.size();

		vector<pair<double,pair<int,int>>> pairs;
		for(int i = 0 ; i < N ; i++){
			for(int j = i+1; j < N ; j++){
				P p1 = problem.trees[i].position[problem.trees[i].root];
				P p2 = problem.trees[j].position[problem.trees[j].root];
				pairs.push_back({norm(p1-p2),{i,j}});
			}
		}
		
		sort(pairs.begin(),pairs.end());
		
		
		vector<ExtendedAnswer> answers;
		
		for(int _ = 0 ; _ < 5 ; _++){
			ExtendedAnswer answer(&problem);
			
			for(int li = 0 ; li < pairs.size(); li++){
				int i = pairs[li].second.first;
				int j = pairs[li].second.second;
				P p1 = problem.trees[i].position[problem.trees[i].root];
				P p2 = problem.trees[j].position[problem.trees[j].root];
				bool separated = false;
				for( auto &l : answer.lines ){
					if( GeomUtils::is_separating(l,p1,p2) ){
						separated = true;
						break;
					}
				}
				
				if( !separated ){
					vector<pair<double,L>> cand;
					double cur = answer.overall_score();
					
					for(int m = 0 ; m < 8000 / problem.trees.size() ; m++){
						double A = 1. * xor64() / (unsigned int)(-1);
						double B = 1. * xor64() / (unsigned int)(-1);
						P mp = p1 + (p2-p1) * A;
					// cerr << cur << " " << NaiveScoring::overall_score_fast(problem,answer) << "|" <<  NaiveScoring::overall_score(problem,answer) << endl;
						complex<double> rot = exp(complex<double>(0,PI*B));
						P vec = (p1-p2) * P(rot.real(),rot.imag());
						L l = L(mp,mp+vec);
						if( l.a != P(-1,-1) and GeomUtils::is_separating(l,p1,p2) ){
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
						if( fix_l.a != null_point and GeomUtils::is_separating(fix_l,p1,p2) ){
								// double loss = answer.overall_score() - answer.overall_score(fix_l);
								// cerr << loss -  cand[0].first << endl;
								answer.add_line(fix_l);
								
								break;
						}
					}
				}
			}
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
						for(int li = 0 ; li < ls.size() ; li++){
							const L &l = ls[li];
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
				}else i++;
			}
			ExtendedAnswer hogehoge(&problem);
			for(int i = 0 ; i < answer.lines.size() ; i++){
				hogehoge.add_line(answer.lines[i]);
			}
			answers.push_back(hogehoge);
			random_shuffle(pairs.begin(),pairs.end());
		}
		vector<int> idx;
		for(int i = 0 ; i < answers.size() ; i++)
			idx.push_back(i);
		sort(idx.begin(),idx.end(),[&](int a,int b){
			
			return answers[a].overall_score() > answers[b].overall_score();
		});
		// if( idx.size() > 3 ) idx.resize(3);
		vector<ExtendedAnswer> res_ans;
		for(int i = 0 ; i < idx.size() ; i++)
			res_ans.push_back(answers[idx[i]]);
		
		return remake_trees(&problem,res_ans);
	}
	
	ExtendedAnswer greedy1(Problem &problem){
		Problem problem_processed = *preprocess(problem);
		int N = problem.trees.size();

		vector<pair<double,pair<int,int>>> pairs;
		for(int i = 0 ; i < N ; i++){
			for(int j = i+1; j < N ; j++){
				P p1 = problem.trees[i].position[problem.trees[i].root];
				P p2 = problem.trees[j].position[problem.trees[j].root];
				pairs.push_back({norm(p1-p2),{i,j}});
			}
		}
		sort(pairs.begin(),pairs.end());
		
		ExtendedAnswer answer(&problem_processed);
		for(int li = 0 ; li < pairs.size(); li++){
			int i = pairs[li].second.first;
			int j = pairs[li].second.second;
			P p1 = problem_processed.trees[i].position[problem_processed.trees[i].root];
			P p2 = problem_processed.trees[j].position[problem_processed.trees[j].root];
			bool separated = false;
			for(int li = 0 ; li < answer.lines.size() ; li++){
				const L &l = answer.lines[li];
				if( GeomUtils::is_separating(l,p1,p2) ){
					separated = true;
					break;
				}
			}
			
			if( !separated ){
				vector<pair<double,L>> cand;
				double cur = answer.overall_score();
				for(int m = 0 ; m < 80000 / problem_processed.trees.size() ; m++){
					double A = 1. * xor64() / (unsigned int)(-1);
					double B = 1. * xor64() / (unsigned int)(-1);
					P mp = p1 + (p2-p1) * A;
					complex<double> rot = exp(complex<double>(0,PI*B));
					P vec = (p1-p2) * P(rot.real(),rot.imag());
					L l = L(mp,mp+vec);
					if( l.a != P(-1,-1) and GeomUtils::is_separating(l,p1,p2) ){
						double loss = cur - answer.overall_score(l);
						cand.emplace_back(loss,l);
					}
				}
				sort(cand.begin(),cand.end(),[&](const pair<double,L> &a,const pair<double,L> &b){
					return a.first < b.first;
				});
				for(int i = 0 ; i < cand.size() ; i++){
					L fix_l = GeomUtils::convert_to_integer_line(cand[i].second,p1,p2);
					if( fix_l.a != null_point and GeomUtils::is_separating(fix_l,p1,p2) ){
							answer.add_line(fix_l);
							// answer.draw_sol();
							break;
					}
				}
			}
		}
		
		// return answer;
		answer = answer.refined_answer();
		ExtendedAnswer tmp(&problem);
		for( auto &l : answer.lines )
			tmp.add_line(l);
		answer = tmp;

		
		while( nowTime() < TIME_LIMIT  ){
			double start = nowTime();
			// if( problem.NR < 60000  ){
				// random_shuffle(pairs.begin(),pairs.end());
			// }
			random_shuffle(answer.lines.begin(),answer.lines.end());
			
			for(int xxx = 0 ; xxx <  answer.lines.size() ; xxx++){
				if( nowTime() > TIME_LIMIT  ) return answer;
				
				ExtendedAnswer train(&problem);
				for(int i = 0 ; i < answer.lines.size() ; i++){
					if( xxx != i  ) train.add_line(answer.lines[i]);
				}
				for(int li = 0 ; li < pairs.size(); li++){
					int i = pairs[li].second.first;
					int j = pairs[li].second.second;
					P p1 = problem.trees[i].position[problem.trees[i].root];
					P p2 = problem.trees[j].position[problem.trees[j].root];
					bool separated = false;
					for( auto &l : train.lines ){
						if( GeomUtils::is_separating(l,p1,p2) ){
							separated = true;
							break;
						}
					}
					
					if( !separated ){
						vector<pair<double,L>> cand;
						double cur = train.overall_score();
						
						for(int m = 0 ; m < 20000 / problem.trees.size() ; m++){
							double A = 1. * xor64() / (unsigned int)(-1);
							
							double B = 1. * xor64() / (unsigned int)(-1);
							P mp = p1 + (p2-p1) * A;
							complex<double> rot = exp(complex<double>(0,PI*B));
							P vec = (p1-p2) * P(rot.real(),rot.imag());
							L l = L(mp,mp+vec);
							if( l.a != P(-1,-1) and GeomUtils::is_separating(l,p1,p2) ){
								double loss = cur - train.overall_score(l);
								cand.emplace_back(loss,l);
							}
						}
						sort(cand.begin(),cand.end(),[&](const pair<double,L> &a,const pair<double,L> &b){
							return a.first < b.first;
						});
						for(int i = 0 ; i < cand.size() ; i++){
							L fix_l = GeomUtils::convert_to_integer_line(cand[i].second,p1,p2);
							if( fix_l.a != null_point and GeomUtils::is_separating(fix_l,p1,p2) ){
									train.add_line(fix_l);
									break;
							}
						}
					}
				}
				train = train.refined_answer();
				if( train.overall_score() > answer.overall_score()  ){
					//cerr << train.overall_score() - answer.overall_score() << "|" << train.overall_score() << "|" << xxx << endl;
					answer = train;
					// if( configuration.visualize_mode ) answer.draw_sol();
					break;
				}
						
			}
			//cerr << nowTime() - start << endl;
			
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
