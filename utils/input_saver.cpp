#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stack>
#include <cassert>
#include <sstream>
#include <vector>
#include <fstream>
using namespace std;

#ifdef LOCAL


int main(int argc,char *argv[]) {
	if( argc != 2 ){
		cerr << "Something is wrong!" << endl;
		return 1;
	}
	ofstream ofs(argv[1]);
	auto getVector = [&](vector<int>& v) {
		for (int i = 0; i < v.size(); ++i){
			cin >> v[i];
			ofs << v[i] << (i+1 == v.size()?"\n":" ");
		}
	};

    int NP;
    cin >> NP;
	ofs << NP << endl;
    int Npoints;
    cin >> Npoints;
	ofs << Npoints << endl;
    vector<int> points(Npoints);
    getVector(points);

    int Nroots;
    cin >> Nroots;
	ofs << Nroots << endl;
    vector<int> roots(Nroots);
    getVector(roots);
}
#endif
