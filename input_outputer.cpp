#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <map>
#include <stack>
#include <cassert>
#include <sstream>
#include <vector>
using namespace std;

#ifdef LOCAL
template<class T> void getVector(vector<T>& v) {
    for (int i = 0; i < v.size(); ++i){
        cin >> v[i];
		cerr << v[i] << " ";
	}
	cerr << endl;
}


int main() {
    int NP;
    cin >> NP;
	cerr << NP << endl;
    int Npoints;
    cin >> Npoints;
	cerr << Npoints << endl;
    vector<int> points(Npoints);
    getVector(points);

    int Nroots;
    cin >> Nroots;
	cerr << Nroots << endl;
    vector<int> roots(Nroots);
    getVector(roots);
}
#endif
