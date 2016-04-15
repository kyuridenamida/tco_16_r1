#include <bits/stdc++.h>
using namespace std;



int main(){
	set<pair<int,int>> S;
	for(int i = 0 ; i <= 1000 ; i++){
		for(int j = 0 ; j <= 1000 ; j++){
			if( max(i,j) ){
				int I = i / __gcd(i,j);
				int J = j / __gcd(i,j);
				S.insert({I,J});
			}
		}
	}
	cout << S.size() << endl;
}