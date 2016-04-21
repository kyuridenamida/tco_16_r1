/* 超注意：　ほとんどのアルゴリズムは多角形が反時計回りであることを仮定する． */
#include <vector>
#include <complex>
#include <algorithm>
#include <cassert>
#include <list>
#include <array>
#include <cmath>
using namespace std;
static const double EPS = 1e-9;
const double INF = 1e12;
const double PI = acos(-1);
#define REP(i, n) for ( int i = 0; i < (n); i++ )

struct P{
	double x,y;
	P(){x=y=0;}
	P(double x,double y) : x(x), y(y) {}
	P operator+(const P &t) const {
		return P(x + t.x, y + t.y);
	}
	P operator-(const P &t) const {
		return P(x - t.x, y - t.y);
	}
	P operator*(double t) const {
		return P(x * t, y * t);
	}
	P operator*(const P &t) const {
		return P(x*t.x-y*t.y,y*t.x+x*t.y);
	}
	P operator/(double t) const {
		return P(x / t, y / t);
	}
	bool operator == (const P &t) const{
		return t.x == x and t.y == y;
	}
	bool operator != (const P &t) const{
		return t.x != x or t.y != y;
	}
};
P operator * (double t,const P &a){
	return P(a.x * t, a.y * t);
}
typedef vector<P> G;

bool operator < (const P& a, const P& b) {
	return a.x != b.x ? a.x < b.x : a.y < b.y;
}

inline double cross(const P &a,const P &b) { return a.x * b.y - a.y * b.x; }
inline double dot(const P &a,P const &b) { return a.x * b.x + a.y * b.y; }
inline double norm(const P &a){ return a.x * a.x + a.y * a.y; }
inline double abs(const P &a){ return sqrt(a.x * a.x + a.y * a.y); }


//直線
struct L{
	P a,b;
	double A,B,C;
	L(const P &a,const P &b) : a(a), b(b) {
		double X0 = a.x;
		double Y0 = a.y;
		double X1 = b.x;
		double Y1 = b.y;
		A = (Y1-Y0);
		B = (X0-X1);
		C = Y0 * (X1-X0) + X0 * (Y0-Y1);
	}
	bool operator == (const L &t) const{
		return t.a == a and t.b == b;
	}
	bool operator != (const L &t) const{
		return t.a != a or t.b != b;
	}

	bool operator < (const L& t) const {
		return a != t.a ? a < t.a : b < t.b;
	}
};


int ccw(P a, P b, P c) {
	b = b - a; c = c - a;
	if (cross(b, c) > 0)   return +1;	// a → b で反時計方向に折れて b → c(？)
	if (cross(b, c) < 0)   return -1;	// a → b で時計方向に折れて b → c(？)
	if (dot(b, c) < 0) return +2;    	// a→bで逆向いてaを通り越してb→c(c--a--b)
	if (norm(b) < norm(c)) return -2;	// a→bでそのままb→c(a--b--c)
	return 0;							// a→bで逆向いてb→c(または b == c)
}

bool intersectLS(const L &l, const L &s) {
	return (l.A*s.a.x+l.B*s.a.y+l.C) * (l.A*s.b.x+l.B*s.b.y+l.C) < EPS;
	//return cross(l.b-l.a, s.a-l.a)*cross(l.b-l.a, s.b-l.a) < EPS;
}



inline double distanceLP_check(const L &l, const P &p,const double &r2) {
	return (l.A*p.x+l.B*p.y+l.C)*(l.A*p.x+l.B*p.y+l.C) > r2 * (l.A*l.A+l.B*l.B) + EPS;
}

inline double distanceLP(const L &l, const P &p) {
	return (l.A*p.x+l.B*p.y+l.C)*(l.A*p.x+l.B*p.y+l.C) / (l.A*l.A+l.B*l.B);
}


 
P crosspoint(const L &l, const L &m) {
	// double y = - (l.C * m.A - l.A * m.C) / (l.B * m.A - l.A * m.B);
	// double x = (-l.B * y - l.C) / l.A;
	// return P(x,y);
	
	double A = cross(l.b - l.a, m.b - m.a);
	double B = cross(l.b - l.a, l.b - m.a);
	if (abs(A) < EPS && abs(B) < EPS) return m.a;
	return m.a + B / A * (m.b - m.a);
}
 
#define curr(P, i) P[(i) % P.size()]
#define next(P, i) P[(i+1)%P.size()]
#define prev(P, i) P[(i+P.size()-1) % P.size()]

// 単純多角形の面積の"2倍"を求める O(n)
double area2(const G& poly) {
	double A = 0;
	REP(i,poly.size()) A += cross(curr(poly, i), next(poly, i));
	return A;
}
// 凸包を求める O(n log n)
vector<P> convex_hull2(vector<P> ps) {
  if( ps.size() < 3 ) return ps;
  int n = ps.size(), k = 0;
  sort(ps.begin(), ps.end());

  vector<P> ch(2*n);
  for (int i = 0; i < n; ch[k++] = ps[i++]) // lower-hull
    while (k >= 2 && ccw(ch[k-2], ch[k-1], ps[i]) <= 0) --k;
  for (int i = n-2, t = k+1; i >= 0; ch[k++] = ps[i--]) // upper-hull
    while (k >= t && ccw(ch[k-2], ch[k-1], ps[i]) <= 0) --k;
  ch.resize(k-1);
  return ch;
}


//http://d.hatena.ne.jp/TobiasGSmollett/20150220/1424445987
//お借りした
struct Circle{
  P p;
  double r2;
  vector<P> ps;
  Circle(){}
  Circle(P p, double r2) : p(p) , r2(r2){}
  Circle(P p, double r2,vector<P> ps) : p(p) , r2(r2), ps(ps) {}
  
  
  bool contain(P a){
    return norm(a-p) <= r2;
  }
  
  static Circle circumCircle(P a,P b){
    P q=(a+b)/2.0;
	return Circle(q,norm(a-q));
  }
  
  static Circle circumscribedCircle(P p, P q, P r){
    P a=(q-p)*2.0,b=(r-p)*2.0;
    P c(dot(p,p)-dot(q,q),dot(p,p)-dot(r,r));
    Circle res;
	double x = a.y*c.y-b.y*c.x;
    double y = b.x*c.x-a.x*c.y;
    res.p = P(x,y)/cross(a,b);
    return Circle(res.p, norm(p-res.p));
  }

  static Circle minEnclosingCircle(vector<P>ps){
    if(ps.size()==0) return Circle(P(0,0),0,{});
    if(ps.size()==1) return Circle(ps[0],0,{ps[0]});
	if(ps.size()==2){
		Circle c = circumCircle(ps[0],ps[1]);
		c.ps = {ps[0],ps[1]};
		return c;
	}
    Circle circle=circumscribedCircle(ps[0],ps[1],ps[2]);
	circle.ps = {ps[0],ps[1],ps[2]};
	
    for(int i=2;i<ps.size();i++){
      if(!circle.contain(ps[i])){
			circle=circumscribedCircle(ps[0],ps[1],ps[i]);
			circle.ps = {ps[0],ps[1],ps[i]};
			for(int j=1;j<i;j++){
			  if(!circle.contain(ps[j])){
				circle=circumscribedCircle(ps[0],ps[j],ps[i]);
				circle.ps = {ps[0],ps[j],ps[i]};
				for(int k=0;k<j;k++){
				  if(!circle.contain(ps[k])){
					circle=circumscribedCircle(ps[i],ps[j],ps[k]);
					circle.ps = {ps[i],ps[j],ps[k]};
				  }
				}
			  }
			}
      }
    }
    return circle;
  }

};