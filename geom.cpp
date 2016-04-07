/* 超注意：　ほとんどのアルゴリズムは多角形が反時計回りであることを仮定する． */
#include <vector>
#include <complex>
#include <algorithm>
#include <cmath>
using namespace std;
static const double EPS = 1e-10;
const double INF = 1e12;
const double PI = acos(-1);
#define REP(i, n) for ( int i = 0; i < (n); i++ ) 

typedef complex<double> P;
typedef vector<P> G;

namespace std {
	bool operator < (const P& a, const P& b) {
		return real(a) != real(b) ? real(a) < real(b) : imag(a) < imag(b);
	} 
}
double cross(P a,P b) { return imag(conj(a)*b); }
double dot(P a,P b) { return real(conj(a)*b);}

//直線
struct L : public vector<P> {
	L(P a,P b) {
		push_back(a); push_back(b);
	}
	
};

//円
struct C {
	P p;
	double r;
	C(const P &p, double r) : p(p), r(r) { }
	C(){} 
};

int ccw(P a, P b, P c) {
	b -= a; c -= a;
	if (cross(b, c) > 0)   return +1;	// a → b で反時計方向に折れて b → c(？)
	if (cross(b, c) < 0)   return -1;	// a → b で時計方向に折れて b → c(？)
	if (dot(b, c) < 0) return +2;    	// a→bで逆向いてaを通り越してb→c(c--a--b)
	if (norm(b) < norm(c)) return -2;	// a→bでそのままb→c(a--b--c)
	return 0;							// a→bで逆向いてb→c(または b == c)
}
 
bool intersectLL(const L &l, const L &m) {
	return abs(cross(l[1]-l[0], m[1]-m[0])) > EPS || 
	abs(cross(l[1]-l[0], m[0]-l[0])) < EPS;
}
bool intersectLS(const L &l, const L &s) {
	return cross(l[1]-l[0], s[0]-l[0])*cross(l[1]-l[0], s[1]-l[0]) < EPS;
}
bool intersectLP(const L &l, const P &p) {
	return abs(cross(l[1]-p, l[0]-p)) < EPS;
}
bool intersectSS(const L &s, const L &t) {	// <= 0を < にすると通ることがあるらしい。
	return	ccw(s[0],s[1],t[0])*ccw(s[0],s[1],t[1]) <= 0 && 
			ccw(t[0],t[1],s[0])*ccw(t[0],t[1],s[1]) <= 0;
}
bool intersectSP(const L &s, const P &p) {
	return abs(s[0]-p)+abs(s[1]-p)-abs(s[1]-s[0]) < EPS;
}
// 垂線の足
P projection(const L &l, const P &p) {
	double t = dot(p-l[0], l[0]-l[1]) / norm(l[0]-l[1]);
	return l[0] + t*(l[0]-l[1]);
}
// 直線に対称な点
P reflection(const L &l, const P &p) {
	return p + 2.0 * (projection(l, p) - p);
}

double distanceLP(const L &l, const P &p) {
	return abs(p - projection(l, p));
}
double distanceLL(const L &l, const L &m) {
	return intersectLL(l, m) ? 0 : distanceLP(l, m[0]);
}
double distanceLS(const L &l, const L &s) {
	if (intersectLS(l, s)) return 0;
	return min(distanceLP(l, s[0]), distanceLP(l, s[1]));
}
double distanceSP(const L &s, const P &p) {
	const P r = projection(s, p);
	if (intersectSP(s, r)) return abs(r - p);
	return min(abs(s[0] - p), abs(s[1] - p));
}
double distanceSS(const L &s, const L &t) {
	if (intersectSS(s, t)) return 0;
	return min(min(distanceSP(s, t[0]), distanceSP(s, t[1])),min(distanceSP(t, s[0]), distanceSP(t, s[1])));
}
 
P crosspoint(const L &l, const L &m) {
	double A = cross(l[1] - l[0], m[1] - m[0]);
	double B = cross(l[1] - l[0], l[1] - m[0]);
	if (abs(A) < EPS && abs(B) < EPS) return m[0];
	return m[0] + B / A * (m[1] - m[0]);
}
 
#define curr(P, i) P[(i) % P.size()]
#define next(P, i) P[(i+1)%P.size()]
#define prev(P, i) P[(i+P.size()-1) % P.size()]
enum { OUT, ON, IN };
// 点-多角形包含判定(凸でなくとも良い) O(n) ロバストなことを確認
int contains(const G& polygon, const P& p) {
	bool in = false;
	for (int i = 0; i < polygon.size(); ++i) {
		P a = curr(polygon,i) - p, b = next(polygon,i) - p;
		if (imag(a) > imag(b)) swap(a, b);
		if (imag(a) <= 0 && 0 < imag(b))
			if (cross(a, b) < 0) in = !in;
		if (cross(a, b) == 0 && dot(a, b) <= 0) return ON;
	}
	return in ? IN : OUT;
}

// 点-多角形包含判定(凸じゃないといけない) O(log n)
int convex_contains(const G &polygon, const P &p) {
  const int n = polygon.size();
  P g = (polygon[0] + polygon[n/3] + polygon[2*n/3]) / 3.0; // inner-point
  int a = 0, b = n;
  while (a+1 < b) { // invariant: c is in fan g-P[a]-P[b]
    int c = (a + b) / 2;
    if (cross(polygon[a]-g, polygon[c]-g) > 0) { // angle < 180 deg
      if (cross(polygon[a]-g, p-g) > 0 && cross(polygon[c]-g, p-g) < 0) b = c;
      else a = c;
    } else {
      if (cross(polygon[a]-g, p-g) < 0 && cross(polygon[c]-g, p-g) > 0) a = c;
      else b = c;
    }
  }
  b %= n;
  if (cross(polygon[a] - p, polygon[b] - p) < 0) return 0;
  if (cross(polygon[a] - p, polygon[b] - p) > 0) return 2;
  return 1;
}

//凸多角形の共通部分(2関数) 怪しいらしい. O(n + m)．
bool intersect_1pt(const P& a, const P& b,
                   const P& c, const P& d, P &r) {
  double D =  cross(b - a, d - c);
  if (fabs(D) < EPS) return false;
  double t =  cross(c - a, d - c) / D;
  double s = -cross(a - c, b - a) / D;
  r = a + t * (b - a);
  return t >= -EPS && t <= 1 + EPS && s >= -EPS && s <= 1 + EPS;
}
G convex_intersect(const G &X, const G &Q) {
  const int n = X.size(), m = Q.size();
  int a = 0, b = 0, aa = 0, ba = 0;
  enum { Xin, Qin, Unknown } in = Unknown;
  G R;
  do {
    int a1 = (a+n-1) % n, b1 = (b+m-1) % m;
    double C = cross(X[a] - X[a1], Q[b] - Q[b1]);
    double A = cross(X[a1] - Q[b], X[a] - Q[b]);
    double B = cross(Q[b1] - X[a], Q[b] - X[a]);
    P r;
    if (intersect_1pt(X[a1], X[a], Q[b1], Q[b], r)) {
      if (in == Unknown) aa = ba = 0;
      R.push_back( r );
      in = B > 0 ? Xin : A > 0 ? Qin : in;
    }
    if (C == 0 && B == 0 && A == 0) {
      if (in == Xin) { b = (b + 1) % m; ++ba; }
      else           { a = (a + 1) % m; ++aa; }
    } else if (C >= 0) {
      if (A > 0) { if (in == Xin) R.push_back(X[a]); a = (a+1)%n; ++aa; }
      else       { if (in == Qin) R.push_back(Q[b]); b = (b+1)%m; ++ba; }
    } else {
      if (B > 0) { if (in == Qin) R.push_back(Q[b]); b = (b+1)%m; ++ba; }
      else       { if (in == Xin) R.push_back(X[a]); a = (a+1)%n; ++aa; }
    }
  } while ( (aa < n || ba < m) && aa < 2*n && ba < 2*m );
  if (in == Unknown) {
    if (convex_contains(Q, X[0])) return X;
    if (convex_contains(X, Q[0])) return Q;
  }
  return R;
}
// 単純多角形の面積の"2倍"を求める O(n)
double area2(const G& poly) {
	double A = 0;
	REP(i,poly.size()) A += cross(curr(poly, i), next(poly, i));
	return A;
}
// 凸包を求める O(n log n)
vector<P> convex_hull(vector<P> ps) {
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
//円同士の交点
vector<P> C_cp(C a,C b){
	vector<P> ret;
	double L = abs(a.p-b.p);

	if(	L-a.r-b.r > EPS || (abs(a.p-b.p)<EPS && fabs(a.r-b.r)<EPS) || 
		abs(a.p-b.p) < abs(a.r-b.r)
	)return ret;
 
	double theta = atan2(b.p.imag()-a.p.imag(),b.p.real()-a.p.real());
	double c = (L*L+a.r*a.r-b.r*b.r)/(2*L*a.r);
	ret.push_back(
		P(a.p.real()+a.r*cos(theta+acos(c)),
		  a.p.imag()+a.r*sin(theta+acos(c)))
	);
	if(fabs(L-a.r-b.r) > EPS)
		ret.push_back(
			P(a.p.real()+a.r*cos(theta-acos(c)),
			  a.p.imag()+a.r*sin(theta-acos(c)))
		);
	return ret;
}


// 垂線の足元を求める。
P getPedal(L l, P p){
	double A;
	if(abs(l[1].real()-l[0].real()) < EPS){
		return P(l[1].real(),p.imag()); // important
	}else{
		A = (l[1].imag()-l[0].imag())/(l[1].real()-l[0].real());
	}
	double a = -A , b = 1 , c = A*l[0].real() - l[0].imag();
	double t = (a*p.real() + b*p.imag() + c)/(a*a+b*b);
	return p-t * P(a,b);
}
 
// 円と直線の交点
vector<P> crosspointCL(const L l, const C c){
	vector<P> ret;
	P p = getPedal(l,c.p);
	if(	abs(p-c.p) > c.r+EPS)return ret;
	P e = P((l[1]-l[0])/abs(l[1]-l[0]));
	double S = sqrt(c.r*c.r-abs(p-c.p)*abs(p-c.p));
	ret.push_back(p+S*e);
	ret.push_back(p-S*e);
	return ret;
}
 
// 外接円の中心を求める。
P getCircumcenter(P a,P b,P c){
	double A1 = 2 * ( b.real() - a.real() );
	double B1 = 2 * ( b.imag() - a.imag() );
	double C1 = pow(a.real(),2)-pow(b.real(),2) + pow(a.imag(),2)-pow(b.imag(),2);
	double A2 = 2 * ( c.real() - a.real() );
	double B2 = 2 * ( c.imag() - a.imag() );
	double C2 = pow(a.real(),2)-pow(c.real(),2) + pow(a.imag(),2)-pow(c.imag(),2);
	double X = (B1 * C2 - B2 * C1) / (A1 * B2 - A2 * B1);
	double Y = (C1 * A2 - C2 * A1) / (A1 * B2 - A2 * B1);
	return P(X,Y);
}
//内接円に関する知識
//・半径=2*面積/周長
//・中心座標=3つの座標(xa,ya),(xb,yb),(xc,yc)について3辺の長さをa,b,cとおくと座標の重みつき平均値の和で求まる。
//a/(a+b+c)*(xa,ya)+...

// 凸多角形の切断。凸多角形をある直線で切断し，
//その左側だけ残す．直線よりも左側にあるものと，交差するときはその交点を出力として吐いている．
G convex_cut(const G& poly, const L& l) {
  G Q;
  for (int i = 0; i < poly.size(); ++i) {
    P A = curr(poly, i), B = next(poly, i);
    if (ccw(l[0], l[1], A) != -1) Q.push_back(A);
    if (ccw(l[0], l[1], A)*ccw(l[0], l[1], B) < 0)
      Q.push_back(crosspoint(L(A, B), l));
  }
  return Q;
}
// 最近点対
pair<P,P> closestPair(vector<P> p) {
  int n = p.size(), s = 0, t = 1, m = 2, S[n]; S[0] = 0, S[1] = 1;
  sort(p.begin(),p.end()); // "p < q" <=> "p.x < q.x"
  double d = norm(p[s]-p[t]);
  for (int i = 2; i < n; S[m++] = i++) for(int j = 0 ; j < m ; j++) {
    if (norm(p[S[j]]-p[i])<d) d = norm(p[s = S[j]]-p[t = i]);
    if (real(p[S[j]]) < real(p[i]) - d) S[j--] = S[--m];
  }
  return make_pair( p[s], p[t] );
}

//ボロノイ図 正しく動いた on dxlib
//全体の領域をRとする。
//ボロノイ領域を求めたい点と、他の点に対して、垂直二等分線でRを切断し、注目点を含む方の領域を新たにRとする。
//これを全ての点に対して繰り返し行えば、ボロノイ領域を得ることができる。
//つまり、gは領域の多角形で、
//voronoi_cellをs=全てのiに対して試せば多角形郡が構築される

L bisector(P a, P b) {
  P A = (a+b)*P(0.5,0);
  return L(A, A+(b-a)*P(0, PI/2));
}
G voronoi_cell(G g, vector<P> v, int s) {
  for(int i = 0 ; i < v.size() ; i++){
    if (i!=s)
      g = convex_cut(g, bisector(v[s], v[i]));
  }
  return g;
}

//共通接線　正しく動いた。on dxlib
//半径0とすれば接線を求めるのに使えることを確認。
bool eq(double a,double b){
	return fabs(a-b) < 1e-9;
}
vector<L> contact(C p, C q){
  vector<L> ret;
  if(p.r < q.r) swap(p, q);
  double d = abs(p.p - q.p);
  P n = q.p - p.p;
  n /= abs(n);

  if(d + EPS < abs(p.r - q.r)){
    ret.clear();
  } else if(eq(d, abs(p.r - q.r))){
    P t, u;
    t = p.p + p.r * n;
    u = t + n * P(0, 1);
    ret.push_back(L(t, u));
  } else {
    if(!eq(p.r, q.r)){
      P t = p.p + (p.r * d / (p.r - q.r)) * n;
      long double theta = asin((p.r - q.r) / d);
      P u = n * P(cos(theta), sin(theta));
      P v = n * P(cos(-theta), sin(-theta));
      u += t;
      v += t;
      ret.push_back(L(t, u));
      ret.push_back(L(t, v));
    } else {
      P t = p.p + n * P(0, 1) * p.r;
      P u = p.p - n * P(0, 1) * p.r;
      ret.push_back(L(t, t+n));
      ret.push_back(L(u, u+n));
    }

    if(eq(d, p.r + q.r)){
      P t, u;
      t = p.p + p.r * n;
      u = t + n * P(0, 1);
      ret.push_back(L(t, u));
    } else if(d > p.r + q.r){
      P t = p.p + (p.r * d / (p.r + q.r)) * n;
      long double theta = asin((p.r + q.r) / d);
      P u = n * P(cos(theta), sin(theta));
      P v = n * P(cos(-theta), sin(-theta));
      u += t;
      v += t;
      ret.push_back(L(t, u));
      ret.push_back(L(t, v));
    }
  }
  return ret;
}
//単純多角形の三角形分割 (耳分解) O(n)
typedef vector<P> triangle;
triangle make_triangle(const P& a, const P& b, const P& c) {
  triangle ret(3);
  ret[0] = a; ret[1] = b; ret[2] = c;
  return ret;
}
bool triangle_contains(const triangle& tri, const P& p) {
  return ccw(tri[0], tri[1], p) >= 0 &&
         ccw(tri[1], tri[2], p) >= 0 &&
         ccw(tri[2], tri[0], p) >= 0;
}
bool ear_Q(int i, int j, int k, const G& P) {
  triangle tri = make_triangle(P[i], P[j], P[k]);
  if (ccw(tri[0], tri[1], tri[2]) <= 0) return false;
  for (int m = 0; m < P.size(); ++m)
    if (m != i && m != j && m != k)
      if (triangle_contains(tri, P[m]))
        return false;
  return true;
}
void triangulate(const G& P, vector<triangle>& t) {
  const int n = P.size();
  vector<int> l, r;
  for (int i = 0; i < n; ++i) {
    l.push_back( (i-1+n) % n );
    r.push_back( (i+1+n) % n );
  }
  int i = n-1;
  while (t.size() < n-2) {
    i = r[i];
    if (ear_Q(l[i], i, r[i], P)) {
      t.push_back(make_triangle(P[l[i]], P[i], P[r[i]]));
      l[ r[i] ] = l[i];
      r[ l[i] ] = r[i];
    }
  }
}