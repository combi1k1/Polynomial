#pragma GCC optimize "-O3"
#include<bits/stdc++.h>

using namespace std;

#define pb      push_back
#define sz(x)   x.size()
#define all(x)  x.begin(),x.end()

const int   N   = (1 << 18);
const int   mod = 998244353;

typedef long long   ll;
typedef vector<int> vi;

int read()  {
	int  x = 0;
	char c = getchar();
	while(c < '0' || c > '9')   c = getchar();
	while(48 <= c && c <= 57)   {
		x = x * 10 + c - '0';
		c = getchar();
	}
	return x;
}

int add(int x,int y)    {
	x += y;
	if(x >= mod)x -= mod;
	return x;
}
int sub(int x,int y)    {
	x -= y;
	if(x <  0)  x += mod;
	return x;
}
int mul(int x,int y)    {
	return (ll)x * y % mod;
}
int qpow(int a,int exp) {
	int ans = 1;
	while (exp) {
		if(exp & 1) ans = mul(ans,a);
		a = mul(a,a);     exp >>= 1;
	}
	return ans;
}

int root[20][2];
int  rev[20][N];

void fft(int n,vi &a,bool inv)      {
	if (n == 1) return;
	for(int i = 0 ; i < n ; ++i)    {
		int z = __builtin_ctz(n);
		if (i < rev[z][i])
			swap(a[i],a[rev[z][i]]);
	}
	for(int k = 1, e = 1 ; k < n ; k <<= 1, ++e)
		for(int i = 0 ; i < n ; i += (k << 1))  {
			int W = 1;
			for(int j = 0 ; j < k ; ++j)    {
				int x = a[i + j], y = mul(a[i + j + k],W);
				a[i + j]     = add(x,y);
				a[i + j + k] = sub(x,y);
				W = mul(W,root[e][inv]);
			}
		}	
	if (inv)    {
		int m = qpow(n,mod - 2);
		for(int i = 0 ; i < n ; ++i)
			a[i] = mul(a[i],m);
	}
}

vi  Mul(vi  P,vi  Q)    {   //multiple 2 polynomials in complexity O(nlog(n))
	int need = sz(P) + sz(Q) - 1;
	int n = 1;
	while(n < need) n <<= 1;
	P.resize(n);    fft(n,P,0);
	Q.resize(n);    fft(n,Q,0);
	for(int i = 0 ; i < n ; ++i)
		P[i] = mul(P[i],Q[i]);
	fft(n,P,1); P.resize(need);
	return P;
}
vi  Inv(vi  A)  {   //find polynomial B which satisfies Mul(A,B) is 1 in modulo x ^ (sz(A))
	int n = sz(A), m = n + 1 >> 1;
	if (n == 1)
		return {qpow(A[0],mod - 2)};
	vi  B = Inv(vi(A.begin(),A.begin() + m));
	int s = 1;
	while(s < 2 * n)
		s <<= 1;
	A.resize(s);    fft(s,A,0);
	B.resize(s);    fft(s,B,0);
	for(int i = 0 ; i < s ; ++i)
		B[i] = mul(sub(2,mul(A[i],B[i])),B[i]);
	fft(s,B,1); B.resize(n);
	return B;
}
vi  Div(vi  A,vi  B)    {   //find polynomial C which satisfies deg(A - B.C) < deg(B)
	int n = sz(A), m = sz(B);
	if (n < m)  return vi();
	reverse(all(A));
	reverse(all(B));
	B.resize(n - m + 1);
	B = Mul(A,Inv(B));
	B.resize(n - m + 1);
	reverse(all(B));
	return B;
}
vi  Mod(vi  A,vi  B)    {   //find polynomial C which satisfies B divides A - C
	int n = sz(A), m = sz(B);
	if (n < m)  return A;
	vi  Q = Mul(Div(A,B),B);
	Q.resize(m - 1);
	for(int i = 0 ; i < m - 1 ; ++i)
		Q[i] = sub(A[i],Q[i]);
	return Q;
}

int main()  {	
	for(int i = 1 ; i < 19 ; ++i)   {
		int e = (mod >> i);
		root[i][0] = qpow(3,e); e = mod - 1 - e;
		root[i][1] = qpow(3,e);
		for(int j = 0 ; j < (1 << i) ; ++j) {
			int x = j;
			for(int k = 0, l = i - 1 ; k < l ; ++k, --l)
				if((x >> k & 1) != (x >> l & 1))    {
					x ^= (1 << k);
					x ^= (1 << l);
				}
			rev[i][j] = x;
		}
	}
	
	int n = read(), m = read(); 
	
	vi  P(n + 1), a(m);
	
	for(int i = 0 ; i <= n ; ++i)   P[i] = read();  //P is the polynomial
	for(int i = 0 ; i <  m ; ++i)   a[i] = read();  //a is the set of points need to be compute
	
	//the following code is for multipoint evaluation in complexity O(nlog(n)^2)
	
	vi  T[m << 2 | 1];
	
	function<void(int,int,int)>     init = [&](int i,int L,int R)   {
		if (L == R) {
			T[i] = {mod - a[L],1};
			return;
		}
		int M = L + R >> 1;
		init(i << 1    ,L    ,M);
		init(i << 1 | 1,M + 1,R);
		T[i] = Mul(T[i << 1],T[i << 1 | 1]);
	};
	function<void(int,int,int,vi)>  calc = [&](int i,int L,int R,vi A)  {
		if (i > 1)  A = Mod(A,T[i]);
		if (L == R) {
			printf("%d ",A[0]);
			return;
		}
		int M = L + R >> 1;
		calc(i << 1    ,L    ,M,A);
		calc(i << 1 | 1,M + 1,R,A);
	};
	
	init(1,0,m - 1);
	calc(1,0,m - 1,P);
}
